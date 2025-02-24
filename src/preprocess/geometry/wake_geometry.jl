"""
    discretize_wake(
        duct_coordinates,
        centerbody_coordinates,
        rotorzloc, # rotor axial locations
        wake_length,
        npanels,
        dte_minus_cbte;
    )

Calculate wake sheet panel node z-coordinates.

# Arguments
- `duct_coordinates::Matrix{Float}` : Array of input duct coordinates
- `centerbody_coordinates::Matrix{Float}` : Array of input centerbody_coordinates coordinates
- `rotorzloc ::Vector{Float}` : rotor axial locations
- `wake_length::Float` : non-dimensional length of wake to extend beyond aft-most body trailing edge.
- `npanels::Vector{Int}` : A vector of the number of panels between each discrete point.  For example: [number of panels between the rotors; number of panels between the stator and the first trailing edge; number of panels between the trailing edges; number of panels between the last trailing edge and the end of the wake]
- `dte_minus_cbte::Float` : indicator as to whether the duct trailing edge minus the centerbody trailing edge is positive, zero, or negative.
"""
function discretize_wake(
    duct_coordinates,
    centerbody_coordinates,
    rotorzloc,
    wake_length,
    npanels,
    dte_minus_cbte;
    le_bracket=5,
    fitscale=1e2,
)

    # try a root finder to get the true leading edge
    duct_lez, duct_minz = findmin(@view(duct_coordinates[:, 1]))

    duct_lesp = FLOWMath.Akima(
        @view(duct_coordinates[(duct_minz - le_bracket):(duct_minz + le_bracket), 2]),
        @view(duct_coordinates[(duct_minz - le_bracket):(duct_minz + le_bracket), 1]),
    )

    duct_ler = Roots.find_zero(
        x -> FLOWMath.derivative(duct_lesp, x) * fitscale,
        (
            duct_coordinates[duct_minz - le_bracket, 2],
            duct_coordinates[duct_minz + le_bracket, 2],
        ),
        Roots.Brent();
        atol=eps(),
    )
    duct_lez = duct_lesp(duct_ler)

    # trailing edge shouldn't need anything fancy
    duct_tez = duct_coordinates[1, 1]

    # extract hub leading and trailing edge location
    cb_lez = centerbody_coordinates[1, 1]
    cb_tez = centerbody_coordinates[end, 1]

    # calculate duct chord
    duct_chord = max(duct_tez, cb_tez) - min(duct_lez, cb_lez)
    # duct_chord = duct_tez - duct_lez

    # dimensionalize wake_length
    wake_length *= duct_chord

    # ensure rotors are ordered properly
    !issorted(rotorzloc) && sort!(rotorzloc)

    # make sure that all rotors are inside the duct
    for rzl in rotorzloc
        @assert rzl > duct_lez "Rotor is in front of duct leading edge."
        @assert rzl < duct_tez "Rotor is behind duct trailing edge."
        @assert rzl > cb_lez "Rotor is in front of centerbody leading edge."
        @assert rzl < cb_tez "Rotor is behind centerbody trailing edge."
    end

    # combine all discrete locations into one ordered array
    if iszero(dte_minus_cbte)
        zd = vcat(rotorzloc, cb_tez, cb_tez + wake_length)
        @assert length(npanels) == length(rotorzloc) + 1 "Length of vector `npanels` should be one more than the length of vector `rotorzloc` when the duct and centerbody trailing edges align."
    elseif dte_minus_cbte < 0 #duct_tez < cb_tez
        zd = vcat(rotorzloc, duct_tez, cb_tez, cb_tez + wake_length)
        @assert length(npanels) == length(rotorzloc) + 2 "Length of vector `npanels` should be two more than the length of vector `rotorzloc` when the duct and centerbody trailing edges do not align."
    else #dte_minus_cbte > 0 # duct_tez < cb_tez
        zd = vcat(rotorzloc, cb_tez, duct_tez, duct_tez + wake_length)
        @assert length(npanels) == length(rotorzloc) + 2 "Length of vector `npanels` should be two more than the length of vector `rotorzloc` when the duct and centerbody trailing edges align."
    end

    # calculate indices for the start of each discrete section
    sectionid = cumsum(vcat(1, npanels))

    # construct x-discretization
    zwake = similar(zd, sum(npanels) + 1) .= 0.0
    for i in 1:(length(zd) - 1)
        zrange = range(zd[i], zd[i + 1]; length=npanels[i] + 1)
        zwake[sectionid[i]:sectionid[i + 1]] .= zrange
    end

    # get rotor location sectionid
    ridx = findall(x -> x in rotorzloc, zwake)

    # return dimensionalized wake x-coordinates
    return zwake, ridx, [duct_lez duct_ler]
end

"""
    initialize_wake_grid(rp_duct_coordinates, rp_centerbody_coordinates, zwake, rwake)

Initialize the wake grid.

# Arguments:
- `rp_duct_coordinates::Matrix{Float}` : The re-paneled duct coordinates
- `rp_centerbody_coordinates::Matrix{Float}` : The re-paneled centerbody coordinates
- `zwake::Vector{Float}` : The axial positions of the wake sheet panel nodes
- `rwake::Vector{Float}` : The radial positions of the blade elements for the foremost rotor

# Returns:
- `wake_grid::Array{Float,3}` : 3D Array of axial and radial wake_grid points
"""
function initialize_wake_grid(rp_duct_coordinates, rp_centerbody_coordinates, zwake, rwake)
    TF = promote_type(
        eltype(rp_duct_coordinates),
        eltype(rp_centerbody_coordinates),
        eltype(zwake),
        eltype(rwake),
    )

    # number of streamwise wake_grid locations
    nx = length(zwake)

    # number or radial wake_grid locations
    nr = length(rwake)

    # initialize outputs
    wake_grid = zeros(TF, 2, nx, nr)

    return initialize_wake_grid!(
        wake_grid, rp_duct_coordinates, rp_centerbody_coordinates, zwake, rwake
    )
end

"""
    initialize_wake_grid!(
        wake_grid, rp_duct_coordinates, rp_centerbody_coordinates, zwake, rwake
    )

In-place version of `initialize_wake_grid`.
"""
function initialize_wake_grid!(
    wake_grid, rp_duct_coordinates, rp_centerbody_coordinates, zwake, rwake
)

    # get dimensions
    _, nx, nr = size(wake_grid)

    # set x-locations
    wake_grid[1, :, :] .= zwake

    # set first rotor radial locations
    wake_grid[2, 1, :] .= rwake

    # initial x-position of wake
    rzl = zwake[1]

    # centerbody trailing edge
    cb_tez = rp_centerbody_coordinates[1, end]

    # duct trailing edge
    duct_tez = rp_duct_coordinates[1, 1]

    # set inner radial locations
    _, icenterbody_rotor = findmin(x -> abs(x - rzl), view(rp_centerbody_coordinates, 1, :))
    _, icenterbody_te = findmin(x -> abs(x - cb_tez), zwake)

    wake_grid[2, 1:icenterbody_te, 1] .= rp_centerbody_coordinates[2, icenterbody_rotor:end]
    wake_grid[2, (icenterbody_te + 1):end, 1] .= rp_centerbody_coordinates[2, end]

    # set outer radial locations
    _, leidx = findmin(view(rp_duct_coordinates, 1, :))
    _, iduct_rotor = findmin(x -> abs(x - rzl), view(rp_duct_coordinates, 1, 1:leidx))
    _, iduct_te = findmin(x -> abs(x - duct_tez), zwake)
    wake_grid[2, 1:iduct_te, end] .= rp_duct_coordinates[2, iduct_rotor:-1:1]
    wake_grid[2, (iduct_te + 1):end, end] .= rp_duct_coordinates[2, 1]

    # --- Define Inner wake_grid Points --- #

    # get the area of the annulus at the first rotor
    area = pi * (wake_grid[2, 1, end]^2 - wake_grid[2, 1, 1]^2)

    # get the area of each blade element on the first rotor
    section_area = pi * (wake_grid[2, 1, 2:end] .^ 2 .- wake_grid[2, 1, 1:(end - 1)] .^ 2)

    # apply conservation of mass at each wake_grid point
    for ix in 2:nx
        xcenterbody = wake_grid[1, ix, 1]
        rcenterbody = wake_grid[2, ix, 1]

        xduct = wake_grid[1, ix, end]
        rduct = wake_grid[2, ix, end]

        # calculate the area of the annulus at the current streamwise location
        local_area = pi * (rduct .^ 2 .- rcenterbody .^ 2)

        for ir in 2:(nr - 1)

            # define radial wake_grid point based on conservation of mass
            local_section_area = local_area * section_area[ir - 1] / area
            wake_grid[2, ix, ir] = sqrt(
                local_section_area / pi + wake_grid[2, ix, ir - 1]^2
            )

            # linearly interpolate centerbody and duct x-positions
            frac = (wake_grid[2, ix, ir] - rcenterbody) / (rduct - rcenterbody)
            wake_grid[1, ix, ir] = xcenterbody + frac * (xduct - xcenterbody)
        end
    end

    return wake_grid
end

"""
    generate_wake_grid(
        problem_dimensions,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        Rhub1,
        Rtip1,
        tip_gap1,
        zwake;
        grid_solver_options=GridSolverOptions(),
        verbose=false,
        silence_warnings=true,
    )

Initialize and solve for elliptic grid on which wake sheets are defined.

# Arguments
- `problem_dimensions::` : A ProblemDimensions object
- `rp_duct_coordinates::` : repaneled duct coordinates
- `rp_centerbody_coordinates::` : repaneled centerbody coordinates
- `Rhub1::` : Hub radius of first rotor
- `Rtip1::` : Tip radius of first rotor
- `tip_gap1::` : Tip gap of first rotor (MUST BE ZERO for now)
- `zwake::` : axial positions of wake sheet panel nodes

# Keyword Arguments
- `grid_solver_options::GridSolverOptionsType=GridSolverOptions()` : options for solving the elliptic grid.
- `verbose::Bool=false` : flag to print verbose statements
- `silence_warnings::Bool=true` : flag to supress warnings

# Returns
- `wake_grid::Array{Float,3}` : 3D Array of axial and radial wake_grid points after solution of elliptic system.
"""
function generate_wake_grid(
    problem_dimensions,
    rp_duct_coordinates,
    rp_centerbody_coordinates,
    Rhub1,
    Rtip1,
    tip_gap1,
    zwake;
    grid_solver_options=GridSolverOptions(),
    verbose=false,
    silence_warnings=true,
)
    TF = promote_type(
        eltype(rp_duct_coordinates),
        eltype(rp_centerbody_coordinates),
        typeof(Rhub1),
        typeof(Rtip1),
    )

    wake_grid = zeros(TF, 2, problem_dimensions.nwsn, problem_dimensions.nws)

    return generate_wake_grid!(
        wake_grid,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        Rhub1,
        Rtip1,
        tip_gap1,
        zwake;
        grid_solver_options=grid_solver_options,
        verbose=verbose,
        silence_warnings=silence_warnings,
    )
end

"""
    generate_wake_grid!(
        wake_grid,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        Rhub1,
        Rtip1,
        tip_gap1,
        zwake;
        grid_solver_options=grid_solver_options,
        verbose=false,
        silence_warnings=true,
    )

In-place version of `generate_wake_grid`.
"""
function generate_wake_grid!(
    wake_grid,
    rp_duct_coordinates,
    rp_centerbody_coordinates,
    Rhub1,
    Rtip1,
    tip_gap1,
    zwake;
    grid_solver_options=grid_solver_options,
    verbose=false,
    silence_warnings=true,
)

    #rotor panel edges
    rpe = range(Rhub1, Rtip1, size(wake_grid, 3))

    # wake sheet starting radius including dummy sheets for tip gap.
    if tip_gap1 == 0.0
        rwake = rpe
    else
        rwake = [rpe; Rtip1 + tip_gap1]
    end

    # Initialize wake "grid"
    initialize_wake_grid!(
        wake_grid, rp_duct_coordinates, rp_centerbody_coordinates, zwake, rwake
    )

    # - Relax "Grid" - #

    relax_grid!(
        grid_solver_options, wake_grid; verbose=verbose, silence_warnings=silence_warnings
    )

    return wake_grid
end

"""
    relax_grid!(
        grid_solver_options::GridSolverOptionsType,
        wake_grid;
        verbose=false,
        silence_warnings=true,
        tabchar="    ",
        ntab=1,
    )

Relax/Solve initial wake grid according to elliptic system of equations.

# Arguments
- `grid_solver_options::GridSolverOptionsType' : options for elliptic grid solver
- `wake_grid::Array{Float,3}` : Initialized wake grid

# Keyword Arguments
- `verbose=false::' : flag for printing verbose statements
- `silence_warnings=true::' : flag for supressing warnings
- `tabchar::String="    "::' : string to use for tabbing over verbose statements.
- `ntab::Int=1' : number of tabs for printing verbose statements
"""
function relax_grid!(
    grid_solver_options::GridSolverOptions,
    wake_grid;
    verbose=false,
    silence_warnings=true,
    tabchar="    ",
    ntab=1,
)
    if grid_solver_options.precondition
        if verbose
            println(tabchar^ntab * "Preconditioning Elliptic Grid System using SLOR")
        end
        # - Relax grid to allow Newton solve a tractable starting point - #
        relax_grid!(
            wake_grid;
            iteration_limit=grid_solver_options.iteration_limit,
            atol=grid_solver_options.atol,
            converged=grid_solver_options.converged,
            verbose=verbose,
            tabchar="\t",
            ntab=1,
        )
    end

    # - Converge grid with NLsolve - #

    if !grid_solver_options.converged[1]

        # reset convergence flag
        grid_solver_options.converged[1] = false

        if verbose
            println(
                tabchar^ntab *
                "Solving Elliptic Grid System using $(grid_solver_options.algorithm) Method",
            )
        end

        # solve
        solve_elliptic_grid!(
            wake_grid;
            algorithm=grid_solver_options.algorithm,
            atol=grid_solver_options.atol,
            iteration_limit=grid_solver_options.iteration_limit,
            grid_solver_options.converged,
            grid_solver_options.residual_value,
            grid_solver_options.iterations,
            verbose=verbose,
        )
    else
        if verbose
            println(
                tabchar^ntab *
                "Preconditioning Elliptic Grid System converged within final tolerance, skipping non-linear solve.",
            )
        end
    end

    return wake_grid
end

function relax_grid!(
    grid_solver_options::SLORGridSolverOptions,
    wake_grid;
    verbose=false,
    silence_warnings=true,
)
    if verbose
        println("Solving Elliptic Grid System using SLOR")
    end
    relax_grid!(
        wake_grid;
        iteration_limit=grid_solver_options.iteration_limit,
        atol=grid_solver_options.atol,
        converged=grid_solver_options.converged,
        verbose=verbose,
        tabchar="\t",
        ntab=1,
    )
    return wake_grid
end

"""
    generate_wake_panels(wake_grid)

Generate paneling for each wake sheet emanating from the rotor blade elements.

# Arguments:
- `wake_grid::Array{Float,3}` : axial and radial locations of each wake_grid point (after relaxation/solution)

# Returns:
- `wake_vortex_panels::NamedTuple` : A named tuple of panel values describing the wake vortex panels
"""
function generate_wake_panels(wake_grid)

    # extract wake_grid size
    _, nz, nr = size(wake_grid)

    # define wake lines
    wake_lines = [[wake_grid[1, :, ir]'; wake_grid[2, :, ir]'] for ir in 1:nr]

    # generate paneling for each wake line
    wake_panels = generate_panels(wake_lines; isbody=false)

    return wake_panels
end

"""
    generate_wake_panels!(wake_panels, wake_grid)

In-place version of `generate_wake_panels`.
"""
function generate_wake_panels!(wake_panels, wake_grid)
    # extract wake_grid size
    _, nz, nr = size(wake_grid)

    # define wake lines
    wake_lines = [[wake_grid[1, :, ir]'; wake_grid[2, :, ir]'] for ir in 1:nr]

    # generate paneling for each wake line
    return generate_panels!(wake_panels, wake_lines; isbody=false)
end

"""
    get_wake_k(r, nwn)

Calculate geometric constant for use in later calculation of wake panel node strengths.

# Arguments
- `r::Vector{Float}` : Vector of wake panel node radial positions

# Returns
- `K::Vector{Float}` : Vector of geometric constants used in calculation of panel node strengths.
"""
function get_wake_k(r)
    # initialize output
    K = similar(r)

    return get_wake_k!(K, r)
end

"""
    get_wake_k!(K, r)

In-place version of `get_wake_k`.
"""
function get_wake_k!(K, r)
    # Loop through panels
    for (iw, wnr) in enumerate(r)
        # check if panel has zero radius
        if wnr < eps()
            K[iw] = 0.0
        else
            K[iw] = -1.0 ./ (8.0 .* pi^2 .* wnr .^ 2)
        end
    end

    return K
end
