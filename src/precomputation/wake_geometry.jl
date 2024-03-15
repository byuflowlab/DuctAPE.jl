"""
    discretize_wake(duct_coordinates, centerbody_coordinates, rotor_parameters, wake_length, nwake_sheets, dte_minus_cbte)

Calculate wake x-coordinates.
npanels is a vector of number of panels between each discrete point.  so something like [number of panels between the rotors; number of panels between the stator and the first trailing edge; number of panels between the trailing edges; number of panels between the last trailing edge and the end of the wake]
"""
function discretize_wake(
    duct_coordinates,
    centerbody_coordinates,
    rotorzloc, # rotor axial locations
    wake_length,
    npanels,
    dte_minus_cbte;
)

    # extract duct leading and trailing edge locations
    duct_lez = minimum(duct_coordinates[:, 1])
    duct_tez = maximum(duct_coordinates[:, 1])

    # extract hub leading and trailing edge location
    cb_lez = minimum(centerbody_coordinates[:, 1])
    cb_tez = maximum(centerbody_coordinates[:, 1])

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

    # # get sectionid of body trailing edges on wake
    # ductTE_index = searchsortedfirst(zwake, duct_tez)
    # centerbodyTE_index = searchsortedfirst(zwake, cb_tez)

    # return dimensionalized wake x-coordinates
    return zwake, ridx#, ductTE_index, centerbodyTE_index
end

"""
    initialize_wake_grid(body_geometry, rzl, rblade)

Intialize the wake wake_grid

# Arguments:
- `body_geometry::DuctAPE.body_geometry` : Duct Geometry Object.
- `rzl::Vector{Float}` : Vector of x-positions for the rotors
- `rwake::Vector{Float}` : Vector of r-positions of the blade elements for the foremost rotor
- `grid_options::DuctAPE.GridOptions` : GridOptions object

# Returns:
- `zgrid::Matrix{Float64,2}` : 2D Array of x wake_grid points
- `rgrid::Matrix{Float64,2}` : 2D Array of r wake_grid points
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

function generate_wake_grid(
    problem_dimensions,
    rp_duct_coordinates,
    rp_centerbody_coordinates,
    Rhub1,
    Rtip1,
    tip_gap1,
    zwake;
    wake_nlsolve_ftol=1e-14,
    wake_max_iter=100,
    max_wake_relax_iter=3,
    wake_relax_tol=1e-14,
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
        wake_nlsolve_ftol=wake_nlsolve_ftol,
        wake_max_iter=wake_max_iter,
        max_wake_relax_iter=max_wake_relax_iter,
        wake_relax_tol=wake_relax_tol,
        verbose=verbose,
        silence_warnings=silence_warnings,
    )
end

function generate_wake_grid!(
    wake_grid,
    rp_duct_coordinates,
    rp_centerbody_coordinates,
    Rhub1,
    Rtip1,
    tip_gap1,
    zwake;
    wake_nlsolve_ftol=1e-14,
    wake_max_iter=100,
    max_wake_relax_iter=3,
    wake_relax_tol=1e-14,
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

    # Relax "Grid"
    relax_grid!(
        wake_grid;
        max_wake_relax_iter=max_wake_relax_iter,
        wake_relax_tol=wake_relax_tol,
        verbose=verbose,
        silence_warnings=silence_warnings,
    )
    solve_elliptic_grid!(
        wake_grid;
        wake_nlsolve_ftol=wake_nlsolve_ftol,
        wake_max_iter=wake_max_iter,
        verbose=verbose,
    )

    return wake_grid
end

"""
    generate_wake_panels(zgrid, rgrid; kwargs...)

Generate paneling for each wake line emanating from the rotor blade elements.

# Arguments:
- `zgrid::Matrix{Float}` : x-location of each wake_grid point
- `rgrid::Matrix{Float}` : r-location of each wake_grid point

# Keyword Arguments:
- `method::FLOWFoil.AxisymmetricProblem` : default = AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true]),

# Returns:
- `wake_panels::Vector{FLOWFoil.AxisymmetricPanel}` : vector of panel objects describing the wake lines
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
"""
function get_wake_k(wake_vortex_panels)
    # initialize output
    K = zeros(eltype(wake_vortex_panels.node), wake_vortex_panels.totnode)

    # Loop through panels
    for (iw, wnr) in enumerate(wake_vortex_panels.node[2, :])
        # check if panel has zero radius
        if wnr < eps()
            K[iw] = 0.0
        else
            K[iw] = -1.0 ./ (8.0 .* pi^2 .* wnr .^ 2)
        end
    end

    return K
end

#---------------------------------#
#    DFDC-like wake relaxation    #
#---------------------------------#
"""
    relax_grid!(xg, rg, nxi, neta; max_wake_relax_iter, wake_relax_tol)

Relax wake_grid using elliptic wake_grid solver.

# Arguments:
 - `xg::Matrix{Float64}` : Initial x wake_grid points guess
 - `rg::Matrix{Float64}` : Initial r wake_grid points guess
 - `nxi::Int` : number of xi (x) stations in the wake_grid
 - `neta::Int` : number of eta (r) stations in the wake_grid

# Keyword Arguments:
 - `max_wake_relax_iter::Int` : maximum number of iterations to run, default=100
 - `wake_relax_tol::Float` : convergence tolerance, default = 1e-9

# Returns:
 - `x_relax_points::Matrix{Float64}` : Relaxed x wake_grid points
 - `r_relax_points::Matrix{Float64}` : Relaxed r wake_grid points
"""
function relax_grid!(
    wake_grid;
    max_wake_relax_iter=100,
    wake_relax_tol=1e-9,
    verbose=false,
    silence_warnings=true,
)
    xr = view(wake_grid, 1, :, :)
    rr = view(wake_grid, 2, :, :)

    TF = eltype(wake_grid)
    _, nxi, neta = size(wake_grid)

    # initialize output and computation arrays
    # xr = [xg[i, j] for i in 1:nxi, j in 1:neta]
    # rr = [rg[i, j] for i in 1:nxi, j in 1:neta]
    C = zeros(TF, nxi)
    D = zeros(TF, 2, nxi)

    #set up relaxation factors
    #TODO: how are these decided?
    if max_wake_relax_iter > 0
        relaxfactor1 = 1.0
        relaxfactor2 = 1.1
        relaxfactor3 = 1.4
    else
        relaxfactor1 = 1.0
        relaxfactor2 = 1.0
        relaxfactor3 = 1.0
        max_wake_relax_iter = 1
    end

    relaxfactor = relaxfactor1

    # convergence tolerances for each phase
    # #TODO: where do these come from?
    dset1 = 1e-1
    dset2 = 5e-3
    dset3 = 5e-7

    #initialize streamline values (used to calculate derivatives later)
    eta = zeros(TF, neta)
    for j in 2:neta
        # DFDC comment for this is roughly: streamline values (axisymmetric streamfunction). eta(j) = (j-1)**2 gives uniform spacing in r
        # not sure why it doesn't match the actual expression though.
        eta[j] = view(rr, 1, j) .^ 2 - view(rr, 1, j - 1) .^ 2 + eta[j - 1]
    end

    # Get the x positions along the x axis (used later for various ratios compared to internal points as they move)
    xi = [xr[i, 1] for i in 1:nxi]
    # Get the r positions along the inner radial boundary
    # centerbodyrstations = [rr[i, 1] for i in 1:nxi]
    centerbodyrstations = view(rr, 1:nxi, 1)

    # used as scaling for tolerance
    dxy = max(
        xr[1, 1],
        rr[1, 1],
        abs(xi[end] - xi[1]),
        abs(centerbodyrstations[end] - centerbodyrstations[1]),
    )

    #next dfdc goes to AXELL function
    #skip over most of the stuff since the intlet and walls don't get relaxed

    for iterate in 1:max_wake_relax_iter
        if verbose
            println("iteration $iterate")
        end

        #initialize max step
        dmax = 0.0

        #loop over interior streamlines:
        for j in 2:(neta - 1)
            #rename indices for convenience
            jm1 = j - 1
            jp1 = j + 1

            ## -- Relax Outlet -- ##

            #rename outlet points for convenience
            xej = xr[end, j]
            xem1j = xr[end - 1, j]
            xem2j = xr[end - 2, j]

            rej = rr[end, j]
            rem1j = rr[end - 1, j]
            rem2j = rr[end - 2, j]

            #x arc distance between j+1 and j-1 radial positions of outlet
            xarc = xr[end, jp1] - xr[end, jm1]
            #r arc distances between j+1 and j-1 radial positions of outlet
            rarc = rr[end, jp1] - rr[end, jm1]

            # Check if sqrt of zero will take place
            if isapprox(xarc, rarc)
                xhat = TF(0.0)
                rhat = TF(0.0)
            else
                # normalized x component of arc length (x direction) of current position on outlet
                xhat = xarc / sqrt(xarc^2 + rarc^2)

                # normalized r component of arc length (r direction) of current position on outlet
                rhat = rarc / sqrt(xarc^2 + rarc^2)
            end

            #xi distance between outlet and second to last xi station before outlet
            dxem1 = xi[end] - xi[end - 1]
            #xi distance along inner bound between second and third to last xi stations before outlet.
            dxem2 = xi[end - 1] - xi[end - 2]

            #ratios of x and r spacing comparing current interior r position vs inner bound x position
            dxem1dxinb = (xej - xem1j) / dxem1
            dxem2dxinb = (xem1j - xem2j) / dxem2
            drem1dxinb = (rej - rem1j) / dxem1
            drem2dxinb = (rem1j - rem2j) / dxem2

            #2nd-order 3-point difference to get tangential velocity
            #TODO: figure out what is happening here.
            rez =
                xhat * (dxem1dxinb + dxem1 * (dxem1dxinb - dxem2dxinb) / (dxem1 + dxem2)) +
                rhat * (drem1dxinb + dxem1 * (drem1dxinb - drem2dxinb) / (dxem1 + dxem2))

            x_xij = xhat * (1.0 / dxem1 + 1.0 / (dxem1 + dxem2))
            x_rij = rhat * (1.0 / dxem1 + 1.0 / (dxem1 + dxem2))
            x_xi = x_xij * xhat + x_rij * rhat

            # outletrelax = 1.0 * relaxfactor
            doutlet = -relaxfactor * rez / x_xi

            xr[end, j] += doutlet * xhat
            rr[end, j] += doutlet * rhat

            ## -- END Outlet Relaxation -- ##

            ## -- Relax points on the jth streamline -- ##
            #TODO:NEED TO RECONCILE WRITTEN THEORY WITH CODE...
            #TODO: can all the derivatives be replaced with FLOWMath or similar?

            #march down xi direction
            for i in 2:(nxi - 1)

                #rename indices for convenience
                im1 = i - 1
                ip1 = i + 1

                #part of the first derivative calculations
                dximinus = xi[i] - xi[i - 1]
                dxiplus = xi[i + 1] - xi[i]
                dxiavg = 0.5 * (dximinus + dxiplus)

                #part of the first derivative calculations
                detaminus = eta[j] - eta[j - 1]
                detaplus = eta[j + 1] - eta[j]
                detaavg = 0.5 * (detaminus + detaplus)

                # - Calculate 1st Derivatives (non-uniform spacing)
                x_eta = (xr[i, j + 1] - xr[i, j - 1]) / (2.0 * detaavg)
                r_eta = (rr[i, j + 1] - rr[i, j - 1]) / (2.0 * detaavg)
                x_xi = (xr[i + 1, j] - xr[i - 1, j]) / (2.0 * dxiavg)
                r_xi = (rr[i + 1, j] - rr[i - 1, j]) / (2.0 * dxiavg)

                #alpha, beta, and gamma as defined in theory doc.
                alpha = x_eta^2 + r_eta^2
                beta = x_eta * x_xi + r_eta * r_xi
                gamma = x_xi^2 + r_xi^2
                # J = r_eta * x_xi - x_eta * r_xi

                # - Calculate 2nd Derivatives

                #Used in calculating second derivatives with non-constant intervals.
                ravgplus = 0.5 * (rr[i, j] + rr[i, j + 1])
                ravgminus = 0.5 * (rr[i, j] + rr[i, j - 1])
                ravg = 0.5 * (ravgplus + ravgminus)

                ximinuscoeff = detaminus * detaplus / (dximinus * dxiavg)
                xipluscoeff = detaminus * detaplus / (dxiplus * dxiavg)
                etaminuscoeff = detaplus / detaavg * ravgminus / ravg
                etapluscoeff = detaminus / detaavg * ravgplus / ravg

                #x_ξξ
                x_xixi =
                    (xr[i + 1, j] - xr[i, j]) * xipluscoeff -
                    (xr[i, j] - xr[i - 1, j]) * ximinuscoeff

                #x_ηη
                x_etaeta =
                    (xr[i, j + 1] - xr[i, j]) * etapluscoeff -
                    (xr[i, j] - xr[i, j - 1]) * etaminuscoeff

                #x_ξη
                x_xieta =
                    detaminus *
                    detaplus *
                    (
                        xr[i + 1, j + 1] - xr[i - 1, j + 1] - xr[i + 1, j - 1] +
                        xr[i - 1, j - 1]
                    ) / (4.0 * dxiavg * detaavg)

                #r_ξξ
                r_xixi =
                    (rr[i + 1, j] - rr[i, j]) * xipluscoeff -
                    (rr[i, j] - rr[i - 1, j]) * ximinuscoeff

                #r_ηη
                r_etaeta =
                    (rr[i, j + 1] - rr[i, j]) * etapluscoeff -
                    (rr[i, j] - rr[i, j - 1]) * etaminuscoeff

                #r_ξη
                r_xieta =
                    detaminus *
                    detaplus *
                    (
                        rr[i + 1, j + 1] - rr[i - 1, j + 1] - rr[i + 1, j - 1] +
                        rr[i - 1, j - 1]
                    ) / (4.0 * dxiavg * detaavg)

                #D's look to be from eqns 92 and 93 in theory doc, but don't actually match up completely.  They don't even match up with the notes in the code...
                D[1, i] =
                    alpha * x_xixi - 2.0 * beta * x_xieta + gamma * x_etaeta -
                    beta * r_xi * x_eta * detaminus * detaplus / ravg

                D[2, i] =
                    alpha * r_xixi - 2.0 * beta * r_xieta + gamma * r_etaeta -
                    beta * r_xi * r_eta * detaminus * detaplus / ravg

                #SLOR Stuff
                #TODO: replace with linearsolve package?
                #TODO: what are A, B, and C?
                A =
                    alpha * (ximinuscoeff + xipluscoeff) +
                    gamma * (etaminuscoeff + etapluscoeff)

                if i == 2
                    B = 0 #TODO: Why?
                else
                    B = -alpha * ximinuscoeff
                end

                C[i] = -alpha * xipluscoeff

                ainv = 1.0 / (A - B * C[i - 1])
                C[i] *= ainv
                D[1, i] = (D[1, i] - B * D[1, i - 1]) * ainv
                D[2, i] = (D[2, i] - B * D[2, i - 1]) * ainv
            end #for i (streamwise stations)

            # Rest of SLOR Stuff
            # run through streamline backwards?
            for ib in 2:(nxi - 1)
                i = nxi - ib + 1

                D[1, i] -= C[i] * D[1, i + 1]
                D[2, i] -= C[i] * D[2, i + 1]

                #over relaxation step?
                xr[i, j] += relaxfactor * D[1, i]
                rr[i, j] += relaxfactor * D[2, i]

                # stuff for convergence checking
                ad1 = abs(D[1, i])
                ad2 = abs(D[2, i])

                dmax = max(dmax, ad1, ad2)
            end #for ib
        end #for j (radial stations)

        # -- Update relaxation factors
        #TODO: need to figure out how the dset numbers and relaxation factors are chosen.  Why are they the values that they are? Does it have something to do with the SLOR setup?
        if dmax < wake_relax_tol * dxy
            if verbose
                println("iterations = ", iterate)
            end
            return wake_grid
        end

        relaxfactor = relaxfactor1

        if dmax < dset1 * dxy
            relaxfactor = relaxfactor2
        end

        if dmax < dset2 * dxy
            relaxfactor = relaxfactor3
        end

        if dmax < dset3 * dxy
            if verbose
                println("iterations = ", iterate)
            end
            return wake_grid
        end
    end

    if verbose
        println("iterations = ", iterate)
    end

    if !silence_warnings
        @warn "Wake grid relaxation did not converge, iteration limit of $(max_wake_relax_iter) met."
    end

    return wake_grid
end
