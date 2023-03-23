#=

Funtions related to generation of the wake "grid"

Authors: Judd Mehr,

=#

# TODO: Use updated version for the wake geometry generation

"""
    initilaize_grid_points(body_geometry, rotor_x_positions, radial_positions, grid_options, debug=false)

Get grid boundary and initial interior points.

**Arguments:**
- `body_geometry::DuctTAPE.body_geometry` : Duct Geometry Object.
- `rotor_x_positions::Vector{Float}` : Vector of x-positions for the rotors
- `radial_positions::Vector{Float}` : Vector of r-positions of the blade elements for the foremost rotor
- `grid_options::DuctTAPE.GridOptions` : GridOptions object

**Returns:**
- `x_grid_points::Matrix{Float64,2}` : 2D Array of x grid points
- `r_grid_points::Matrix{Float64,2}` : 2D Array of r grid points
"""
function initialize_grid_points(
    body_geometry, rotor_x_positions, radial_positions; wake_length=1.0, debug=false
)

    # --- RENAME THINGS FOR CONVENIENCE
    TF = eltype(body_geometry.duct_inner_spline.ydata)

    duct_range = body_geometry.duct_range
    hub_range = body_geometry.hub_range
    hub_spline = body_geometry.hub_spline
    duct_inner_spline = body_geometry.duct_inner_spline

    nrotors = length(rotor_x_positions)
    nr = length(radial_positions)

    #check that rotor/stator are in the correct order
    for i in 1:(nrotors - 1)
        @assert rotor_x_positions[i] < rotor_x_positions[i + 1]
    end

    # --- SETUP/FIND BOUNDARY GEOMETRY

    # -- Get rear boundary
    # calculate duct trailing edge
    duct_te = max(duct_range[2], hub_range[2])
    # calculate duct chord
    duct_chord = duct_te - min(duct_range[1], hub_range[1])

    #wake starts at body_geometry trailing edge and extends the input wake_length relative to body_geometry chord.
    wake_te = duct_te * (1.0 + wake_length)

    # --- DEFINE X GRID SPACING PIECES

    # choose x-step to be same length as radial step length, assumes that radial positions are linearly spaced
    x_step = radial_positions[2] - radial_positions[1]

    # if there are more than 1 rotors, need to account for following rotor locations so that grid aligns with rotors
    if length(rotor_x_positions) > 1

        # find the lengths between rotors using an approximate step size of the blade element spacing
        x_length = zeros(Int, nrotors - 1)
        for i in 1:(nrotors - 1)
            x_length[i] = ceil(
                Int, (rotor_x_positions[i + 1] - rotor_x_positions[i]) / x_step
            )
            if x_length[i] == 1
                x_length[i] += 1
            end
        end

        #need to save the x-index of the rotor positions for later use
        rotoridxs = [1; cumsum(x_length)]

        #get the length from the last rotor to the end of the wake
        last_length = round(Int, (wake_te - rotor_x_positions[end]) / x_step)

        # put all the x stations together
        wake_x_stations = unique(
            reduce(
                vcat,
                [
                    [
                        range(
                            rotor_x_positions[i],
                            rotor_x_positions[i + 1];
                            length=x_length[i],
                        ) for i in nrotors - 1
                    ]
                    range(rotor_x_positions[end], wake_te; length=last_length)
                ],
            ),
        )

    else

        # rotor index is only the first index if there is only one rotor.
        rotoridxs = [1]

        # - If only one rotor,rotor_x_positions- no need for complicated stuff.
        wake_x_stations = range(rotor_x_positions[1], wake_te; step=x_step)
    end

    #get number of x stations
    nx = length(wake_x_stations)

    # --- DEFINE R STATION LOCATIONS ALONG INNER/OUTER BOUNDARIES
    # -- Get radial locations of outer boundary
    router = zeros(TF, nx)
    for i in 1:nx
        if wake_x_stations[i] < duct_range[2]
            router[i] = duct_inner_spline(wake_x_stations[i])
        else
            router[i] = duct_inner_spline(duct_range[2])
        end
    end

    # -- Get radial locations of inner boundary
    rinner = zeros(TF, nx)
    for i in 1:nx
        if wake_x_stations[i] < hub_range[2]
            rinner[i] = hub_spline(wake_x_stations[i])
        else
            rinner[i] = hub_spline(hub_range[2])
        end
    end

    # --- INITIALIZE GRID POINTS
    #initialize grid matrix
    x_grid_points = zeros(TF, nx, nr)
    r_grid_points = zeros(TF, nx, nr)

    #first column of radial grid points are radial stations of foremost rotor
    r_grid_points[1, :] = radial_positions

    #first column of x grid points are all the x position of the foremost rotor
    x_grid_points[1, :] .= rotor_x_positions[1]

    #first and last rows of points from setup above
    x_grid_points[:, 1] = wake_x_stations
    x_grid_points[:, end] = wake_x_stations
    r_grid_points[:, 1] = rinner
    r_grid_points[:, end] = router

    # --- DEFINE INNER GRID POINTS

    # -- Set up grid definition Loop
    #Get area of foremost rotor annulus
    rotorarea = pi * (router[1]^2 - rinner[1]^2)

    #get the annulus area of the radial stations
    drotorarea = [
        pi * (r_grid_points[1, j + 1]^2 - r_grid_points[1, j]^2) for j in 1:(nr - 1)
    ]

    #get areas of subsequent annuli
    localarea = [pi * (r_grid_points[i, end]^2 - r_grid_points[i, 1]^2) for i in 1:nx]

    # March conservation of mass along x-direction to get grid points.
    for i in 2:nx #for each x station after the rotor plane

        #get the x,r coordinates at the x stations
        xwall = x_grid_points[i, end]
        xhub = x_grid_points[i, 1]
        rwall = r_grid_points[i, end]
        rhub = r_grid_points[i, 1]

        # get the area of the annulus at the current x station
        # localarea = pi * (rwall^2 - rhub^2)

        for j in 2:(nr - 1) #for each radial station between the hub and wall

            #get local change in area using conservation of mass
            dlocalarea = (drotorarea[j - 1] * localarea[i]) / rotorarea

            #define radial grid point based on conservation of mass
            r_grid_points[i, j] = sqrt(dlocalarea / pi + r_grid_points[i, j - 1]^2)

            #NOTE that the rest of the loop only does something if xs aren't the same like they are now.
            #calculate ratio for weighted average
            frac = (r_grid_points[i, j] - rhub) / (rwall - rhub)

            # place x grid point based on weighted ratio
            x_grid_points[i, j] = xhub + frac * (xwall - xhub)
        end # for radial stations
    end # for x stations

    return x_grid_points, r_grid_points, nx, nr, rotoridxs
end

"""
    relax_grid!(xg, rg, nxi, neta; max_iterations, tol)

Relax grid using elliptic grid solver.

**Arguments:**
 - `xg::Matrix{Float64}` : Initial x grid points guess
 - `rg::Matrix{Float64}` : Initial r grid points guess
 - `nxi::Int` : number of xi (x) stations in the grid
 - `neta::Int` : number of eta (r) stations in the grid

**Keyword Arguments:**
 - `max_iterations::Int` : maximum number of iterations to run, default=100
 - `tol::Float` : convergence tolerance, default = 1e-9

**Returns:**
 - `x_relax_points::Matrix{Float64}` : Relaxed x grid points
 - `r_relax_points::Matrix{Float64}` : Relaxed r grid points
"""
function relax_grid!(xr, rr, nxi, neta; max_iterations=100, tol=1e-9, verbose=false)
    TF = eltype(rr)

    # initialize output and computation arrays
    # xr = [xg[i, j] for i in 1:nxi, j in 1:neta]
    # rr = [rg[i, j] for i in 1:nxi, j in 1:neta]
    C = zeros(TF, nxi)
    D = zeros(TF, 2, nxi)

    #set up relaxation factors
    #TODO: how are these decided?
    if max_iterations > 0
        relaxfactor1 = 1.0
        relaxfactor2 = 1.1
        relaxfactor3 = 1.4
    else
        relaxfactor1 = 1.0
        relaxfactor2 = 1.0
        relaxfactor3 = 1.0
        max_iterations = 1
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
    # hubrstations = [rr[i, 1] for i in 1:nxi]
    hubrstations = view(rr, 1:nxi, 1)

    # used as scaling for tolerance
    dxy = max(
        xr[1, 1], rr[1, 1], abs(xi[end] - xi[1]), abs(hubrstations[end] - hubrstations[1])
    )

    #next dfdc goes to AXELL function
    #skip over most of the stuff since the intlet and walls don't get relaxed

    for iterate in 1:max_iterations
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
        if dmax < tol * dxy
            if verbose
                println("iterations = ", iterate)
            end
            return nothing
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
            return nothing
        end
    end

    if verbose
        println("iterations = ", iterate)
    end
    return nothing
end

"""
    generate_wake_panels(x_grid_points, r_grid_points, nr; kwargs)

Generate vector of Axisymmetric panel objects for the wake lines emanating from the rotor blade elements.

**Arguments:**
- `x_grid_points::Matrix{Float}` : x-location of each grid point
- `r_grid_points::Matrix{Float}` : r-location of each grid point
- `nr::Int` : number of radial grid points

**Keyword Arguments:**
- `method::FLOWFoil.AxisymmetricProblem` : default = AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true]),

**Returns:**
- `wake_panels::Vector{FLOWFoil.AxisymmetricPanel}` : vector of panel objects describing the wake lines
"""
function generate_wake_panels(
    x_grid_points,
    r_grid_points,
    nr;
    method=ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true]),
)

    # - Organize Coordinates into wake lines - #
    wake_line_coordinates = [[x_grid_points[:, i] r_grid_points[:, i]] for i in 1:nr]

    # - Use FLOWFoil to Generate Vortex Panels - #
    wake_panels = ff.generate_panels(method, wake_line_coordinates)

    return wake_panels
end

"""

    generate_wake_panels(x_grid_points, r_grid_points, nr; kwargs)

Generate vector of Axisymmetric panel objects for the wake lines emanating from the rotor blade elements.

**Arguments:**
- `body_geometry::BodyGeometry` : BodyGeometry object describing the duct and hub
- `rotor_x_positions::Vector{Float}` : Vector of x-positions for the rotors
- `radial_positions::Vector{Float}` : Vector of r-positions of the blade elements for the foremost rotor

**Keyword Arguments:**
- `wake_length::Float` : length of wake (non-dimensional based on maximum duct chord) to extend past the furthest trailing edge.
- `method::FLOWFoil.AxisymmetricProblem` : default = AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true]),

**Returns:**
- `wake_panels::Vector{FLOWFoil.AxisymmetricPanel}` : vector of panel objects describing the wake lines

"""
function generate_wake_grid(
    body_geometry,
    rotor_x_positions,
    radial_positions;
    wake_length=1.0,
    method=ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true]),
    debug=false,
)

    # - Initialize Grid Points - #
    x_grid_points, r_grid_points, nx, nr, rotoridxs = initialize_grid_points(
        body_geometry, rotor_x_positions, radial_positions; wake_length=1.0, debug=false
    )

    # - Relax Grid Points - #
    relax_grid!(x_grid_points, r_grid_points, nx, nr)

    # - Generate Panels - #
    # These are the grid centers
    # TODO: need to determine where the wake panels actually lie. probably not like this.
    wake_panels = generate_wake_panels(
        (x_grid_points[:, 1:(end - 1)] .+ x_grid_points[:, 2:end]) / 2.0,
        (r_grid_points[:, 1:(end - 1)] .+ r_grid_points[:, 2:end]) / 2.0,
        nr - 1;
        method=method,
    )

    # - Put the points together in one grid - #
    #TODO: not sure if this is going to be used elsewhere, may remove this later
    # wake_grid = [[x_grid_points[j, i] r_grid_points[j, i]] for i in 1:nr, j in 1:nx]

    return wake_panels, length(wake_panels), rotoridxs, r_grid_points
end
