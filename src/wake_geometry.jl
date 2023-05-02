"""
    discretize_wake(duct_coordinates, hub_coordinates, rotor_parameters, wake_length, nwake_sheets)

Calculate wake x-coordinates.
npanels is a vector of number of panels between each discrete point.  so something like [number of panels between the rotors; number of panels between the stator and the first trailing edge; number of panels between the trailing edges; number of panels between the last trailing edge and the end of the wake]
"""
function discretize_wake(
    duct_coordinates,
    hub_coordinates,
    xrotors, # rotor x locations
    wake_length,
    npanels;
)

    # extract duct leading and trailing edge locations
    xduct_le = minimum(duct_coordinates[:, 1])
    xduct_te = maximum(duct_coordinates[:, 1])

    # extract hub leading and trailing edge location
    xhub_le = minimum(hub_coordinates[:, 1])
    xhub_te = maximum(hub_coordinates[:, 1])

    # calculate duct chord
    duct_chord = max(xduct_te, xhub_te) - min(xduct_le, xhub_le)

    # dimensionalize wake_length
    wake_length *= duct_chord

    # ensure rotors are ordered properly
    @assert issorted(xrotors)

    # make sure that all rotors are inside the duct
    for xrotor in xrotors
        @assert xrotor > xduct_le "Rotor is in front of duct leading edge."
        @assert xrotor < xduct_te "Rotor is behind duct trailing edge."
        @assert xrotor > xhub_le "Rotor is in front of hub leading edge."
        @assert xrotor < xhub_te "Rotor is behind hub trailing edge."
    end

    # combine all discrete locations into one ordered array
    if isapprox(xhub_te, xduct_te)
        xd = vcat(xrotors, xhub_te, xhub_te + wake_length)
        @assert length(npanels) == length(xrotors) + 1 "Length of vector `npanels` should be one more than the length of vector `xrotors` when the duct and hub trailing edges align."
    elseif xduct_te < xhub_te
        xd = vcat(xrotors, xduct_te, xhub_te, xhub_te + wake_length)
        @assert length(npanels) == length(xrotors) + 2 "Length of vector `npanels` should be two more than the length of vector `xrotors` when the duct and hub trailing edges align."
    else
        xd = vcat(xrotors, xhub_te, xduct_te, xduct_te + wake_length)
        @assert length(npanels) == length(xrotors) + 2 "Length of vector `npanels` should be two more than the length of vector `xrotors` when the duct and hub trailing edges align."
    end

    # # find (approximate) hub and tip locations
    # _, leidx = findmin(view(duct_coordinates, :, 1))
    # _, ihub = findmin(x -> abs(x - xrotors[1]), view(hub_coordinates, :, 1))
    # _, iduct = findmin(x -> abs(x - xrotors[1]), view(duct_coordinates, 1:leidx, 1))
    # Rhub = hub_coordinates[ihub, 2]
    # Rtip = duct_coordinates[iduct, 2]

    # calculate number of panels per unit length
    # panel_density = wake_refinement * nwake_sheets / (Rtip - Rhub)

    # calculate number of panels for each discrete section
    # npanels = [ceil(Int, (xd[i + 1] - xd[i]) * panel_density) for i in 1:(length(xd) - 1)]

    # calculate indices for the start of each discrete section
    indices = cumsum(vcat(1, npanels))

    # construct x-discretization
    xwake = similar(xd, sum(npanels) + 1)
    for i in 1:(length(xd) - 1)
        xrange = range(xd[i], xd[i + 1]; length=npanels[i] + 1)
        xwake[indices[i]:indices[i + 1]] .= xrange
    end

    # get rotor location indices
    ridx = findall(x -> x in xrotors, xwake)
    ductTE_index = findall(x -> x == xduct_te, xwake)
    hubTE_index = findall(x -> x == xhub_te, xwake)

    # return dimensionalized wake x-coordinates
    return xwake, ridx, ductTE_index, hubTE_index
end

"""

    generate_wake_grid(body_geometry, xrotor, rblade, wake_length=1.0; kwargs...)

Generate the wake grid.

# Arguments:
- `body_geometry::BodyGeometry`: Duct and hub geometry
- `xrotor::Vector{Float}`: x-position of each rotor
- `rblade::Vector{Float}` : r-position of each blade element (for the first rotor)
- `wake_length::Float` : non-dimensional length (based on maximum duct chord) that the wake
    extends past the furthest trailing edge.

# Keyword Arguments:
- `method = AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true])`

"""
function generate_wake_grid(xwake, rwake, wake_length=1.0)

    # initialize the wake grid
    xgrid, rgrid = initialize_wake_grid(body_geometry, xrotor, rblade, wake_length)

    # align the wake grid with streamlines
    update_grid!(xgrid, rgrid)

    return xgrid, rgrid, rotor_indices
end

"""
    generate_wake_panels(xgrid, rgrid; kwargs...)

Generate paneling for each wake line emanating from the rotor blade elements.

# Arguments:
- `xgrid::Matrix{Float}` : x-location of each grid point
- `rgrid::Matrix{Float}` : r-location of each grid point

**Keyword Arguments:**
- `method::FLOWFoil.AxisymmetricProblem` : default = AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true]),

**Returns:**
- `wake_panels::Vector{FLOWFoil.AxisymmetricPanel}` : vector of panel objects describing the wake lines
"""
function generate_wake_panels(
    xgrid, rgrid; method=ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true])
)
    @assert size(xgrid) == size(rgrid)

    # extract grid size
    nx, nr = size(xgrid)

    # define wake lines
    wake_lines = [[xgrid[:, ir] rgrid[:, ir]] for ir in 1:nr]

    # generate paneling for each wake line
    wake_panels = ff.generate_panels(method, wake_lines)

    return wake_panels
end

"""
    initialize_wake_grid(body_geometry, xrotor, rblade)

Intialize the wake grid

**Arguments:**
- `body_geometry::DuctTAPE.body_geometry` : Duct Geometry Object.
- `xrotor::Vector{Float}` : Vector of x-positions for the rotors
- `rwake::Vector{Float}` : Vector of r-positions of the blade elements for the foremost rotor
- `grid_options::DuctTAPE.GridOptions` : GridOptions object

**Returns:**
- `xgrid::Matrix{Float64,2}` : 2D Array of x grid points
- `rgrid::Matrix{Float64,2}` : 2D Array of r grid points
"""
function initialize_wake_grid(duct_coordinates, hub_coordinates, xwake, rwake)

    # number of streamwise grid locations
    nx = length(xwake)

    # number or radial grid locations
    nr = length(rwake)

    # initialize outputs
    xgrid = similar(duct_coordinates, nx, nr)
    rgrid = similar(duct_coordinates, nx, nr)

    # set x-locations
    xgrid .= xwake

    # set first rotor radial locations
    rgrid[1, :] .= rwake

    # initial x-position of wake
    xrotor = xwake[1]

    # hub trailing edge
    xhub_te = hub_coordinates[end, 1]

    # duct trailing edge
    xduct_te = duct_coordinates[1, 1]

    # set inner radial locations
    _, ihub_rotor = findmin(x -> abs(x - xrotor), view(hub_coordinates, :, 1))
    _, ihub_te = findmin(x -> abs(x - xhub_te), xwake)

    rgrid[1:ihub_te, 1] .= hub_coordinates[ihub_rotor:end, 2]
    rgrid[(ihub_te + 1):end, 1] .= hub_coordinates[end, 2]

    # set outer radial locations
    _, leidx = findmin(view(duct_coordinates, :, 1))
    _, iduct_rotor = findmin(x -> abs(x - xrotor), view(duct_coordinates, 1:leidx, 1))
    _, iduct_te = findmin(x -> abs(x - xduct_te), xwake)
    rgrid[1:iduct_te, end] .= duct_coordinates[iduct_rotor:-1:1, 2]
    rgrid[(iduct_te + 1):end, end] .= duct_coordinates[1, 2]

    # --- Define Inner Grid Points --- #

    # get the area of the annulus at the first rotor
    area = pi * (rgrid[1, end]^2 - rgrid[1, 1]^2)

    # get the area of each blade element on the first rotor
    section_area = pi * (rgrid[1, 2:end] .^ 2 .- rgrid[1, 1:(end - 1)] .^ 2)

    # apply conservation of mass at each grid point
    for ix in 2:nx
        xhub = xgrid[ix, 1]
        rhub = rgrid[ix, 1]

        xduct = xgrid[ix, end]
        rduct = rgrid[ix, end]

        # calculate the area of the annulus at the current streamwise location
        local_area = pi * (rduct .^ 2 .- rhub .^ 2)

        for ir in 2:(nr - 1)

            # define radial grid point based on conservation of mass
            local_section_area = local_area * section_area[ir - 1] / area
            rgrid[ix, ir] = sqrt(local_section_area / pi + rgrid[ix, ir - 1]^2)

            # linearly interpolate hub and duct x-positions
            frac = (rgrid[ix, ir] - rhub) / (rduct - rhub)
            xgrid[ix, ir] = xhub + frac * (xduct - xhub)
        end
    end

    return xgrid, rgrid
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
function relax_grid!(xr, rr; max_iterations=100, tol=1e-9, verbose=false)
    TF = eltype(rr)
    nxi, neta = size(xr)

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
