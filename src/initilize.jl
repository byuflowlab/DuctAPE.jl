#=
Functions for Grid Initialization

Author: Judd Mehr,

Process Notes:
1. Define paneling on duct and hub
2. Evaluate for panel control points
3. Assemble abarij, augmented with Kutta condition
4. Solve for gammabar0i
5. Set initial streamfunction grid
6. Evaluate vx0i, vr0i at grid boundaries
7. use vx0i, vr0i to get equipotential lines
8. Set Q1 = 0
9. Relax grid with SLOR (InterativeSolvers) using vx0i, vr0i for Neumann BCs

=#

"""
    generate_grid_points(duct, rotors, grid_options, debug=false)

Get grid boundary and initial interior points.

**Arguments:**
 - `duct::DuctTAPE.Duct` : Duct Object.
 - `rotors::Array{DuctTAPE.Rotor}` : Array of Rotor objects
 - `grid_options::DuctTAPE.GridOptions` : GridOptions object

**Returns:**
 - `x_grid_points::Matrix{Float64,2}` : 2D Array of x grid points
 - `r_grid_points::Matrix{Float64,2}` : 2D Array of r grid points
"""
function generate_grid_points(duct, rotors, grid_options; debug=false)

    # --- RENAME THINGS FOR CONVENIENCE

    # Rename HUB Geometry for convenience
    hubx = duct.hubxcoordinates
    hubr = duct.hubrcoordinates
    hubspline = fm.Akima(hubx, hubr)
    hubLEx = hubx[1]
    hubLEr = hubr[1]
    hubTEx = hubx[end]
    hubTEr = hubr[end]

    # Rename WALL Geometry for convenience
    wallx = duct.wallinnerxcoordinates
    wallr = duct.wallinnerrcoordinates
    wallspline = fm.Akima(wallx, wallr)
    wallLEx = wallx[1]
    wallLEr = wallr[1]
    wallTEx = wallx[end]
    wallTEr = wallr[end]

    # Rename options for convenience
    nr = grid_options.num_radial_stations

    # --- SETUP/FIND BOUNDARY GEOMETRY

    # -- Get foremost boundaries
    if isempty(rotors)
        #find foremost x coordinate
        xfront = duct.LEx - grid_options.inlet_length * duct.chord
        #define radius of disk
        R = wallLEr - hubLEr
        #TODO: add functionality for user explicit definition of radial positions at some point.
        radial_stations_front = range(hubLEr, wallLEr; length=nr)
    else
        #find foremorst x coordinate
        xfront = minimum([rotors[i].xlocation for i in length(rotors)])
        #get annulus radius
        if xfront < wallx[1]
            rtip = wallLEr
        elseif xfront > wallx[end]
            rtip = wallTEr
        else
            rtip = wallspline(xfront)
        end

        if xfront < hubx[1] || xfront > hubx[end]
            rhub = 0.0
        else
            rhub = hubspline(xfront)
        end
        R = rtip - rhub
        #define radial stations for foremost rotor TODO: add functionality for explicit user definition
        radial_stations_front = range(rhub, rtip; length=nr)
    end

    #get radial station step size for later use
    dr = radial_stations_front[2] - radial_stations_front[1]

    # -- Get rear boundary
    #wake starts at duct trailing edge and extends the input wake_length relative to duct chord.
    xrear = duct.TEx + duct.chord * grid_options.wake_length

    # -- Get x-ranges for inner and outer boundaries
    inlet_range = [xfront; duct.LEx]
    duct_range = [max(duct.LEx, xfront); duct.TEx]
    wake_range = [duct.TEx; xrear]

    # --- DEFINE X GRID SPACING PIECES

    # -- Get inlet spacing if applicable
    if inlet_range[1] < inlet_range[2] && grid_options.inlet_length > 0.0 && isempty(rotors)
        #number of inlet x stations
        numxi = round(Int, (inlet_range[2] - inlet_range[1]) / (1.0 * dr))

        #inlet x stations
        xi = range(inlet_range[1], inlet_range[2]; length=numxi)
    else
        numxi = 0
        xi = []
    end

    # -- Get duct spacing

    #number of duct x stations
    nd = round(Int, (duct_range[2] - duct_range[1]) / (1.0 * dr))
    #duct x stations
    xd = range(duct_range[1], duct_range[2]; length=nd)

    # -- Get wake spacing
    xdsize = deepcopy(dr) #length of last element on duct spacing
    xw = [xdsize + duct_range[2]] #make first wake element same size as last duct element
    xwidx = 1 #initialize index
    #TODO: this method doesn't guarentee the actual prescribed wake length, but should get close
    while xw[end] < wake_range[2]
        xdsize *= grid_options.wake_expansion_factor
        push!(xw, xw[xwidx] + xdsize)
        xwidx += 1
    end
    numxw = length(xw)

    ##alternate option for wake spacing: equal spacing
    #numxw = round(Int, (wake_range[2] - wake_range[1]) / dr)
    #xw = range(wake_range[1], wake_range[2]; length=numxw)[2:end]

    # put all the x stations together
    xstations = [xi[1:(end - 1)]; xd; xw[1:end]]
    #get number of x stations
    nx = length(xstations)

    # --- DEFINE R STATION LOCATIONS ALONG INNER/OUTER BOUNDARIES
    # -- Get radial locations of outer boundary
    router = [0.0 for i in 1:nx]
    for i in 1:nx
        if xstations[i] < wallLEx
            router[i] = wallLEr
        elseif xstations[i] < wallTEx
            router[i] = wallspline(xstations[i])
        else
            router[i] = wallTEr
        end
    end

    # -- Get radial locations of inner boundary
    rinner = [0.0 for i in 1:nx]
    for i in 1:nx
        if xstations[i] < hubLEx
            rinner[i] = hubLEr
        elseif xstations[i] < hubTEx
            rinner[i] = hubspline(xstations[i])
        else
            rinner[i] = hubTEr
        end
    end

    # --- INITIALIZE GRID POINTS
    #initialize grid matrix
    x_grid_points = [0.0 for i in 1:length(xstations), j in 1:nr]
    r_grid_points = [0.0 for i in 1:length(xstations), j in 1:nr]

    #first column of radial grid points are radial stations of foremost rotor
    r_grid_points[1, :] = collect(radial_stations_front)

    #first column of x grid points are all the x position of the foremost rotor
    x_grid_points[1, :] .= xfront

    #first and last rows of points from setup above
    x_grid_points[:, 1] = xstations
    x_grid_points[:, end] = xstations
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

    return x_grid_points, r_grid_points, nx, nr
end

"""
    relax_grid(xg, rg)

Relax grid using elliptic grid solver.

**Arguments:**
- `xg::Matrix{Float64}` : Initial x grid points guess
 - `rg::Matrix{Float64}` : Initial r grid points guess

**Returns:**
 - `x_relax_points::Matrix{Float64}` : Relaxed x grid points guess
 - `r_relax_points::Matrix{Float64}` : Relaxed r grid points guess
"""
function relax_grid(xg, rg, nxi, neta; max_iterations=100, tol=1e-9, verbose=false)

    # initialize output and computation arrays
    xr = [xg[i, j] for i in 1:nxi, j in 1:neta]
    rr = [rg[i, j] for i in 1:nxi, j in 1:neta]
    C = [0.0 for i in 1:nxi]
    D = [0.0 for i in 1:2, j in 1:nxi]

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
    eta = [0.0 for i in 1:neta]
    for j in 2:neta
        # DFDC comment for this is roughly: streamline values (axisymmetric streamfunction). eta(j) = (j-1)**2 gives uniform spacing in r
        # not sure why it doesn't match the actual expression though.
        eta[j] = rg[1, j]^2 - rg[1, j - 1]^2 + eta[j - 1]
    end

    # Get the x positions along the x axis (used later for various ratios compared to internal points as they move)
    xi = [xr[i, 1] for i in 1:nxi]
    # Get the r positions along the inner radial boundary
    hubrstations = [rr[i, 1] for i in 1:nxi]

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

            rej = rg[end, j]
            rem1j = rg[end - 1, j]
            rem2j = rg[end - 2, j]

            #x arc distance between j+1 and j-1 radial positions of outlet
            xarc = xr[end, jp1] - xr[end, jm1]
            #r arc distances between j+1 and j-1 radial positions of outlet
            rarc = rg[end, jp1] - rg[end, jm1]
            # normalized x component of arc length (x direction) of current position on outlet
            xhat = xarc / sqrt(xarc^2 + rarc^2)
            # normalized r component of arc length (r direction) of current position on outlet
            rhat = rarc / sqrt(xarc^2 + rarc^2)

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

            #march down xi direction
            for i in 2:(nxi - 1)

                #rename indices for convenience
                im1 = i - 1
                ip1 = i + 1

                #TODO: What are all these and why are they needed for derivatives?
                ravgplus = 0.5 * (rr[i, j] + rr[i, j + 1])
                ravgminus = 0.5 * (rr[i, j] + rr[i, j - 1])
                ravg = 0.5 * (ravgplus + ravgminus)

                dximinus = xi[i] - xi[i - 1]
                dxiplus = xi[i + 1] - xi[i]
                dxiavg = 0.5 * (dximinus + dxiplus)

                detaminus = eta[j] - eta[j - 1]
                detaplus = eta[j + 1] - eta[j]
                detaavg = 0.5 * (detaminus + detaplus)

                # - Calculate 1st Derivatives (non-uniform spacing)
                x_eta = (xr[i, j + 1] - xr[i, j - 1]) / (2.0* detaavg)
                r_eta = (rr[i, j + 1] - rr[i, j - 1]) / (2.0* detaavg)
                x_xi = (xr[i + 1, j] - xr[i - 1, j]) / (2.0*dxiavg)
                r_xi = (rr[i + 1, j] - rr[i - 1, j]) / (2.0*dxiavg)

                #alpha, beta, and gamma as defined in theory doc.
                alpha = x_eta^2 + r_eta^2
                beta = x_eta * x_xi + r_eta * r_xi
                gamma = x_xi^2 + r_xi^2
                # J = r_eta * x_xi - x_eta * r_xi

                # - Calculate 2nd Derivatives

                #Used in calculating second derivatives with non-constant intervals.
                ximinuscoeff = detaminus * detaplus / (dximinus * dxiavg)
                xipluscoeff = detaminus * detaplus / (dxiplus * dxiavg)
                etaminuscoeff = detaplus / detaavg * ravgminus / ravg
                etapluscoeff = detaminus / detaavg * ravgplus / ravg


                x_xixi =
                    (xr[i + 1, j] - xr[i, j]) * xipluscoeff -
                    (xr[i, j] - xr[i - 1, j]) * ximinuscoeff

                x_etaeta =
                    (xr[i, j + 1] - xr[i, j]) * etapluscoeff -
                    (xr[i, j] - xr[i, j - 1]) * etaminuscoeff

                x_xieta =
                    detaminus *
                    detaplus *
                    (
                        xr[i + 1, j + 1] - xr[i - 1, j + 1] - xr[i + 1, j - 1] +
                        xr[i - 1, j - 1]
                    ) / (4.0 * dxiavg * detaavg)

                r_xixi =
                    (rr[i + 1, j] - rr[i, j]) * xipluscoeff -
                    (rr[i, j] - rr[i - 1, j]) * ximinuscoeff

                r_etaeta =
                    (rr[i, j + 1] - rr[i, j]) * etapluscoeff -
                    (rr[i, j] - rr[i, j - 1]) * etaminuscoeff



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
            return xr, rr
        end

        relaxfactor = relaxfactor1

        if dmax < dset1 * dxy
            relaxfactor = relaxfactor2
        end

        if dmax < dset2 * dxy
            relaxfactor = relaxfactor3
        end

        if dmax < dset3 * dxy
            return xr, rr
        end
    end

    return xr, rr
end

"""
    initialize_grid(duct, rotors, grid_options; max_iterations=-1, tol=1e-9)

Initialize grid via zero-thrust, unit freestream solution.

**Arguments:**
- `duct::DuctTAPE.Duct` : Duct Geometry Object
- `rotors::Array{DuctTAPE.Rotor}` : Array of rotor objects
- `grid_options::DuctTAPE.GridOptions` : Grid options object

**Returns:**
 - `grid::DuctTAPE.Grid` : Grid Object
"""
function initialize_grid(duct, rotors, grid_options; max_iterations=-1, tol=1e-9)

    # get initial grid points
    xg, rg, nx, nr = generate_grid_points(duct, rotors, grid_options)

    # relax grid
    xr, rr = relax_grid(xg, rg, nx, nr; max_iterations=max_iterations, tol=tol)

    return Grid(xr, rr, nx, nr)
end
