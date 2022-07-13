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

    # --- Define Grid Edge Locations

    # - Get forward x edge location of forward most rotor plane
    #TODO: need to figure out how to define rotor plane (quarter chord after rotation?)
    # if duct only, include inlet
    if isempty(rotors)
        xfront = duct.LEx - grid_options.inlet_length * duct.chord

        #else front of grid is foremost rotor plane
    else
        xfront = minimum([rotors[i].xlocation for i in length(rotors)])
    end

    # - Get rear x edge on wake_length input and duct chord
    #wake starts at duct trailing edge and extends the input wake_length relative to duct chord.
    xrear = duct.TEx + duct.chord * grid_options.wake_length

    # - Get inner r edges using cobination of hub geometry and TE.
    hubx = duct.hubxcoordinates
    hubr = duct.hubrcoordinates
    hubspline = fm.Akima(hubx, hubr)
    hubLEx = hubx[1]
    hubLEr = hubr[1]
    hubTEx = hubx[end]
    hubTEr = hubr[end]

    # - Get outer r edges using combination of wall geometry and TE
    wallx = duct.wallinnerxcoordinates
    wallr = duct.wallinnerrcoordinates
    wallspline = fm.Akima(wallx, wallr)
    wallLEx = wallx[1]
    wallLEr = wallr[1]
    wallTEx = wallx[end]
    wallTEr = wallr[end]

    # --- Get x-ranges
    inlet_range = [xfront; duct.LEx]
    duct_range = [duct.LEx; duct.TEx]
    wake_range = [duct.TEx; xrear]

    # --- Define Grid Spacings
    #TODO: start with simple linear spacing, add refinement capabilities later.
    #radial spacing
    nr = grid_options.num_radial_stations

    # inlet spacing
    numxi = grid_options.num_xinlet_stations
    if inlet_range[1] < inlet_range[2] && grid_options.inlet_length > 0.0
        xi = range(inlet_range[1], inlet_range[2]; length=numxi)
    else
        xi = []
    end

    #duct spacing
    numxd = grid_options.num_xduct_stations
    xd = range(duct_range[1], duct_range[2]; length=numxd)

    # wake spacing
    xdsize = xd[end] - xd[end - 1] #length of last element on duct spacing
    xw = [xdsize + wake_range[1]] #make first wake element same size as last duct element
    xwidx = 1 #initialize index
    #TODO: this method doesn't guarentee the actual prescribed wake length, but should get close
    while xw[end] < wake_range[2]
        xdsize *= grid_options.wake_expansion_factor
        push!(xw, xw[xwidx] + xdsize)
        xwidx += 1
    end
    numxw = length(xw)

    #alternate option: equal spacing
    # xw = range(wake_range[1], wake_range[2]; length=numxw)

    # put all the x stations together
    numxall = numxi + numxd + numxw
    xall = [xi; xd; xw]

    # --- Define Grid Points

    # -- Get radial locations of outer boundary
    router = [0.0 for i in 1:numxall]
    for i in 1:numxall
        if xall[i] < wallLEx
            router[i] = wallLEr
        elseif xall[i] < wallTEx
            router[i] = wallspline(xall[i])
        else
            router[i] = wallTEr
        end
    end

    # -- Get radial locations of inner boundary
    rinner = [0.0 for i in 1:numxall]
    for i in 1:numxall
        if xall[i] < hubLEx
            rinner[i] = hubLEr
        elseif xall[i] < hubTEx
            rinner[i] = hubspline(xall[i])
        else
            rinner[i] = hubTEr
        end
    end

    # --- Define Grid Points
    #initialize grid matrix
    x_grid_points = [0.0 for i in 1:length(xall), j in 1:nr]
    r_grid_points = [0.0 for i in 1:length(xall), j in 1:nr]

    #Start with points on rotor plane from hub to wall
    x_grid_points[1, :] .= xfront
    r_grid_points[1, :] = collect(range(hubr[1], wallr[1]; length=nr))

    x_grid_points[:, 1] = xall
    x_grid_points[:, end] = xall
    r_grid_points[:, 1] = rinner
    r_grid_points[:, end] = router

    # March conservation of mass along x-direction to get grid points.
    #get the annulus area at the rotor plane
    rotorarea = pi * (wallr[1]^2 - hubr[1]^2)
    for i in 2:length(xall) #for each x station after the rotor plane

        #get the x,r coordinates at the x stations
        xwall = xall[i]
        xhub = xall[i]
        rwall = router[i]
        rhub = rinner[i]

        # get the area of the annulus at the current x station
        xarea = pi * (rwall^2 - rhub^2)

        for j in 2:(nr - 1) #for each radial station between the hub and wall

            # get the annulus area at the current x station between the current radial stations
            darea = pi * (r_grid_points[1, j]^2 - r_grid_points[1, j - 1]^2)

            #calculate ratio of current annulus and rotor annulus areas
            arearatio = darea / rotorarea

            #define radial grid point based on conservation of mass
            r_grid_points[i, j] = sqrt(arearatio * xarea / pi + r_grid_points[i, j - 1]^2)

            #NOTE that the rest of the loop only does something if xs aren't the same like they are now.
            #calculate ratio for weighted average
            frac = (r_grid_points[i, j] - rhub) / (rwall - rhub)

            # place x grid point based on weighted ratio
            x_grid_points[i, j] = xhub + frac * (xwall - xhub)
        end # for radial stations
    end # for x stations

    return x_grid_points, r_grid_points
end

"""
    initialize_grid()

Initialize grid via zero-thrust, unit freestream solution.

**Arguments:**

**Returns:**
 - `grid::DuctTAPE.Grid` : Grid Object
"""
function initialize_grid()

    # get initial grid points
    #
    # relax grid
    #
    return grid
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
function relax_grid(xg, rg, nx, nr; max_iterations=-1)

    # initialize output and computation arrays
    xr = [xg[i, j] for i in 1:nx, j in 1:nr]
    rr = [rg[i, j] for i in 1:nx, j in 1:nr]
    D = [0.0 for i in 1:2, j in 1:nx]

    #set up relaxation factors
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

    #dfdc code does this, not sure how it's used yet:
    rpos = [0.0 for i in 1:nr]
    for j in 2:nr
        rpos[j] = rg[1, j]^2 - rg[1, j - 1]^2 + rpos[j - 1]
    end

    # Get the x positions along the wall (used later for various ratios compared to internal points as they move)
    xall = xg[:, 1]

    #next dfdc goes to AXELL function
    #skip over most of the stuff since the intlet and walls don't get relaxed

    for iterate in 1:max_iterations

        #loop over interior streamlines:
        for j in 2:(nr - 1)
            #rename indices for convenience
            jm1 = j - 1
            jp1 = j + 1

            ## -- Relax Outlet -- ##
            #outlet may or may not be a vertical line after relaxation...

            xej = xg[end, j]
            xem1j = xg[end - 1, j]
            xem2j = xg[end - 2, j]

            rej = rg[end, j]
            rem1j = rg[end - 1, j]
            rem2j = rg[ende - 2, j]

            xarc = xg[end, jp1] - xg[end, jm1] #x arc distance between j+1 and j-1 radial positions
            rarc = rg[end, jp1] - rg[end, jm1] #r arc distances between j+1 and j-1 radial positions

            xhat = vectx / sqrt(vectx^2 + vectr^2) # normalized x component of arc length (x direction)
            rhat = dy / sqrt(vectx^2 + vectr^2) # normalized r component of arc length (r direction)

            dxem1wall = xall[end] - xall[end - 1]
            dxem2wall = xall[end - 1] - xall[end - 2]

            #ratios of x and r spacing comparing current interior r position vs wall x position
            dxem1dxwall = (xej - xem1j) / dxem1wall
            dxem2dxwall = (xem1j - xem2j) / dxem2wall
            drem1drwall = (rej - rem1j) / dxem1wall
            drem2drwall = (rem1j - rem2j) / dxem2wall

            #2nd-order 3-point difference to get tangential velocity
            rez =
                xhat * (
                    dxem1dxwall +
                    dxem1wall * (dxem1dxwall - dxem2dxwall) / (dxem1wall + dxem2wall)
                ) +
                rhat * (
                    drem1dxwall +
                    dxem1wall * (drem1dxwall - drem2dxwall) / (dxem1wall + dxem2wall)
                )

            xi_xij = xhat * (1.0 / dxem1wall + 1.0 / (dxem1wall + dxem2wall))
            xi_rij = rhat * (1.0 / dxem1wall + 1.0 / (dxem1wall + dxem2wall))
            xi_arc = xi_xij * xhat + xi_rij * rhat

            arcrelax = 1.0 * relaxfactor
            darc = -arcrelax * rez / xi_arc

            ## -- Relax points on the jth streamline -- ##
            #TODO:NEED TO RECONCILE WRITTEN THEORY WITH CODE...

            for i in 2:(nx - 1)
                #rename indices for convenience
                im1 = i - 1
                ip1 = i + 1

                ravgplus = 0.5 * (r[i, j] + r[i, r + 1])
                ravgminus = 0.5 * (r[i, j] + r[i, j - 1])
                ravg = 0.5 * (ravgplus + ravgminus)

                dximinus = xall[i] - xall[i - 1]
                dxiplus = xall[i + 1] - xall[i]
                dxiavg = 0.5 * (dxijm1 + dxijp1)

                detaminus = rpos[j] - rpos[j - 1]
                detaplus = rpos[j + 1] - rpos[j]
                detaavg = 0.5(detaminus + detaplus)

                x_eta = 0.5 * (xg[i, j + 1] - xg[i, j - 1]) / detaavg
                r_eta = 0.5 * (rg[i, j + 1] - rg[i, j - 1]) / detaavg
                x_xi = 0.5 * (xg[i + 1, j] - xg[i - 1, j]) / dxiavg
                r_xi = 0.5 * (rg[i + 1, j] - rg[i - 1, j]) / dxiavg

                alpha = x_eta^2 + r_eta^2
                beta = x_eta * x_xi + r_eta * r_xi
                gamma = x_xi^2 + r_xi^2
                J = r_eta * x_xi - x_eta * r_xi

                ximinuscoeff = detaminus * detaplus / (dximinus * dxiavg)
                xipluscoeff = detaminus * detaplus / (dxiplus * dxiavg)
                etaminuscoeff = detaplus / detavg * ravgminus / ravg
                etapluscoeff = detaminus / detavg * ravgplus / ravg

                A =
                    alpha * (ximinuscoeff + xipluscoeff) +
                    gamma * (etaminuscoeff + etapluscoeff)

                if i == 2
                    B = 0
                else
                    B = -alpha * ximinuscoeff
                end

                C[i] = -alpha * xipluscoeff

                dxdxieta =
                    detaminus *
                    detaplus *
                    (
                        xg[i + 1, j + 1] - xg[i - 1, j + 1] - xg[i + 1, j - 1] +
                        xg[i - 1, j - 1]
                    ) / (4.0 * dxiavg * detaavg)

                drdxieta =
                    detaminus *
                    detaplus *
                    (
                        rg[i + 1, j + 1] - rg[i - 1, j + 1] - rg[i + 1, j - 1] +
                        rg[i - 1, j - 1]
                    ) / (4.0 * dxiavg * detaavg)

                x_xixi =
                    (xg[i + 1, j] - xg[i, j]) * xipluscoeff -
                    (xg[i, j] - xg[i - 1, j]) * ximinuscoeff

                x_etaeta =
                    (xg[i, j + 1] - xg[i, j]) * etapluscoeff -
                    (xg[i, j] - xg[i, j - 1]) * etaminuscoeff

                D[1, i] =
                    alpha * x_xixi - 2.0 * beta * dxdxieta + gamma * x_etaeta -
                    beta * r_xi * x_eta * detaminus * detaplus / ravg

                r_xixi =
                    (rg[i + 1, j] - rg[i, j]) * xipluscoeff -
                    (rg[i, j] - rg[i - 1, j]) * ximinuscoeff

                r_etaeta =
                    (rg[i, j + 1] - rg[i, j]) * etapluscoeff -
                    (rg[i, j] - rg[i, j - 1]) * etaminuscoeff

                D[2, i] =
                    alpha * r_xixi - 2.0 * beta * drdxieta + gamma * r_etaeta -
                    beta * r_xi * r_eta * detaminus * detaplus / ravg

                ainv = 1.0 / (A - B * C[i - 1])
                C[i] += ainv
                D[1, i] = (D[1, i] - B * D[1, i - 1]) * ainv
                D[2, i] = (D[2, i] - B * D[2, i - 1]) * ainv
            end #for i

            D[1, end] = 0.0
            D[2, end] = 0.0

            for ib in 2:(xn - 1)
                i = xn - ib + 1
                D[1, i] = D[1, i] - C[i] * D[1, i + p]
                D[2, i] = D[2, i] - C[i] * D[2, i + p]

                xr[i, j] = xg[i, j] + relaxfactor * D[1, i]
                rr[i, j] = rg[i, j] + relaxfactor * D[2, i]

                ad1 = abs(D[1, i])
                ad2 = abs(D[2, i])
                dmax = maximum(dmax, ad1, ad2)
            end #for i
        end #for j

        if dmax < tol * dxy
            return xr, rr
        else
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
    end

    return xr, rr
end
