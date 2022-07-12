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
    #TODO add expansion factor for wake
    numxw = grid_options.num_xwake_stations
    xw = range(wake_range[1], wake_range[2]; length=numxw)

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
    relax_grid(x_grid_points, r_grid_points)

Relax grid using elliptic grid solver.

**Arguments:**
 - `x_grid_points::Matrix{Float64}` : Initial x grid points guess
 - `r_grid_points::Matrix{Float64}` : Initial r grid points guess

**Returns:**
 - `x_relax_points::Matrix{Float64}` : Relaxed x grid points guess
 - `r_relax_points::Matrix{Float64}` : Relaxed r grid points guess
"""
function relax_grid(x_grid_points, r_grid_points)

    # initialize arrays

    return nothing
end
