#=

Solver Functions

Authors: Judd Mehr,

Procedure:
Iteration setup:
gengeom does the paneling definitions
1. Define paneling on grid streamsurfaces, define source drag panels

This could be the gsys function in solve.f (line 330ish)
2. Evaluate at body panels, rotor blade stations, source drag panels: axij , arij , aij, bxij , brij ,bij

convgthbg appears to cover at least the rest of this including the newton iteration:
3. Set initial guess for Γk
4. Set corresponding Γ ̃ and H ̃ fields
5. Set initial guess for γi using (41) and (42)
6. Set initial σi = 0

One Newton iteration:
1. Using current Γk, σi, evaluate γi, vxi , vri , vθi , V⃗i and derivatives w.r.t. Γk, σi 2. Evaluate residuals of equations (75), (73), (74), and derivatives w.r.t. Γk, σi
3. Solve Newton system for δΓk, δσi
4. Update Γk, σi

---------------------------------------------------------
dfdc steps:
- Load file
- go to oper menu does the following:
    - preallocates arrays, save some convenience stuff (IGNORE THIS FOR NOW)
    - calls gengeom function (see below)
    - initializes plotting stuff (IGNORE)
    - waits for user input, if exec: (see line 961 in oper.f)
        - call rotinitbld
        - call setgrdflow
        - call convgthbg (solver for gamma_theta)
        - call tqcalc (solves thrust and torque and stuff)
        - call rotrprt (saves rotor state, probably IGNORE)

GENGEOM:
- adjusts walls paneling as needed (IGNORE)
- sets drag objects (IGNORE for now)
- initializes rotors [done]
- initializes grid [done]
- sets up wake elements on grid [done]
- sets up control points and "pointers" (this seems to simply be book keeping) TODO: need to figure out best way to do book keeping in julia since big global arrays/indices are not going to be good.

ROTORINTBLD:
- not too many pieces.  Just sets up reasonable initial blade circulations using momentum theory (does iterate to get self-induced stuff)

SETGRDFLOW:
- simply sets up grid flow data from circulation and entropy on rotors

CONVGTHBG:
- much more involved.  Lots of steps including linear system setup and solution.

TQCALC:
=#

#############################
##### ----- TYPES ----- #####
#############################

"""
    OperatingConditions{TVI,TVR,TF}

**Fields:**
 - `vinf::Array{Float}` : Array of freestream velocities
 - `vref::Array{Float}` : Array of reference velocities
 - `rho::Float` : air density value
 - `vso::Float` : speed of sound value
 - `mu::Float` : air viscosity value
"""
struct OperatingConditions{TVI,TVR,TF}
    vinf::TVI
    vref::TVR
    rho::TF
    vso::TF
    mu::TF
    # boundarylayer::TB
end

"""
"""
struct Panels{TPEx,TPEr,TPC,TPT}
    panel_edges_x::TPEx
    panel_edges_r::TPEr
    panel_centers::TPC
    panel_type::TPT
end

"""
"""
struct Outputs{TTt,TTd,TTr,TPt,TPr,TQt,TQr}
    totalthrust::TTt
    ductthrust::TTd
    rotorthrusts::TTr #array
    totalpower::TPt
    rotorpowers::TPr #array
    totaltorque::TQt
    rotortorques::TQr #array
end

#######################################
##### ----- ITERATION SETUP ----- #####
#######################################

"""
"""
function generate_paneling(ductgeometry, ductsplines, rotors, wakegrid)

    ## -- INITIALIZE ARRAYS -- ##
    # i.e. Count number of panels

    # - number of hub geometry panels
    #take number of xcoordinates for the hub and subtract the last one to get the number of panels
    numhubpan = length(ductgeometry.hubxcoordinates) - 1

    #if the hub has a blunt trailing edge. need to add another panel.
    if ductgeometry.hubbluntTE
        numhubpan += 1
    end

    # - number of wall geometry panels
    #take the number of x coordinates for the inner and outer wall geometry and subtract the last point from each to get the number of panels for the whole airfoil
    numwallpan =
        length(ductgeometry.wallinnerxcoordinates) +
        length(ductgeometry.wallouterxcoordinates) - 2

    #add another panel if the duct wall has a blunt trailing edge.
    #TODO: need to add bluntTE stuff to gridgeometry object creating function.
    if ductgeometry.wallbluntTE
        numwallpan += 1
    end

    # - number of rotor source panels
    #rename num_radial_stations for convenience
    numrs = length(rotors[1].radialstations)

    #number of rotor panels is going to be the number of rotors multiplied by one fewer than the number of radial stations
    numrotorsourcepan = length(rotors) * (numrs - 1)

    # - number of wake panels
    # TODO: need to make sure that panels don't double up for multiple rotors
    #each internal (not hub/tip) rotor radial station sheds a vortex sheet that extends to the end of the grid.
    numrotorwakepan = (numrs - 2) * (wakegrid.nx - 1)

    #the hub sheds a vortex sheet from its trailing edge to the end of the grid.
    numhubTEwakepan = wakegrid.nx - wakegrid.hubTEidx

    #the duct wall also sheds a vortex sheet from its trailing edge to the end of the grid.
    numwallTEwakepan = wakegrid.nx - wakegrid.wallTEidx

    # - number of drag object source panels (TODO LATER)

    ## -- GET PANEL EDGES, CENTERS AND TYPES -- ##

    # - wall panels
    #initialize from counts above.
    wall_panel_edge_x = [(0.0, 0.0) for i in 1:numwallpan]
    wall_panel_edge_r = [(0.0, 0.0) for i in 1:numwallpan]
    wall_panel_center = [(0.0, 0.0) for i in 1:numwallpan]

    # need to combine the inner and outer geometry coordinates
    wallxcoordinates = [
        ductgeometry.wallouterxcoordinates[1:(end - 1)]
        ductgeometry.wallinnerxcoordinates
    ] #TODO: need to make sure this is the correct direction. (what is the correct direction??)
    #TODO: also need to make sure that the inner and outer coordinates match up correctly (is LE and LE for both?) THIS IS A USER INPUT!!! PROBABLY NEED TO ADD CHECKS/CORRECTIONS FOR INCORRECT INPUT

    #similar for r coordinates
    wallrcoordinates = [
        ductgeometry.wallouterrcoordinates[1:(end - 1)]
        ductgeometry.wallinnerrcoordinates
    ]

    # loop through the total number of wall coordinates
    for i in 1:numwallpan
        if ductgeometry.wallbluntTE && i == numwallpan

            #create the trailing edge panel if needed
            wall_panel_edge_x[i] = (wallxcoordinates[end], wallxcoordinates[1])

            wall_panel_edge_r[i] = (wallrcoordinates[end], wallrcoordinates[1])

        else
            #get the x coordinates from the geometry and save as the panel edges.
            wall_panel_edge_x[i] = (wallxcoordinates[i], wallxcoordinates[i + 1])

            #similar for r coordinates
            wall_panel_edge_r[i] = (wallrcoordinates[i], wallrcoordinates[i + 1])
        end

        #no matter what case, get the average of the x and r coordinates to find the center point on the panel.
        wall_panel_center[i] = (
            sum(wall_panel_edge_x[i]) / 2.0, sum(wall_panel_edge_r[i]) / 2.0
        )
    end

    # create wall panels object
    wall_panels = Panels(wall_panel_edge_x, wall_panel_edge_r, wall_panel_center, "w")

    # - hub panels
    #initialize from counts above
    hub_panel_edge_x = [(0.0, 0.0) for i in 1:numhubpan]
    hub_panel_edge_r = [(0.0, 0.0) for i in 1:numhubpan]
    hub_panel_center = [(0.0, 0.0) for i in 1:numhubpan]

    #loop through hub panel count
    for i in 1:numhubpan
        if ductgeometry.hubbluntTE && i == numhubpan
            #define trailing edge panel if needed
            #
            hub_panel_edge_x[i] = (
                                   ductgeometry.hubxcoordinates[end], ductgeometry.hubxcoordinates[end]
            )

            hub_panel_edge_r[i] = (
                ductgeometry.hubrcoordinates[end], 0.0)

        else
            # get the x geometry coordinates and set as the panel edges
            hub_panel_edge_x[i] = (
                ductgeometry.hubxcoordinates[i], ductgeometry.hubxcoordinates[i + 1]
            )

            #similar for r coordinates
            hub_panel_edge_r[i] = (
                ductgeometry.hubrcoordinates[i], ductgeometry.hubrcoordinates[i + 1]
            )
        end

        #no mater the case, get the average of the x and r coordinates to set the center points of the panels
        hub_panel_center[i] = (
            sum(hub_panel_edge_x[i]) / 2.0, sum(hub_panel_edge_r[i]) / 2.0
        )
    end

    # create the hub panels object
    hub_panels = Panels(hub_panel_edge_x, hub_panel_edge_r, hub_panel_center, "w")

    # - rotor source panels
    #initialze based on the counts above
    rotor_panel_edge_x = [(0.0, 0.0) for i in 1:numrotorsourcepan]
    rotor_panel_edge_r = [(0.0, 0.0) for i in 1:numrotorsourcepan]
    rotor_panel_center = [(0.0, 0.0) for i in 1:numrotorsourcepan]

    #since there may be more than 1 rotor, set up a custom index
    rpanidx = 1

    #rename pertinent values for convenience
    rotoridxs = wakegrid.rotoridxs
    gridxs = wakegrid.x_grid_points
    gridrs = wakegrid.r_grid_points

    #loop through each rotor
    for i in 1:length(rotors)

        #create a blade object for each rotor (to get dimensional radial data)
        blade = initialize_blade_dimensions(ductgeometry, ductsplines, rotors[i])

        #loop through each of the radial stations
        for j in 1:(numrs - 1)

            # for now, the xlocation is the x coordinate of the panels
            #TODO: this will change when rake is added.
            rotor_panel_edge_x[rpanidx] = (
                gridxs[rotoridxs[i], j], gridxs[rotoridxs[i], j + 1]
            )

            #use radial stations to define panel edges.
            rotor_panel_edge_r[rpanidx] = (
                gridrs[rotoridxs[i], j], gridrs[rotoridxs[i], j + 1]
            )

            #as before, average the edges to get the centers
            rotor_panel_center[rpanidx] = (
                sum(rotor_panel_edge_x[rpanidx]) / 2.0,
                sum(rotor_panel_edge_r[rpanidx]) / 2.0,
            )

            #increment the custom index
            rpanidx += 1
        end

        if i > 1
            #TODO: Decide if it's better to do this here, or elsewhere.  Needs to be after any grid relaxation happens.
            # if more than one rotor, the rotor radial stations have more than likely changed, for aft rotors, and rotor information needs to be reinterpolated accordingly
            reinterpolate_rotor!(wakegrid, rotors[i], rotoridxs[i])
        end
    end

    # create the rotor source panels object
    rotor_source_panels = Panels(
        rotor_panel_edge_x, rotor_panel_edge_r, rotor_panel_center, "s"
    )

    # - wake panels

    #add together all the various wake panel counts from above
    numwakepan = numrotorwakepan + numhubTEwakepan + numwallTEwakepan

    #initialize from the total count
    wake_panel_edge_x = [(0.0, 0.0) for i in 1:numwakepan]
    wake_panel_edge_r = [(0.0, 0.0) for i in 1:numwakepan]
    wake_panel_center = [(0.0, 0.0) for i in 1:numwakepan]

    #loop through the hub trailing edge vortex sheet panels first.
    for i in 1:(numhubTEwakepan)

        #rename for convenience: position starts at the TE index
        hidx = wakegrid.hubTEidx + i - 1

        #use the grid coordinates from the hub TE to the end of the grid at the first radial position of the grid as panel edges.
        wake_panel_edge_x[i] = (
            wakegrid.x_grid_points[hidx, 1], wakegrid.x_grid_points[hidx + 1, 1]
        )

        #same for r positions
        wake_panel_edge_r[i] = (
            wakegrid.r_grid_points[hidx, 1], wakegrid.r_grid_points[hidx + 1, 1]
        )

        #again, average for centers.
        wake_panel_center[i] = (
            sum(wake_panel_edge_x[i]) / 2.0, sum(wake_panel_edge_r[i]) / 2.0
        )
    end

    #create a custom index for going forward since there's not an easier way to count.
    wpanidx = numhubTEwakepan + 1

    # loop through the grid radial stations
    for i in 2:(wakegrid.nr - 1)

        #loop through the grid x stations (from first rotor to end of wake)
        for j in 1:(wakegrid.nx - 1)

            # get the panel edges from the grid points
            wake_panel_edge_x[wpanidx] = (
                wakegrid.x_grid_points[j, i], wakegrid.x_grid_points[j + 1, i]
            ) #TODO: check that the indexing here is correct.  is it [r,x] or [x,r]?

            # same for r points
            wake_panel_edge_r[wpanidx] = (
                wakegrid.r_grid_points[j, i], wakegrid.r_grid_points[j + 1, i]
            ) #TODO: check that the indexing here is correct.  is it [r,x] or [x,r]?

            #again, average for centers
            wake_panel_center[wpanidx] = (
                sum(wake_panel_edge_x[wpanidx]) / 2.0, sum(wake_panel_edge_r[wpanidx]) / 2.0
            )

            #update custom index
            wpanidx += 1
        end
    end

    #loop through wall trailing edge vortex sheet panels
    for i in 1:(numwallTEwakepan)

        #rename index for convenience
        widx = wakegrid.wallTEidx + i - 1

        #get grid points at top radial position starting at duct wall trailing edge index.
        wake_panel_edge_x[wpanidx] = (
            wakegrid.x_grid_points[widx, end], wakegrid.x_grid_points[widx + 1, end]
        )

        #similar for r coordinates
        wake_panel_edge_r[wpanidx] = (
            wakegrid.r_grid_points[widx, end], wakegrid.r_grid_points[widx + 1, end]
        )

        #average for centers
        wake_panel_center[wpanidx] = (
            sum(wake_panel_edge_x[wpanidx]) / 2.0, sum(wake_panel_edge_r[wpanidx]) / 2.0
        )

        # increment custom index
        wpanidx += 1
    end

    # create wake (vortex sheet) panels object
    wake_panels = Panels(wake_panel_edge_x, wake_panel_edge_r, wake_panel_center, "v")

    #return the 3 types of panel objects.
    return wall_panels, hub_panels, wake_panels, rotor_source_panels
end

"""
need to get all the a's and b's (coefficient matrices) for the walls and rotors.
Not quite sure where this is done in dfdc, but the functions talking about pointers might be a good place to start (e.g. dfdcsubs.f line 1770)
"""
function initialize_linear_system() end

####################################
##### ----- NEWTON SOLVE ----- #####
####################################

"""
probably don't need to do things the way they are done in dfdc. There are better ways in julia.  Probably use the LinearSolve.jl package and whatever other packages convenient to get the derivatives as needed.
"""
function solve_linear_system() end
