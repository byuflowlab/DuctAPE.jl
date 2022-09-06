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
- Set Up
    - initialize velocities (VMAVGINIT)
    - generate gamma_theta solution (GTHCALC)
    - update wake gamma values
    - Solve system for initialized right hand side (GAMSOLV)
- for number of iterations
    - start at line 631 in rotoper.f

VMAVGINIT:
- very simple, just put guesses for vm_average on grid (dfdc doesn't do this very well, see todos in your implementation.

VMAVGCALC:

GTHCALC:
- calculate gamma_i's on vortex sheet wakes using eqn 45 from dfdc docs

GAMSOLV:
- CVPGEN (set up book keeping?)
- QAIC (this should be where all the coefficients are generated)
- SYSP (more book keeping?)
- GSYS (setup and factor system)
- SETROTORSRC (set rotor drag source strengths)
- SETDRGOBJSRC (same for drag objects, ignore for now)
- VMAVGINIT
- GSOLVE (solve based on Qinf, sigmas, and gamma_theta)
- QCSUM (get velocities at control points)
- QCPFOR (calculate cp's and forces on surfaces)
- SETGRDFLW (put data onto grid
- STGFIND (find stagnation points)

UPDROTVEL:

TQCALC:
=#

###############################
##### ----- EXPORTS ----- #####
###############################

## -- TYPES

export Outputs

## -- FUNCTIONS

#######################################
##### ----- COMPOSITE TYPES ----- #####
#######################################

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
need to get all the a's and b's (coefficient matrices) for the walls and rotors.
Not quite sure where this is done in dfdc, but the functions talking about pointers might be a good place to start (e.g. dfdcsubs.f line 1770)
"""
function initialize_system()

    # Initialize Rotor and Wake Aerodynamics

    system_aero, rotor_velocities, initial_vaxial_average = initialize_system_aerodynamics(
        rotors, blades, wakegrid, freestream; niter=10, rlx=0.5
    )

    # Initialize average V_m
    vm_average = initialize_vm_average(initial_vaxial_average)

    #see calculate gamma_theta_i values on wake vortex sheet panesl
    gamma_theta_wakes = calculate_gamma_theta(system_aero, vm_average)

    #NOTE: gamth = gth after this point in dfdc
    #
    #see GAMSOLV
    #
    # Create some sort of system object, make sure convergence flag is defaulted to false.
    #TODO: figure out what to return
    return system
end

"""
this is gamsolv at line 32 in solve.f in dfdc
"""
function solve_inviscid_panel()

    #CVPGEN (unneeded?)

    #QAIC (AIC matrix for velocities at control points)

    #SYSP (unneeded?)

    #GSYS (set up and factor system

    #SETROTORSRC (set rotor drag source strengths from velocities)

    #SETDRGOBJSRC (ignore for now)

    #GSOLVE (solve system)

    #QCSUM (get velocities at control points)

    #SCPFOR (get cps on surfaces and forces

    #SETGRDFLW (place flow information on wake grid)

    #STGFIND (find stagnation points

    return nothing
end

####################################
##### ----- NEWTON SOLVE ----- #####
####################################

"""
probably don't need to do things the way they are done in dfdc. There are better ways in julia.  Probably use the LinearSolve.jl package and whatever other packages convenient to get the derivatives as needed.
"""
function solve_system(; niter=100)

    #initialize system
    system = initialize_system()

    #iterate
    for iter in 1:niter

        #return with true convergence flag if converged
        if converged
            #update rotor velocities, see UPDROTVEL
            system.converged = true
            return system
        end
    end

    #return current system including a false convergence flag.
    return system
end
        rotor_panel_center = [(0.0, 0.0) for j in 1:numrotorsourcepan]

        #create a blade object for each rotor (to get dimensional radial data)
        blade = initialize_blade_dimensions(ductgeometry, ductsplines, rotors[i])

        #loop through each of the radial stations
        for j in 1:numrotorsourcepan

            # for now, the xlocation is the x coordinate of the panels
            #TODO: this will change when rake is added.
            rotor_panel_edge_x[j] = (gridxs[rotoridxs[i], j], gridxs[rotoridxs[i], j + 1])

            #use radial stations to define panel edges.
            rotor_panel_edge_r[j] = (gridrs[rotoridxs[i], j], gridrs[rotoridxs[i], j + 1])

            #as before, average the edges to get the centers
            rotor_panel_center[j] = (
                sum(rotor_panel_edge_x[j]) / 2.0, sum(rotor_panel_edge_r[j]) / 2.0
            )
        end

        if i > 1
            #TODO: Decide if it's better to do this here, or elsewhere.  Needs to be after any grid relaxation happens.
            # if more than one rotor, the rotor radial stations have more than likely changed, for aft rotors, and rotor information needs to be reinterpolated accordingly
            reinterpolate_rotor!(wakegrid, rotors[i], rotoridxs[i])
        end

        # create the rotor source panels object
        rotor_source_panels[i] = Panels(
            rotor_panel_edge_x, rotor_panel_edge_r, rotor_panel_center, "s"
        )
    end

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
    generate_panel_system(ductgeometry, ductsplines, rotors, wakegrid)

Put all the various panel objects together for convenience.

**Arguments:**
 - `ductgeometry::DuctTAPE.DuctGeometry` : Duct Geometry object
 - `ductsplines::DuctTAPE.DuctSplines` : Duct Splines object
 - `rotors::Array{DuctTAPE.Rotor}` : Array of rotor objects
 - `wakegrid::DuctTAPE.WakeGridGeometry` : Wake Grid object

**Returns:**
 - `panelsystem::DuctTAPE.PanelSystem` : All System Panels
"""
function generate_panel_system(ductgeometry, ductsplines, rotors, wakegrid)

    #generate individual panel objects
    wall_panels, hub_panels, wake_panels, rotor_source_panels = generate_paneling(
        ductgeometry, ductsplines, rotors, wakegrid
    )

    #return panel system containing all panel objects
    return PanelSystem(wall_panels, hub_panels, wake_panels, rotor_source_panels)
end

