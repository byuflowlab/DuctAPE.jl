#=
Functions for initialization of the ducted rotor system before the newton solve

Authors: Judd Mehr,

Procedure:
Iteration setup:
gengeom does the paneling definitions
1. Define paneling on grid streamsurfaces, define source drag panels

This could be the gsys function in solve.f (line 330ish)
Probably QAIC though
2. Evaluate at body panels, rotor blade stations, source drag panels: axij , arij , aij, bxij , brij ,bij

convgthbg appears to cover at least the rest of this including the newton iteration:
3. Set initial guess for Γk
4. Set corresponding Γ ̃ and H ̃ fields
5. Set initial guess for γi using (41) and (42)
6. Set initial σi = 0

One Newton iteration: DFDC doesn't actually do a newton solve, but rather a successive over relaxation
1. Using current Γk, σi, evaluate γi, vxi , vri , vθi , V⃗i and derivatives w.r.t. Γk, σi 2. Evaluate residuals of equations (75), (73), (74), and derivatives w.r.t. Γk, σi
3. Solve Newton system for δΓk, δσi
4. Update Γk, σi

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

TQCALC
=#

"""

Generates Wake Grid and System Paneling
"""
function initialize_geometry(ductgeometry, rotors, gridoptions)

    #TODO: update ductgeometry to include splines throughout.

    # Create Wake Grid
    wakegrid, updatedrotors = generate_wake_grid(ductgeometry, rotors, gridoptions)

    # Create Panels
    systempanels = generate_panel_system(ductgeometry, updatedrotors, wakegrid)

    return wakegrid, systempanels, updatedrotors
end

"""

Sets initial flow conditions for the system.

Runs a panel method on mirrored duct geometry to get initial guess for flow inside duct.
Runs CCBlade using panel method inflow to get first guess for rotor performance.
Then initializes guesses for the blade circulation, wake vorticity, and rotor source strengths.
"""
function initialize_flow_data(ductgeometry, systempanels, rotors, freestream)

    # Run Panel Method to get first guess on flow field.
    initial_rotor_velocities, Avinf, bvinf = first_guess_flowfield(systempanels, freestream)

    # Run CCBlade to get W and cl values
    W, cl, blades = first_guess_rotor_data(rotors, wall_panel_strengths, freestream)

    # Set Blade Circulation
    Gamma = set_blade_circulation(W, cl, blades)

    # Set Wake Enthalpy
    Delta_H = set_wake_enthalpy(Gamma, systempanels)

    # Set Wake Vorticity
    wake_vortex_strengths = set_wake_vorticity(Gamma, Delta_H)

    # Set Rotor Source Strengths
    rotor_source_strengths = set_rotor_panel_strengths()

    # Generate Grid Flow Data Object
    gridflowdata = GridFlowData(Gamma, Delta_H)

    # Generate Panel Strengths Object
    panelstrengths = PanelStrenghts(wall_panel_strengths, wake_vortex_strengths, rotor_source_strengths)

    return gridflowdata, panelstrengths
end

"""
    ### --- DFDC Process
    #CVPGEN (unneeded?)

    #QAIC (AIC matrix for velocities at control points)

    #SYSP (unneeded?)

    #GSYS (set up system)

    #LUDCMP (factor system to LU format)

    #SETROTORSRC (set rotor drag source strengths from velocities)

    #SETDRGOBJSRC (ignore for now)

    #GSOLVE (solve system)

    #QCSUM (get velocities at control points)

    #SCPFOR (get cps on surfaces and forces

    #SETGRDFLW (place flow information on wake grid)

    #STGFIND (find stagnation points
"""
function initialize_problem(ductgeometry, rotors, gridflowdata, freestream)

    wakegrid, systempanels, updatedrotors =  initialize_geometry(ductgeometry, rotors, gridoptions)

    gridflowdata, panelstrengths = initialize_flow_data(ductgeometry, systempanels, rotors, freestream)

    return problem
end
