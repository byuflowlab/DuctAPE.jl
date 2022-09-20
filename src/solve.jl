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

#######################################
##### ----- ITERATION SETUP ----- #####
#######################################

"""
this is gamsolv at line 32 in solve.f in dfdc
"""
function solve_inviscid_panel()

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

    #use nlsolve to solve newton system.
    nlsolve(f!, initial_x, autodiff = :forward)

    #return current system including a false convergence flag.
    return system
end
