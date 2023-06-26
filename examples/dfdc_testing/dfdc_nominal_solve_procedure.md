# DFDC NOMINAL SOLVE PROCEDURE:

## At 'oper' command:
1. call DFOPPER see dfdc.f
2. generate geometry (GENGEOM)
    - repanels bodies
    - gets control points and normals for body panels
    - generates rotor source panels
    - generates wake panels
        - NOTE: NEED TO LOOK INTO WAKE TE PANELS, what are these? what do they do?
    - puts together all the ridiculous pointer arrays for bookkeeping.
3. wait for exec command.

## Upon 'exec' command:
1. ROTINITBLD
    - "Sets reasonable initial circulation using current rotor blade geometry (chord, beta)."
    - in other words, initializes rotor.
2. SETGRDFLW
    - "Sets grid flow data from circulation and entropy on rotor lines"
    - in other words, initializes the wake.
3. CONVGTHBG
    - "Basic solver for GTH for a defined blade geometry. Uses underrelaxed iteration to converge GTH, BGAM"
    - in other words, this is where the magic happens (see below)
4. TQCALC
    - "Sets Thrust, Torque and their sensitivities wrt  QINF, OMEGA, BETA, chord(i), VA,VT, BGAM"
    - in other words, post processing

---

## CONVGTHBG DETAILS:
### PRE-ITERATION
1. VMAVGINIT
    - "Initializes VMAVG using axial velocity estimate"
2. GTHCALC
    - "Calculates GTH at rotor wake points using zero pressure jump relation. Note: this formulation sets GTH on endpoints using center velocity (not averaged VM at endpoints)"
3. GAMSOLV
    - "Generates inviscid panel solution with current flow condition"
    - in other words, this is the linear solve
        - note that they actually set up the system LHS matrix once when GAMSOLV is first called, then LU decompose that and solve the system with updated RHS's using backsubstitution.
    - if not called already calls
        - CVPGEN: set control point and put together associated "pointer" array
        - QAIC: sets induced velocities at control points
        - SYSP: sets up linear system
        - GSYS: defines linear system LHS
        - LUDCMP: LU decomposition of system LHS
    - after above, or after first call, then first calls SETROTORSRC
        - "Sets source strength for rotor profile drag"
    - then calls GSOLVE
        - "Does a direct solve from current RHS (knowns). Generates GAM(.),SIG(.),QNDOF(.)"
        - i.e. this is the actual linear solve
    - then calls QCSUM
        - "Computes velocities at control points for current panel and point-singularity strengths"
    - then calls QCPFOR
        - "Computes velocities, Cp, forces, etc. for current singularity distribitions"
    - then calls SETGRDFLW (see above) and STGFIND (finds stagnation point)

### INSIDE ITERATION
1. GAMSOLV
2. UPDROTVEL
    - "Update blade or disk velocities based on current solution.  Velocities updated include:
        - induced        velocities  VIND
        - absolute frame velocities  VABS
        - relative frame velocities  VREL"
3. for each of the rotors:
    1. do all the blade element and rotor circulation calculations
    2. update blade circulation using over relaxation ("CSOR", whatever that means.)
    2. SETGRDFLW (see above)
4. VMAVGCALC
    - Calculates VMAVG at rotor wake points using current center point velocities QC(1,.), QC(2,.)
5. GTHCALC (see above)

### AFTER ITERATION
1. UPDROTVEL

---

## GSOLVE Details:
- Again, the LHS is set up on the first call, LU decomposed, and used for all the linear solves (using back substitution), thereafter
- The RHS vector starts out as the solution of the previous solve (though it is originally initailzed with the other intial values).
    - the values are overwritten each iteration, starting over with the freestream normal velocity on each panel control point, then vortex and source influences are added.
- The system is N+1xN+1 where N is the sum of the number of hub and duct nodes, plus an additional 2 internal "panels" used for a zero internal tangnet-velocity constraint.
    - this gives N columns, but there are N-3 body panels associated with the rows.
    - there is a kutta condition applied to the duct
    - there is an internal-tangential velocity constraint applied to both the hub and the duct.
        - I am unclear on why this is required, or what it does.  The comments also talk about this being a QNDOF (normal velocity degree of freedom?)

---

## Notes:
- DFDC only converges based on the rotor(s) circulation, everything else is a function of those.
- DFDC doesn't use a Newton solve, as indicated in their documentation, it takes advantage of being able to use the previous iteration's values for the various singularity strengths, velocities, linear system RHS, etc., which greatly reduces the non-linear system size.
