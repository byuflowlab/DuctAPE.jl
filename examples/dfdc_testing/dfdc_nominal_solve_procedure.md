# DFDC NOMINAL SOLVE PROCEDURE:

## at oper command:
1. call DFOPPER see dfdc.f
2. generate geometry
    - repanels bodies
    - gets control points and normals for body panels
    - generates source panels
    - generates wake panels
    - puts together all the ridiculous pointer arrays for bookkeeping.
3. wait for exec command.

## upon exec command:
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

## CONVGTHBG DETAILS:
1. 