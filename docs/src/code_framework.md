# Full Setup and Solution Overview

## Setup

1. [Define Duct/Hub Geometry](#define-ducthub-geometry)
1. [Define Rotor non-dimensional Geometry](#define-rotor-non-dimensional-geometry)
1. [Define Wake Grid Geometry](#define-wake-grid-geometry)
1. [Define Paneling](#define-paneling)
1. Set first guess for rotor flow conditions using freestream
1. From flow conditions, get first guess rotor circulations
1. From rotor circulations, get first guess wake circulations ($\widetilde{\Gamma}$), enthalpies ($\widetilde{H}$), and average $V_m$ at each wake grid point.
1. from wake circulations and average $V_m$, get first guess wake vortex sheet strengths ($\gamma_\theta$) at each wake panel center.
1. From wall and wake geometry, calculate influence coefficients of wall and wake panels on wall panels. (qaic, qaic.f)
1. From $\gamma_\theta$ and influence coefficients (and freestream), calculate velocities at wake control points (also qaic, qaic.f)
1. From control point velocities and wake circulations, calculate blade source strengths. (setrotorsrc, rotoper.f line 1615)
1. get initial solution for $\overline{\gamma}$
1. update control point velocities (sum influence of wall panels, wake panels, and freestream)
1. update grid flow data (tilde values: H, S ,$\Gamma$)


## Solution

 - In loop, until converged or hit iteration limit:
    1. Solve for $\overline{\gamma}$
    1. update control point velocities
    1. update wake grid flow data (all the tilde things in dfdc theory (H, S, $\Gamma$)
    1. Update rotor velocities
    1. find new airfoil data based on updated velocities
    1. find changes in blade circulation
    1. get relaxation factors based on changes in blade circulation
    1. update blade circulation with relaxed values
    1. update $V_m$ values on wakes
    1. update $\gamma_\theta$ values with relaxed values

## Post Process

1. calculate surface pressure and forces (qcpfor in vels.f line 32)
    - Note that this is actually done in every iteration in dfdc, but it doesn't look like it needs to be called every time, since it doesn't update anything, just post processes.

------------------------------------------------------------



# Details

## Setup

DFDC: start in GENGEOM, (line 608 in dfdcsubs.f)

### 1. Define Duct/Hub Geometry

DuctTAPE: Duct wall and hub (walls) are defined using `defineDuctGeometry` (in walls.jl)


### 2. Define Rotor non-dimensional Geometry

DuctTAPE: Rotors are non-dimensionally defined using `initialize_rotor_geometry` (in rotors.jl)


### 3. Define Wake Grid Geometry

DFDC: see INIGRD in inigrd.f (line 32)

DuctTAPE: The rotor wake grid is defined using `initialize_wake_grid` (in wakegrid.jl) using more or less the same methodology used in DFDC.

Notes:
 - Ideally this step should be done using an inviscid panel solution from the duct/hub geometry, which is what the DFDC theory document indicates, but instead a simple constant freestream and conservation of mass is used in practice.
 - If wanting to use an inviscid panel solution for this initialization step, need to figure out how to handle the open hub wall geometry. 


### 4. Define Paneling

DFDC: The rotor source panel and wake vortex panel portions of this is done in DFDC in SETROTWAK in dfdcsubs.f (line 1509) which calls ADDWAKE (line 1637 in the same file).
The wall panel portion of this is started in DFDC in ADJPANL in adjpanl.f (line 32) then finised in CVPGEN in geom.f (line 117)

DuctTAPE: The wall and wake panels are defined in `generate_paneling` (in panels.jl)


### 5. Set first guess for rotor flow conditions using freestream

DFDC:

DuctTAPE:

Notes:
- Could probably use CCBlade directly for this
- would probably be more accurate to also use the inviscid panel solution to get distributed, rather than average velocities across blade.


### 6. From flow conditions, get first guess rotor circulations
DFDC:

DuctTAPE:

Notes:
 - Easy to get these from CCBlade outputs (W and cl)

### 7. From rotor circulations, get first guess wake circulations ($\widetilde{\Gamma}$), enthalpies ($\widetilde{H}$), and average $V_m$ at each wake grid point.
DFDC:

DuctTAPE:

Notes:

### 8. from wake circulations and average $V_m$, get first guess wake vortex sheet strengths ($\gamma_\theta$) at each wake panel center.
DFDC:

DuctTAPE:

Notes:

### 9. From wall and wake geometry, calculate influence coefficients of wall and wake panels on wall panels. (qaic, qaic.f)
DFDC:

DuctTAPE:

Notes:


### 10. From $\gamma_\theta$ and influence coefficients (and freestream), calculate velocities at wake control points (also qaic, qaic.f)
DFDC:

DuctTAPE:

Notes:

### 11. From control point velocities and wake circulations, calculate blade source strengths. (setrotorsrc, rotoper.f line 1615)
DFDC:

DuctTAPE:

Notes:

### 12. get initial solution for $\overline{\gamma}$
DFDC:

DuctTAPE:

Notes:

### 13. update control point velocities (sum influence of wall panels, wake panels, and freestream)
DFDC:

DuctTAPE:

Notes:


### 14. update grid flow data (tilde values: H, S ,$\Gamma$)
DFDC:

DuctTAPE:

Notes:
- Note that the S values aren't used until the post processing calculation of forces.





## Solution

 - In loop, until converged or hit iteration limit:
    1. Solve for $\overline{\gamma}$
    1. update control point velocities
    1. update wake grid flow data (all the tilde things in dfdc theory (H, S, $\Gamma$)
    1. Update rotor velocities
    1. find new airfoil data based on updated velocities
    1. find changes in blade circulation
    1. get relaxation factors based on changes in blade circulation
    1. update blade circulation with relaxed values
    1. update $V_m$ values on wakes
    1. update $\gamma_\theta$ values with relaxed values

## Post Process

1. calculate surface pressure and forces (qcpfor in vels.f line 32)
    - Note that this is actually done in every iteration in dfdc, but it doesn't look like it needs to be called every time, since it doesn't update anything, just post processes.