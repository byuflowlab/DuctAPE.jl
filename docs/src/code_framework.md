# Solver Overview

## Setup

1. [Define Duct/Hub Geometry](#1-define-ducthub-geometry)
1. [Define Rotor non-dimensional Geometry](#2-define-rotor-non-dimensional-geometry)
1. [Define Wake Grid Geometry](#3-define-wake-grid-geometry)
1. [Define Paneling](#4-define-paneling)
1. [Set first guess for rotor flow conditions](#5-set-first-guess-for-rotor-flow-conditions)
1. [Get first guess rotor circulations](#6-get-first-guess-rotor-circulations)
1. [Get first guess grid flow data](#7-get-first-guess-grid-flow-data)
1. [Get first guess wake vortex sheet strengths](#8-get-first-guess-wake-vortex-sheet-strengths)
1. [Calculate velocities and coefficients](#9-calculate-velocities-and-coefficients)
1. [Calculate influence matrix (LHS of system)](#10-calculate-influence-matrix-lhs-of-system)
1. [Calculate blade section source strengths](#11-calculate-blade-section-source-strengths)
1. [Get initial solution for $\overline{\gamma}$](#12-get-initial-solution-for-overlinegamma)
1. [Update control point velocities](#13-update-control-point-velocities-sum-influence-of-wall-panels-wake-panels-and-freestream)
1. [Update grid flow data](#14-update-grid-flow-data)

Once we have an initial solution set, we're ready to start the interative solution process.

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

_DFDC: start in GENGEOM, (line 608 in dfdcsubs.f). This function is called right away after heading into the OPER menu. (the geometry is loaded from the case file first though, which is basically just reading in the data.)_

### **1. Define Duct/Hub Geometry**

DuctTAPE: Duct wall and hub (walls) are defined using `defineDuctGeometry` (in walls.jl)


### **2. Define Rotor non-dimensional Geometry**

DuctTAPE: Rotors are non-dimensionally defined using `initialize_rotor_geometry` (in rotors.jl)


### **3. Define Wake Grid Geometry**

_DFDC: see INIGRD in inigrd.f (line 32)_

DuctTAPE: The rotor wake grid is defined using `initialize_wake_grid` (in wakegrid.jl) using more or less the same methodology used in DFDC.

Notes:
 - Ideally this step should be done using an inviscid panel solution from the duct/hub geometry, which is what the DFDC theory document indicates, but instead a simple constant freestream and conservation of mass is used in practice.
 - If wanting to use an inviscid panel solution for this initialization step, need to figure out how to handle the open hub wall geometry. 


### **4. Define Paneling**

_DFDC: The rotor source panel and wake vortex panel portions of this is done in DFDC in SETROTWAK in dfdcsubs.f (line 1509) which calls ADDWAKE (line 1637 in the same file).
The wall panel portion of this is started in DFDC in ADJPANL in adjpanl.f (line 32) then finised in CVPGEN in geom.f (line 117)_

DuctTAPE: The wall and wake panels are defined in `generate_paneling` (in panels.jl)


### **5. Set first guess for rotor flow conditions**

_DFDC: see ROTINITBLD (rotor initialize blade) in rotoper.f (line 174)_

DuctTAPE: As in the DFDC code, the function that initializes the rotor flow conditions does a few other things, thus we called our version `initialize_system_aerodynamics` in system.jl

Notes:
- Could probably use CCBlade directly for this
- Currently the freestream velocity (constant across the blade) is used to get rotor operating point.
    - Would probably be a more accurate starting point to use an inviscid panel solution to get distributed, rather than average velocities across blade.
    - Also wouldn't be too much extra work, but still have the problem of not knowing how to get an invisid panel solution with the open hub geometry.


### **6. Get first guess rotor circulations**

_DFDC: see ROTINITBLD as well_

DuctTAPE: This is also just part of `initialize_system_aerodynamics`

Notes:
 - It would be easy to get these from CCBlade outputs (W and cl) if using CCBlade for the initialization process.


### **7. Get first guess grid flow data**

Specifically we want to get the circulations ($\widetilde{\Gamma}$), enthalpies ($\widetilde{H}$), and average $V_m$ at each wake grid point.

_DFDC: see ROTBG2GRD in inigrd.f (line 434)  (note this is very similar to SETGRDFLW starting in line 357 in the same file.)_

DuctTAPE: As just noted above, there are 2 nearly identical DFDC functions that sets the grid values from the rotor blade element values.  We have chosen to combine all these using multiple dispatch under `set_grid_aero!`, which depending on the number of inputs will select whether or not to update the entropy values (ROTBG2GRD does not, SETGRDFLW does) which requires more information than is available at this point of the initialization, so the entropy information is set up later.  

Notes:
- The $V_m$ calculations are just part of the `initialize_system_aerodynamics` function.  Only the grid circulation and enthalpy values are defined with `set_grid_aero`.


### **8. Get first guess wake vortex sheet strengths**

Specifically, from grid circulations ($\widetilde{\Gamma}$) and average $V_m$, get first guess wake vortex sheet strengths ($\gamma_\theta$) at each wake panel center.

_DFDC: see GTHCALC in rotoper.f (line 1878)_

DuctTAPE: `calculate_gamma_theta`

Notes:
- DFDC does something odd in defining ($\gamma_\theta$) on the walls. It basically does an interpolated taper from the ($\gamma_\theta$) at the trailing edge to zero at the first rotor station. 
I am not sure why it does that at this point, but it's probably to help some setup later.
It could simply be that filler values are required for the calculations, or perhaps it makes things easier since the grid technically lies on the walls since there are wakes eminating from the wall trailing edges.


### **9. Calculate velocities and coefficients**

_DFDC: see QAIC in qaic.f starting on line 31_

DuctTAPE: TODO

Notes:

### **10. Calculate influence matrix (LHS of system)**

_DFDC: see GSYS in system.f (line 328)_

DuctTAPE: TODO

Notes:


### **11. Calculate blade section source strengths**

_DFDC: see SETROTORSRC in rotoper.f (line 1615)_

DuctTAPE: TODO

Notes:
- Inputs: control point velocities and wake circulations

### **12. Get initial solution for $\overline{\gamma}$**

_DFDC: see GSOLVE in system.f (line 1504)_

DuctTAPE: TODO

Notes:


### **13. Update control point velocities (sum influence of wall panels, wake** panels, and freestream)

_DFDC: see QCSUM in system.f (line 1795)_

DuctTAPE: TODO

Notes:


### **14. Update grid flow data**

($\widetilde{H}$, $\widetilde{S}$ ,$\widetilde{\Gamma}$)

_DFDC: see SETGRDFLW in inigrd.f (line 357)  (see comments above about similarities to [ROTBG2GRD](#7-get-first-guess-grid-flow-data))_

DuctTAPE: TODO

Notes:
- Note that the S values aren't used until the post processing calculation of forces.


----------------------------------------------------------------------------------



## Solution
The following steps are all preformed in a loop, until converged or the iteration limit is met:

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


----------------------------------------------------------------------------------



## Post Process

1. Calculate surface pressure and forces

DFDC: see QCPFOR in vels.f (line 32)

Notes:
- this is actually done in every iteration in dfdc, but it doesn't look like it needs to be called every time, since it doesn't update anything, just post processes.