#=

Solver Functions

Author: Judd Mehr,

Procedure:
Iteration setup:
1. Define paneling on grid streamsurfaces, define source drag panels
2. Evaluate at body panels, rotor blade stations, source drag panels:
bij
3. Set initial guess for Γk
4. Set corresponding Γ ̃ and H ̃ fields
5. Set initial guess for γi using (41) and (42)
6. Set initial σi = 0
One Newton iteration:
axij , arij , aij, bxij , brij ,
1. Using current Γk, σi, evaluate γi, vxi , vri , vθi , V⃗i and derivatives w.r.t. Γk, σi 2. Evaluate residuals of equations (75), (73), (74), and derivatives w.r.t. Γk, σi
3. Solve Newton system for δΓk, δσi
4. Update Γk, σi
Streamline grid update (after Newton convergence)
1. Evaluate V⃗ at grid cells, using grid derivatives
2. Evaluate Γ ̃, H ̃, S ̃ for all grid cells downstream of blade stations and drag-producing objects (source panels)
3. Evaluate Q1 for all grid cells
4. Relax grid via SLOR
5. Begin again at “Iteration setup”

=#


#######################################
##### ----- ITERATION SETUP ----- #####
#######################################


