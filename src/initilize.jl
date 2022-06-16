#=
Functions for Grid Initialization

Author: Judd Mehr,

Process Notes:
1. Define paneling on duct and hub (FLOWFoil)
2. Evaluate for panel control points (FLOWFoil)
3. Assemble abarij, augmented with Kutta condition (FLOWFoil)
4. Solve for gammabar0i (FLOWFoil)
5. Set initial streamfunction grid (FLOWFoil)
6. Evaluate vx0i, vr0i at grid boundaries (FLOWFoil)
7. use vx0i, vr0i to get equipotential lines
8. Set Q1 = 0
9. Relax grid with SLOR (InterativeSolvers) using vx0i, vr0i for Neumann BCs

=#


