#=
in order to generate all the required coefficient matrices, we need to first make sure each geometric item has an associated panels object.
we then need to assemble the required meshes.
we can then put together the coefficient matrices.

needed matrices include:
- body to body (default in FLOWFoil), for system linear solve
- wake to body, for full system linear solve
- rotor to body, for full system linear solve
- body to rotor, for induced velocity calculation
- wake to rotor, for induced velocity calculation
- rotor to rotor, for induced velocity calculation

we only need flowfoil boundary conditions for the body to body case, so we don't want to use the system generation functions from flowfoil, but rather just the vortex and source coefficient functions.


=#


