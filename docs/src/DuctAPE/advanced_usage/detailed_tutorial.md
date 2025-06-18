
!!! note "Airfoils"
    Airfoil types for DuctAPE are currently contained in the C4Blade (Cascade Compatible [CCBlade](https://flow.byu.edu/CCBlade.jl/stable/)) sub-module of DuctAPE which is exported as `c4b` and also contains the various airfoil evaluation functions used for the blade element lookups.
    The available airfoil types include all the airfoil types from CCBlade, as well as `DFDCairfoil` which is an [XROTOR](https://web.mit.edu/drela/Public/web/xrotor/)-like parametric cascade polar used in DFDC.
    In addition there are untested cascade types with similar structure to CCBlades airfoil types called `DTCascade`.
    Furthermore, there is an experimental actuator disk model implemented via the `ADM` airfoil type in C4Blade.

