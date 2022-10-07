```@meta
CurrentModule = DuctTAPE
```

# DuctTAPE ([Duct](#)ed [T](#)wo-dimensional [A](#)ero[P](#)ropulsor [E](#)valuation)

DuctTAPE is a code for the aerodynamic evaluation of 2D, axisymmetric, ducted propulsors design for incompressible aerodynamic applications (although hydrodynamic applications could also apply, but "Aero" made the acronym work).

Currently, this code is a simple coupler from FLOWFoil (an inviscid panel code for both 2D and axisymmetric applications) to [CCBlade](https://flow.byu.edu/CCBlade.jl/stable/).
Future plans include full rotor-duct coupling using methods similar to those used in [DFDC](http://web.mit.edu/drela/Public/web/dfdc/) with modifications to enable ease of use in gradient-based optimization applications.


!!! note "In Order To Use"
    In order to make this documentation work, the non-registered FLOWFoil dependency had to be commented out of the `src/DuctTAPE.jl` and `Project.toml` files.

    In order to use this code without errors, those lines need to be uncommented.

    This is also why none of the automated testing currently passes, and would not anyway with an unregistered dependency.
