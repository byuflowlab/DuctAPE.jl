# Airfoil Polar Corrections

In some cases various airfoil polar corrections may be required.
Of specific note are modifications to airfoil polars for post-stall behavior.
Thus far, DuctAPE is much more robust if the post-stall behavior in the lift polars does not exhibit a decrease in lift at angles of attack beyond that of the maximum lift coefficient.
Therefore a function is provided to help modify polars as needed:

```@docs
DuctAPE.C4Blade.stall_limiters
```

Various other correction methods are available, including the cascade corrections inherent in the [DuctAPE.C4Blade.DFDCairfoil](@ref) type.
The following methods are in addition to the various corrections available alongside the [CCBlade Airfoil Types](@ref).


```@autodocs
Modules = [DuctAPE.C4Blade]
Pages = ["C4Blade/airfoil_corrections.jl"]
```
