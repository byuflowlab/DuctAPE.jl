# CCBlade Airfoil Types

DuctAPE includes stripped-down versions of most of the airfoil types as well as most of the polar correction methods available in [CCBlade](https://flow.byu.edu/CCBlade.jl/stable/reference/#Airfoil-Evaluation).
We have modified several of the CCBlade airfoil types to be compatible with the preallocated caching methods used in DuctAPE.
Note that any of these methods can be constructed directly or using file inputs, and they will still be compatible with the preallocation and internal reconstruction methods.

The first is an airfoil based only on angle of attack.  Its associated evaluation function gives the lift and drag as a function of angle of attack only. Thus the lift and drag fields in this type are vectors.

```@docs
DuctAPE.C4Blade.AlphaAF
```

The next two are provide lift and drag not only as a function of angle of attack but also Reynolds number or Mach number (note that these methods are not quite identical, and therefore they cannot be used interchangeably). The lift and drag fiedls in these types are matrices.

```@docs
DuctAPE.C4Blade.AlphaReAF
```

```@docs
DuctAPE.C4Blade.AlphaMachAF
```

The last CCBlade-like airfoil type combines the previous ones together. The evaluation function for this type give the lift and drag relative to the angle of attack, Reynolds number, and Mach number.  The lift and drag fields for this type are 3D arrays.

```@docs
DuctAPE.C4Blade.AlphaReMachAF
```

As mentioned we also have kept most of the airfoil correction methods and airfoil writing methods from CCBlade.  The [CCBlade documentation](https://flow.byu.edu/CCBlade.jl/stable/reference) is a good source for how to use them if desired.
