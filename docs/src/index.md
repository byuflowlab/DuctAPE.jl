# DuctAPE.jl [[Duct](#)ed [A](#)xisymmetric [P](#)ropulsor [E](#)valuation]

Authors: Judd Mehr,

Contributers: Taylor McDonnell,

DuctAPE is a code for the aerodynamic evaluation of axisymmetric ducted propulsors designed for incompressible (low mach) applications.
It is strongly influenced by the underlying [theory](https://web.mit.edu/drela/Public/web/dfdc/DFDCtheory12-31.pdf) of Ducted Fan Design Code [(DFDC)](https://web.mit.edu/drela/Public/web/dfdc/), utilizing a linear axisymmetric vortex panel method for duct and center body, blade element lifting line rotor representation, and wake model axisymmetrically smeared onto an elliptic grid for efficient computation.
DuctAPE has been developed specifically for applications in gradient-based optimization settings.


### Installation

```julia
pkg> add DuctAPE
```

### Documentation

- Start with [Getting Started](@ref) to get up and running.
- The Advanced Usage tab includes several pages of additional information for customizing your usage.
- The API tab contains public and private method descriptions.
- The Theory tab contain several pages on the underlying theory of DuctAPE.
- The C$^4$Blade tab contains documentation for the C$^4$Blade submodule used for airfoil/cascade management within DuctAPE as well as state initialization.

### Citing

Mehr, J. and Ning, A., "DuctAPE: A steady-state, axisymmetric ducted fan analysis code designed for gradient-based optimization.," AIAA Aviation Forum, July 2024.
