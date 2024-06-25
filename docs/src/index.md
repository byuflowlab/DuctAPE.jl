# DuctAPE.jl [[Duct](#)ed [A](#)xisymmetric [P](#)ropulsor [E](#)valuation]

Authors: Judd Mehr,

Contributers: Taylor McDonnell,

DuctAPE is a code for the aerodynamic evaluation of axisymmetric ducted propulsors designed for incompressible (low mach) applications.
It is strongly influenced by the underlying [theory](https://web.mit.edu/drela/Public/web/dfdc/DFDCtheory12-31.pdf) of Ducted Fan Design Code [(DFDC)](https://web.mit.edu/drela/Public/web/dfdc/), utilizing a linear axisymmetric vortex panel method for duct and center body, blade element lifting line rotor representation, and psuedo wake-screw wake model axisymmetrically smeared onto an elliptic grid for efficient computation.

DuctAPE has been developed specifically for applications in gradient-based optimization settings.
The selected solver methods have been chosen to balance code efficiency as well as robustness while simultaneously allowing for efficient automatic differentiation through DuctAPE employing [ImplicitAD.jl](https://flow.byu.edu/ImplicitAD.jl/dev/).
At the same time, the basic functionality of a DFDC-like solve approach has been maintained for the interested user.

## Package Features

-

## Installation

As DuctAPE is not yet a registered package, if you have access to the development repository at this time you can add the package through Julia's package manager as:

```julia
pkg> add "https://github.com/byuflowlab/DuctAPE.jl"
```
