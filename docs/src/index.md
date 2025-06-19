# DuctAPE.jl [[Duct](#)ed [A](#)xisymmetric [P](#)ropulsor [E](#)valuation]

 [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://flow.byu.edu/DuctAPE.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/DuctAPE.jl/dev)
[![Build Status](https://github.com/byuflowlab/DuctAPE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/byuflowlab/DuctAPE.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

DuctAPE is a code for the aerodynamic evaluation of axisymmetric ducted rotors designed for incompressible (low mach) applications.
It is strongly influenced by the underlying [theory](https://web.mit.edu/drela/Public/web/dfdc/DFDCtheory12-31.pdf) of Ducted Fan Design Code [(DFDC)](https://web.mit.edu/drela/Public/web/dfdc/), utilizing a linear axisymmetric vortex panel method for duct and center body, blade element lifting line rotor representation, and psuedo wake-screw wake model axisymmetrically smeared onto an elliptic grid for efficient computation.

DuctAPE has been developed specifically for applications in gradient-based optimization settings.
The default solver methods have been chosen to balance code efficiency as well as robustness while simultaneously allowing for efficient automatic differentiation through DuctAPE employing [ImplicitAD.jl](https://flow.byu.edu/ImplicitAD.jl/dev/).


## Installation

```julia
pkg> add DuctAPE
```

## Publications

[Mehr, J., and Ning, A., “DuctAPE: A steady-state, axisymmetric ducted fan analysis code designed for gradient-based optimization.,” AIAA Aviation Forum, Las Vegas, Jul. 2024. doi:10.2514/6.2024-4297](https://scholarsarchive.byu.edu/facpub/7214/)
