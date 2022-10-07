# DuctTAPE.jl ([Duct]()ed [T]()wo-dimensional [A]()ero [P]()ropulsor [E]()valuation)

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://flow.byu.edu/DuctTAPE.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/DuctTAPE.jl/dev)
[![Build Status](https://github.com/byuflowlab/DuctTAPE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/byuflowlab/DuctTAPE.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

DuctTAPE is a code for the aerodynamic evaluation of 2D, axisymmetric, ducted propulsors design for incompressible aerodynamic applications (although hydrodynamic applications could also apply, but "Aero" made the acronym work).

Currently, this code is a simple coupler from FLOWFoil (an inviscid panel code for both 2D and axisymmetric applications) to [CCBlade](https://flow.byu.edu/CCBlade.jl/stable/).
Future plans include full rotor-duct coupling using methods similar to those used in [DFDC](http://web.mit.edu/drela/Public/web/dfdc/) with modifications to enable ease of use in gradient-based optimization applications.


## Warning: "Unregistered Dependencies"
Currently, FLOWFoil is an unregistered dependency of FLOWFoil.
For this reason, all automated github testing fails.

In addition, in order to allow the documentation to be automatically deployed, we commented out the FLOWFoil line in both `src/DuctTAPE.jl` and `Project.toml`.  In order for DuctTAPE to function, those line will need to be uncommented by the user.
