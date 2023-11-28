# DuctAPE.jl [[Duct](#)ed [T](#)wo-dimensional ([A](#)xisymmetric) [P](#)ropulsor [E](#)valuation]

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://flow.byu.edu/DuctAPE.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/DuctAPE.jl/dev)
[![Build Status](https://github.com/byuflowlab/DuctAPE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/byuflowlab/DuctAPE.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

DuctAPE is a code for the aerodynamic evaluation of axisymmetric ducted propulsors designed for incompressible *(low Mach) applications.


## Warning: "Unregistered Dependencies"
Currently, FLOWFoil is an unregistered dependency of DuctAPE.
For this reason, all automated github testing fails.

In addition, in order to allow the documentation to be automatically deployed, we commented out the FLOWFoil line in both `src/DuctAPE.jl` and `Project.toml`.  In order for DuctAPE to function, those line will need to be uncommented by the user.
