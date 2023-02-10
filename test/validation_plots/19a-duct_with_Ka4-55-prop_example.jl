#=

Validation using MARIN experimental data with a 19a Duct and a Ka 4-55 Propeller based on plots from Lewis and Ryan & Glover.

This script also serves as an example of how to run the DuctTAPE code for ducted propeller analysis.

=#

using DuctTAPE
const dt = DuctTAPE

include("../../plots_default.jl")

# - Duct Geometry - #
include("../data/marin_19a_duct.jl")

duct_coordinates = [full_x full_r]

# - Hub Geometry - #
# TODO: need to figure out hub geometry...

# - Rotor Geometry - #
# TODO: need to figure out rotor geometry...


