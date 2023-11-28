#=

Validation using MARIN experimental data with a 19a Duct and a Ka 4-55 Propeller based on plots from Lewis and Ryan & Glover.

This script also serves as an example of how to run the DuctAPE code for ducted propeller analysis.

=#

using DuctAPE
const dt = DuctAPE

include("../../plots_default.jl")

# - Duct Geometry - #
include("../data/marin_19a_duct.jl")

duct_coordinates = [full_x full_r]

# - Hub Geometry - #
# TODO: need to figure out hub geometry...

# - Rotor Geometry - #
# TODO: need to figure out rotor geometry...

# - Operating Conditions - #

# advance_ratio = 0.224 #advance ratio for duct surface pressure
advance_ratio = 0.36 #advance ratio for duct surface pressure
# advance_ratio = 0.551 #advance ratio for duct surface pressure
Js = range(0.2, 0.6; step=0.1) #vector of advance ratios for thrust

# Choose reasonable rotation rate

# use advance ratio and rotor radius to calculate correct inflow velocity

# put together freestream
rho = 1.225
mu = 1.81e-5
asound = 340.0
dt.Freestream(Vinf, rho, mu, asound) #note: this is marine prop, maybe use water values?

##### ----- Surface Pressure Validation ----- #####

## - Assemble Rotor Parameters - ##

# rotor_parameters = [(
#     rotor_x_position=rotor_x_position,
#     radial_positions=rotor_radial_positions,
#     chords=chords,
#     twists=twists,
#     airfoils=airfoils,
#     num_radial_stations=num_radial_stations,
#     num_blades=num_blades,
#     Omega=Omega,
# )]

## - Run Solver - ##

state_variables, converged, params = dt.analyze_propulsor(
    duct_coordinates, hub_coordinates, rotor_parameters
)

## - Post Process Solution - ##

##### ----- Thrust Validation ----- #####

# loop through advance ratios
