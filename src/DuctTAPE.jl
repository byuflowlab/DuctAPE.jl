module DuctTAPE

# - DEPENDENCIES
# using FLOWFoil
# ff = FLOWFoil #panels for walls/hub
using FLOWMath
fm = FLOWMath #akima splines
using Statistics
using CCBlade
ccb = CCBlade
# using ForwardDiff
# fd = ForwardDiff #newton method
# using NLsolve

# - INCLUDED FILES
include("types.jl")
include("walls.jl")
include("rotors.jl")
include("wakegrid.jl")
include("panels.jl")
include("system.jl")

include("utils.jl")
include("airfoils.jl")

# include("initailize.jl")
# include("solve.jl")

#######################
##### - EXPORTS - #####
#######################

##  TYPES

export Freestream
export DuctGeometry, DuctSplines
export GridOptions, WakeGridGeometry
export PanelSystem, Panels
export RotorGeometry, BladeDimensions
export Outputs

##  FUNCTIONS

export defineDuctGeometry, split_wall
export defineGridOptions, initialize_wakegrid
export generate_paneling, generate_panel_system
export initialize_system_aerodynamics
export initialize_rotor_geometry, initialize_blade_dimensions

end
