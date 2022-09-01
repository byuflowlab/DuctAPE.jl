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
# using IterativeSolvers
# is = IterativeSolvers #SOR solver

# - INCLUDED FILES
include("types.jl")
include("utils.jl")
include("walls.jl")
include("airfoils.jl")
include("rotors.jl")
include("wakegrid.jl")
include("panels.jl")
include("solve.jl")
end
