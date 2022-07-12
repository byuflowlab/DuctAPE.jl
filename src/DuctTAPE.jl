module DuctTAPE

# - DEPENDENCIES
using FLOWFoil
ff = FLOWFoil #panels for walls/hub
using FLOWMath
fm = FLOWMath #akima splines
# using ForwardDiff
# fd = ForwardDiff #newton method
# using IterativeSolvers
# is = IterativeSolvers #SOR solver

# - EXPORTS

#TYPES

#FUNCTIONS

# - INCLUDED FILES
include("types.jl")
include("initilize.jl")
include("duct.jl")

end
