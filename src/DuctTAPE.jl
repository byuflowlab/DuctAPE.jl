module DuctTAPE

## -- DEPENDENCIES

# using FLOWFoil #UNCOMMENT THIS LINE TO ENABLE FLOWFOIL
using FLOWMath
using CCBlade
ccb = CCBlade #rename ccblade for convenience

## -- INCLUDES

include("utils.jl")
include("types.jl")
include("couple.jl")

end
