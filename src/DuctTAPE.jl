module DuctTAPE

## -- DEPENDENCIES

# using FLOWFoil
using FLOWMath
using CCBlade
ccb = CCBlade #rename ccblade for convenience

## -- INCLUDES

include("utils.jl")
include("types.jl")
include("couple.jl")
include("post_process.jl")

end
