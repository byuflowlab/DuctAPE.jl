#=
Duct Setup Functions
=#

"""
    split_wall(x,r)

Splits full airfoil coordinates into upper and lower halves.

Only works based geometry.  Splits at lowest x-value, does not split based on stagnation point.

**Arguments:**
 - `x::Array{Float}` : Array of x-coordinates, assumed to start at the bottom trailing edge and proceed clockwise
 - `r::Array{Float}` : Array of r-coordinates, assumed to start at the bottom trailing edge and proceed clockwise

**Returns:**
 - `xlower::Array{Float}` : Array of lower x-coordinates
 - `xupper::Array{Float}` : Array of upper x-coordinates
 - `rlower::Array{Float}` : Array of lower r-coordinates
 - `rupper::Array{Float}` : Array of upper r-coordinates
"""
function split_wall(x, r)

    #find minimum x (LE)
    _, xminidx = findmin(x)

    #return split coordinates
    return x[1:xminidx], x[xminidx:end], r[1:xminidx], r[xminidx:end]
end
