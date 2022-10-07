"""
    ismonotonic(A)

Checks if array, A is monotonic

**Arguments:**
 - `A::Array` : Array in question

**Returns:**
 - `true::Bool` : true or false.
"""
function ismonotonic(A)
    current = A[1]
    for i in 2:length(A)
        newval = A[i]
        isless(newval, current) && return false
        current = newval
    end
    return true
end

"""
    lintran(rb1, rbend, ra1, raend, ra)

Linear transfrom of ra from range (ra1, raend) to (rb1, rbend)

**Arguments:**
 - `rb1::Float` : lower bound of range b
 - `rbend::Float` : upper bound of range b
 - `ra1::Float` : lower bound of range a
 - `raend::Float` : upper bound of range a
 - `ra::Array{Float}` : array of values on range a

**Returns:**
 - `rb::Array{Float}` : array of transformed values on range b
"""
function lintran(rb1, rbend, ra1, raend, ra)
    return rb1 .+ (rbend - rb1) / (raend - ra1) .* (ra .- ra1)
end

"""
    get_omega(rpm)

Calculates rad/s from RPM

**Arguments:**
 - `rpm::Float` : RPM

**Returns:**
 - `omega::Float` : radians per second
"""
function get_omega(rpm)
    return rpm * pi / 30.0
end

"""
    split_wall(x,r)

Splits full airfoil coordinates into upper and lower halves.

Only works based on geometry.  Splits at lowest x-value, does not split based on stagnation point.

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
