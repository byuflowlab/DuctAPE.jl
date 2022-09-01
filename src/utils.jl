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
    sinespace(n)

Calculates sine spaced points from 0 to pi/2 (viz. sin(theta), theta ∈ [0, pi/2])

**Arguments:**
 - `n::Int` : Number of points to calculate

**Returns::
 - `pts::Array{Float}` : Array of sine spaced points
"""
function sinespace(n)
    theta = range(0, pi / 2; length=n)
    return sin.(theta)
end

"""
    cosinespace(n)

Calculates cosine spaced points from 0 to pi/2 (viz. cosin(theta), theta ∈ [0, pi/2])

**Arguments:**
 - `n::Int` : Number of points to calculate

**Returns::
 - `pts::Array{Float}` : Array of cosine spaced points
"""
function cosinespace(n)
    theta = range(0, pi / 2; length=n)
    return cos.(theta)
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
