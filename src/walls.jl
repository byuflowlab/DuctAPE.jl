#=
Types and Functions pertaining to duct walls (hub, duct)

Authors: Judd Mehr,
=#

"""
    DuctGeometry{TA,TF}

**Fields:**
 - `wallinnerxcoordinates::Array{Float}` : x coordinates of inner (lower) wall geometry
 - `wallouterxcoordinates::Array{Float}` : x coordinates of outer (upper) wall geometry
 - `wallinnerrcoordinates::Array{Float}` : r coordinates of inner (lower) wall geometry
 - `wallouterrcoordinates::Array{Float}` : r coordinates of outer (upper) wall geometry
 - `hubxcoordinates::Array{Float}` : x coordinates of hub geometry
 - `hubrcoordinates::Array{Float}` : r coordinates of hub geometry
 - `LEx::Float` : x-position of leading edge
 - `TEx::Float` : x-position of trailing edge
 - `chord::Float` : chord length
 - `wallbluntTE::Bool` : flag for blunt trailing edge on wall
 - `hubbluntTE::Bool` : flag for blunt trailing edge on hub
"""
struct DuctGeometry{TA,TF,TB}
    wallinnerxcoordinates::TA
    wallinnerrcoordinates::TA
    wallouterxcoordinates::TA
    wallouterrcoordinates::TA
    hubxcoordinates::TA
    hubrcoordinates::TA
    LEx::TF
    TEx::TF
    chord::TF
    wallbluntTE::TB
    hubbluntTE::TB
end

"""
"""
struct DuctSplines{TSD,TSH}
    wallinnerspline::TSD
    hubspline::TSH
end

"""
    defineDuctGeometry(
        wallinnerxcoordinates,
        wallinnerrcoordinates,
        wallouterxcoordinates,
        wallouterrcoordinates,
        hubxcoordinates=nothing,
        hubrcoordinates=nothing;
        LEx=nothing,
        TEx=nothing,
        chord=nothing,
    )

Constructor function for the DuctGeometry Object.

**Required Arguments:**
 - `wallinnerxcoordinates::Array{Float}` : x coordinates of inner (lower) wall geometry
 - `wallouterxcoordinates::Array{Float}` : x coordinates of outer (upper) wall geometry
 - `wallinnerrcoordinates::Array{Float}` : r coordinates of inner (lower) wall geometry
 - `wallouterrcoordinates::Array{Float}` : r coordinates of outer (upper) wall geometry

**Optional Arguments:**
 - `hubxcoordinates::Array{Float}` : x coordinates of hub geometry
 - `hubrcoordinates::Array{Float}` : r coordinates of hub geometry
Note, if hub x and r coordinates are not set, the x coordinates for the inner wall will be used and the r coordinates will be set to zero.  Also note that if one of these is unset, the other must also be unset.

**Keyword Arguments:**
 - `LEx::Float` : x-position of manually defined leading edge.  Set to foremost x-coordinate of duct and hub geometry otherwise.
 - `TEx::Float` : x-position of mannually defined trailing edge.  Set to the rear-most x-coordinate of duct and hub geometry otherwise.
 - `chord::Float` : manuall defined chord length.  Set to difference between leading and trailing edges otherwise.
"""
function defineDuctGeometry(
    wallinnerxcoordinates,
    wallinnerrcoordinates,
    wallouterxcoordinates,
    wallouterrcoordinates,
    hubxcoordinates=nothing,
    hubrcoordinates=nothing;
    LEx=nothing,
    TEx=nothing,
    chord=nothing,
    bluntTEtol=1e-9,
)
    if TEx == nothing
        TEx = maximum([hubxcoordinates; wallinnerxcoordinates])
    end

    if LEx == nothing
        LEx = minimum([hubxcoordinates; wallinnerxcoordinates])
    end

    if chord == nothing
        chord = TEx - LEx
    end

    #Check that geometry will work!
    @assert chord >= 0.0 "Chord length must be greater than zero."
    if maximum(wallinnerxcoordinates) <= minimum(hubxcoordinates) ||
        minimum(wallinnerxcoordinates) >= maximum(hubxcoordinates)
        @warn(
            "there is no overlap between the hub and wall.  functionality for this configuration is not supported."
        ) #TODO: originally had this as a warning, but the wake caclulation currently hangs in an infinite loop if there is no overlap.  This could probably be fixed, but also probably doesn't need to be if optimization constraints are placed properly
    end

    if hubxcoordinates == nothing
        if hubrcoordinates != nothing
            @warn("no x coordinates defined for hub, overwriting r coordinates to zero")
        end
        hubxcoordinates = wallinnerxcoordinates
        hubrcoordinates = [0.0 for i in 1:length(wallinnerxcoordinates)]
    end

    if hubrcoordinates == nothing
        @warn("no r coordinates defined for hub, setting r coordinates to zero")
        hubrcoordinates = [0.0 for i in 1:length(hubxcoordinates)]
    end

    #check of blunt trailing edges
    wallTEdist = sqrt(
        (wallinnerxcoordinates[end] - wallouterxcoordinates[1])^2 +
        (wallinnerrcoordinates[end] - wallouterrcoordinates[1])^2,
    )
    hubTEdist = hubrcoordinates[end]

    if wallTEdist >= bluntTEtol
        wallbluntTE = true
    else
        wallbluntTE = false
    end

    if hubTEdist >= bluntTEtol
        hubbluntTE = true
    else
        hubbluntTE = false
    end

    #TODO: create wall spline fields
    wallinnerspline = FLOWMath.Akima(wallinnerxcoordinates, wallinnerrcoordinates)
    hubspline = FLOWMath.Akima(hubxcoordinates, hubrcoordinates)

    return DuctGeometry(
        wallinnerxcoordinates,
        wallinnerrcoordinates,
        wallouterxcoordinates,
        wallouterrcoordinates,
        hubxcoordinates,
        hubrcoordinates,
        LEx,
        TEx,
        chord,
        wallbluntTE,
        hubbluntTE,
    ),
    DuctSplines(wallinnerspline, hubspline)
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
