#=
Types and Functions pertaining to duct walls (hub, duct)

Authors: Judd Mehr,
=#

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
        bluntTEtol=1e-9,
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
 - `chord::Float` : manually defined chord length.  Set to difference between leading and trailing edges otherwise.
 - `bluntTEtol::Float` : tolerance for how close trailing edge points need to be before being considered a blunt trailing edge.

# Notes
The wall inner and outer and hub coordinates are all saved from leading to trailing edge (this helps with splining things later).
It is assumed that the wall is similar to an airfoil that may or may not have a blunt trailing edge.
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

    #flip things as needed:
    if wallinnerxcoordinates[1] > wallinnerxcoordinates[end]
        wallinnerxcoordinates = reverse(wallinnerxcoordinates)
        wallinnerrcoordinates = reverse(wallinnerrcoordinates)
    end

    if wallouterxcoordinates[1] > wallouterxcoordinates[end]
        wallouterxcoordinates = reverse(wallouterxcoordinates)
        wallouterrcoordinates = reverse(wallouterrcoordinates)
    end

    #if not hub coordinates are provided match the x coordinates with the wall
    if hubxcoordinates == nothing
        if hubrcoordinates != nothing
            @warn("no x coordinates defined for hub, overwriting r coordinates to zero")
        end
        hubxcoordinates = wallinnerxcoordinates
        hubrcoordinates = [0.0 for i in 1:length(wallinnerxcoordinates)]
    end

    #add zeros for hub r coordinates if nothing is provided
    if hubrcoordinates == nothing
        @warn("no r coordinates defined for hub, setting r coordinates to zero")
        hubrcoordinates = [0.0 for i in 1:length(hubxcoordinates)]
    end

    #check if blunt trailing edges
    wallTEdist = sqrt(
        (wallinnerxcoordinates[end] - wallouterxcoordinates[end])^2 +
        (wallinnerrcoordinates[end] - wallouterrcoordinates[end])^2,
    )
    hubTEdist = hubrcoordinates[end]

    #set blunt trailing edge flags
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

    # create wall spline fields
    wallinnerspline = FLOWMath.Akima(wallinnerxcoordinates, wallinnerrcoordinates)
    wallouterspline = FLOWMath.Akima(wallouterxcoordinates, wallouterrcoordinates)
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
    DuctSplines(wallinnerspline, wallouterspline, hubspline)
end
