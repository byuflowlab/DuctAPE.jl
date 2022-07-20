#=
Set up type drafts until better locations can be found
=#

# TODO: need to decide on units vs non-dimensional approach

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
"""
struct DuctGeometry{TA,TF}
    wallinnerxcoordinates::TA
    wallinnerrcoordinates::TA
    wallouterxcoordinates::TA
    wallouterrcoordinates::TA
    hubxcoordinates::TA
    hubrcoordinates::TA
    LEx::TF
    TEx::TF
    chord::TF
end

"""
    DuctGeometry(
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
    )
end

"""
    Airfoil{}

TODO: need to implement methods for allowing range of Re and Ma as well as AoA.

**Fields:**
 - `alpha::Array{Float}` : Array of angles of attack in degrees
 - `cl::Array{Float}` : Array of lift coefficients at each angle of attack
 - `cd::Array{Float}` : Array of drag coefficients at each angle of attack
"""
struct Airfoil{TA}
    alpha::TA
    cl::TA
    cd::TA
end

"""
    Rotor{TF,TI,TA,TAF}

TODO: need to add any other geometry that might be needed for rotors.

**Fields:**
 - `xlocation::Float` : x location of rotor plane
 - `numblades::Int` : number of rotor blades
 - `radialstations::Array{Float}` : array of radial stations defining rotor blade
 - `stationchords::Array{Float}` : array of chord lengths at radial stations defining rotor blade
 - `stationtwists::Array{Float}` : array of twist values (in degrees) at radial stations defining rotor blade
 - `stationairfoils::Array{Airfoil}` : array of airfoil data objects at radial stations defining rotor blade
"""
struct Rotor{TF,TI,TA,TAF}
    xlocation::TF
    numblades::TI
    radialstations::TA
    stationchords::TA
    stationtwists::TA
    stationairfoils::TAF
end

"""
    GridOptions{TF,TI}

**Fields:**
 - `num_radial_stations::Integer` : Number of radial stations (equal to number of rotor blade elements used in analysis)
 - `inlet_length::Float` : inlet length (unused)
 - `wake_length::Float` : length of wake behind duct relative to chord length
 - `wake_expansion_factor::Float` : expansion factor to apply to wake grid generation
"""
struct GridOptions{TF,TI}
    num_radial_stations::TI
    inlet_length::TF
    wake_length::TF
    wake_expansion_factor::TF
    # wakerelax::TB
end

"""
    defineGridOptions(
        num_radial_stations;
        inlet_length=0.5,
        wake_length=1.0,
        wake_expansion_factor=1.1
    )

Constructor function for the GridOptions object.

**Required Argument:**
 - `num_radial_stations::Integer` : Number of radial stations (equal to number of rotor blade elements used in analysis)

**Keyword Arguments:**
 - `inlet_length::Float` : inlet length (unused)
 - `wake_length::Float` : length of wake behind duct in terms of chord length
 - `wake_expansion_factor::Float` : expansion factor to apply to wake grid generation
"""
function defineGridOptions(
    num_radial_stations; inlet_length=0.5, wake_length=1.0, wake_expansion_factor=1.1
)
    return GridOptions(
        num_radial_stations, inlet_length, wake_length, wake_expansion_factor
    )
end

"""
    Grid{TF,TI}

Grid Object

**Fields:**
 - `x_grid_points::Matrix{Float}` : 2D Array of x grid points
 - `r_grid_points::Matrix{Float}` : 2D Array of radial grid points
 - `nx::Int` : number of x stations
 - `nr::Int` : number of radial stations
"""
struct Grid{TF,TI}
    x_grid_points::TF
    r_grid_points::TF
    nx::TI
    nr::TI
end

# """
#     OperatingCondition{TVI,TVR,TR,TF}

# **Fields:**
#  - `vinf::Array{Float}` : Array of freestream velocities
#  - `vref::Array{Float}` : Array of reference velocities
#  - `rpms::Array{Float}` : array of RPMs for rotor(s)
#  - `rho::Float` : air density value
#  - `vso::Float` : speed of sound value
#  - `mu::Float` : air viscosity value
# """
# struct OperatingCondition{TVI,TVR,TR,TF}
#     vinf::TVI
#     vref::TVR
#     rpms::TR
#     rho::TF
#     vso::TF
#     mu::TF
#     # viscous::TB
# end

# """
# """
# struct Output{TF} end

# """
# """
# struct ActuatorDisk{TF}
#     radialstations::TA
#     stationgammas::TA
# end

# """
# """
# struct Drag{TC,TA}
#     coordinates::TC
#     areas::TA
# end
