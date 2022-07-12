#=
Set up type drafts until better locations can be found
=#

# TODO: need to decide on units vs non-dimensional approach

"""
    Duct{TM}

**Fields:**
 - `wallinnerxcoordinates::Array{Float} : x coordinates of inner (lower) wall geometry
 - `wallouterxcoordinates::Array{Float} : x coordinates of outer (upper) wall geometry
 - `wallinnerrcoordinates::Array{Float} : r coordinates of inner (lower) wall geometry
 - `wallouterrcoordinates::Array{Float} : r coordinates of outer (upper) wall geometry
 - `hubxcoordinates::Array{Float} : x coordinates of hub geometry
 - `hubrcoordinates::Array{Float} : r coordinates of hub geometry
"""
struct Duct{TA,TF}
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

function Duct(
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
        LEx = maximum([hubxcoordinates; wallinnerxcoordinates])
    end

    if chord == nothing
        chord = TEx - LEx
    end

    if hubxcoordinates == nothing
        @assert hubrcoordinates == nothing
        hubxcoordinates = wallinnerxcoordinates
        hubrcoordinates = [0.0 for i in 1:length(wallinnerxcoordinates)]
    end

    return Duct(
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
    OperatingCondition{TVI,TVR,TR,TF}

**Fields:**
 - `vinf::Array{Float}` : Array of freestream velocities
 - `vref::Array{Float}` : Array of reference velocities
 - `rpms::Array{Float}` : array of RPMs for rotor(s)
 - `rho::Float` : air density value
 - `vso::Float` : speed of sound value
 - `mu::Float` : air viscosity value
"""
struct OperatingCondition{TVI,TVR,TR,TF}
    vinf::TVI
    vref::TVR
    rpms::TR
    rho::TF
    vso::TF
    mu::TF
    # viscous::TB
end

"""
    GridOptions{TF,TI,TM}

**Fields:**
 - `wake_length::Float` : length of wake behind duct
"""
struct GridOptions{TF,TI}
    num_xinlet_stations::TI
    num_xduct_stations::TI
    num_xwake_stations::TI
    num_radial_stations::TI
    inlet_length::TF
    wake_length::TF
    # wakerelax::TB
end

function GridOptions(
    num_xinlet_stations,
    num_xduct_stations,
    num_xwake_stations,
    num_radial_stations;
    inlet_length=0.5,
    wake_length=1.0,
)
    return GridOptions(
        num_xinlet_stations,
        num_xduct_stations,
        num_xwake_stations,
        num_radial_stations,
        inlet_length,
        wake_length,
    )
end

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
