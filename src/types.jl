#=
Custom Composite Type Definitions
=#

"""
    RotorGeometry{TF, TI, TA, TC, TT, TAF, TR, TM}

**Fields:**
 - `xlocation::Float` : x location of rotor plane, non-dimensional based on duct chord (max TE location - min LE location of hub/wall)
 - `Rtip::Float` : blade tip radius (dimensional)
 - `numblades::Int` : number of rotor blades
 - `radialstations::Array{Float}` : array of radial stations defining rotor blade, non-dimensional with hub=0, tip=1
 - `chords::Array{Float}` : array of chord lengths at radial stations defining rotor blade, non-dimensional based on blade tip radius
 - `twists::Array{Float}` : array of twist values (in degrees) at radial stations defining rotor blade
 - `airfoils::Array{Airfoil}` : array of airfoil data objects at radial stations defining rotor blade
 - `rpm::Float` : RPM of rotor
"""
struct RotorGeometry{TF,TRt,TRh,TI,TR,TC,TT,TAF,TRpm}
    xlocation::TF #decide if this should be non-dim or not
    Rtip::TRt
    Rhub::TRh
    numblades::TI
    nref::TI #number of refinement stations
    radialstations::TR
    chords::TC
    twists::TT
    airfoils::TAF
    # solidities::TSo
    RPM::TRpm
end
# TODO: Add to rotor object later
#  - `tipgap::Float` : gap between blade tip and duct wall (not implemented yet)
#  - `solidity:Array{Float}` : array of rotor solidity at radial stations defining rotor blade, chord/distance between blade sections

"""
    Freestream{TVI,TVR,TF}

**Fields:**
 - `vinf::Float` : Freestream velocity
 - `rho::Float` : Air density value
 - `mu::Float` : Air dynamic viscosity value
 - `asound::Float` : Speed of sound value
"""
struct Freestream{TVI,TF}
    vinf::TVI
    rho::TF
    mu::TF
    asound::TF
end
