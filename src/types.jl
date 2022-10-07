#=
Custom Composite Type Definitions
=#

"""
    RotorGeometry{TF, TI, TA, TC, TT, TAF, TR, TM}

**Fields:**
 - `xlocation::Float` : x location of rotor plane, non-dimensional based on duct chord (max TE location - min LE location of hub/wall)
 - `numblades::Int` : number of rotor blades
 - `radialstations::Array{Float}` : array of radial stations defining rotor blade, non-dimensional with hub=0, tip=1
 - `chords::Array{Float}` : array of chord lengths at radial stations defining rotor blade, non-dimensional based on blade tip radius
 - `twists::Array{Float}` : array of twist values (in degrees) at radial stations defining rotor blade
 - `airfoils::Array{Airfoil}` : array of airfoil data objects at radial stations defining rotor blade
 - `rpm::Float` : RPM of rotor
"""
struct RotorGeometry{TF,TI,TR,TC,TT,TAF,TRpm}
    xlocation::TF #decide if this should be non-dim or not
    numblades::TI
    nref::TI #number of refinement stations
    radialstations::TR
    # tipgap::TG
    chords::TC
    twists::TT
    # skews::TSk
    # rakes::TRa
    airfoils::TAF
    # solidities::TSo
    RPM::TRpm
end
# TODO: Add to rotor object later
#  - `tipgap::Float` : gap between blade tip and duct wall (not implemented yet)
#  - `skews::Array{Float}` : array of skew values (similar to sweep) at radial stations defining rotor blade, non-dimensional based on rotor tip radius. (note: this is for reference only, the solver can't use this information)
#  - `rakes::Array{Float}` : array of rake values (similar to dihedral) at radial stations defining rotor blade, non-dimensional based on rotor tip radius. (note: this is for reference only right now. it may be implemented into the grid initialization functions later.)
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
