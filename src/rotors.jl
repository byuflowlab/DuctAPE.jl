#=
Types and Functions related to rotors
=#


"""
    Rotor{TF, TI, TA, TC, TT, TAF, TR, TM}

TODO: need to add any other geometry that might be needed for rotors.

**Fields:**
 - `xlocation::Float` : x location of rotor plane
 - `numblades::Int` : number of rotor blades
 - `radialstations::Array{Float}` : array of radial stations defining rotor blade, non-dimensional with hub=0, tip=1
 - `stationchords::Array{Float}` : array of chord lengths at radial stations defining rotor blade
 - `stationtwists::Array{Float}` : array of twist values (in degrees) at radial stations defining rotor blade
 - `stationskews::Array{Float}` : array of skew values (similar to sweep) at radial stations defining rotor blade
 - `stationrakes::Array{Float}` : array of rake values (similar to dihedral) at radial stations defining rotor blade
 - `stationairfoils::Array{Airfoil}` : array of airfoil data objects at radial stations defining rotor blade
 - `stationreynolds::Array{Float}` : array of reynolds numbers at radial stations defining rotor blade
 - `stationcascade::Array{Float}` : array of cascade parameter (TODO, need to find out what this is) at radial stations defining rotor blade
 - `stationmach::Array{Float}` : array of mach numbers at radial stations defining rotor blade (currently unsupported).
"""
struct Rotor{TF,TI,TR,TC,TT,TS,TRa,TAF,TRe,TCa,TM}
    xlocation::TF
    numblades::TI
    radialstations::TR
    stationchords::TC
    stationtwists::TT
    stationskews::TS
    stationrakes::TRa
    stationairfoils::TAF
    stationreynolds::TRe
    stationcascade::TCa
    stationmach::TM #TODO: need to add compressibility corrections to solver before this can be used.
end
