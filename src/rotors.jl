#=
Types and Functions related to rotors

Authors: Judd Mehr,
=#

#############################
##### ----- TYPES ----- #####
#############################

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
struct Rotor{TF,TI,TR,TC,TT,TS,TRa,TAF,TRe,TCa,TM,TRpm}
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
    rpm::TRpm
end

"""
will probably need a blade obeject to hold specific dimensional information concerning the rotor blades.  possibly will put some of the rotor object fields here, and then have a blade field in the rotor object instead.
"""
struct Blade{} end

######################################
##### ----- INITIALIZATION ----- #####
######################################

"""
see dfdcsubs.f line 1257 for a good idea of where to start.

Need to place rotor inside duct, find dimensional hub and tip radius based on hub and wall radial positions at rotor x position (also consider adding in user specified tip gap parameter as well)

also get dimensional rotor blade section locations. These should inform the wake grid generation function.

if rotor rake is present, need to redo parts of grid initialization (and check that the unused portions of the code work now) to account for different x-locations of start of wake.
"""
function initialize_rotor_geometry() end

"""
need to set up source and vortex lines representing starting points of rotor wakes.
"""
function initialize_rotor_wake() end

"""
"""
function initialize_rotor()
    #call the other initialize functions so you just call one function later.
end

##################################
##### ----- FOR SOLVER ----- #####
##################################

"""
"""
function blade_section_gamma(W, chord, cl)
    return 0.5 * W * chord * cl #eqn 75 in dfdc docs
end

