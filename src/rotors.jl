#=
Types and Functions related to rotors

Authors: Judd Mehr,
=#

#############################
##### ----- TYPES ----- #####
#############################

"""
    Rotor{TF, TI, TA, TC, TT, TAF, TR, TM}

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
struct Rotor{TF,TI,TR,TG,TC,TT,TS,TRa,TAF,TRe,TCa,TM,TRpm}
    xlocation::TF
    numblades::TI
    radialstations::TR
    tipgap::TG
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
    Blade{TF, TA}

**Fields:**
 - `rhub::Float` : hub radius (dimensional)
 - `rtip::Float` : tip radius (dimensional)
 - `rdim::Array{Float}` : array of dimensional radial stations
 - `sweptannulus::Float` : area of blade swept annulus
 - `sweptarea::Float` : area of blade tip swept disk
"""
struct Blade{TF,TA}
    hubr::TF
    tipr::TF
    rdim::TA
    sweptannulus::TF
    sweptarea::TF
end

######################################
##### ----- INITIALIZATION ----- #####
######################################

"""

also get dimensional rotor blade section locations. These should inform the wake grid generation function. NOTE: definitely need to start wake generation with rotor locations as well as LE/TE location of hub/duct.  Grid should line up with each of those.

if rotor rake is present, need to redo parts of grid initialization (and check that the unused portions of the code work now) to account for different x-locations of start of wake.  NOTE: not sure if rake will work with this solver. Need to think about that more before implementing.
"""
function initialize_rotor(
    xlocation,
    numblades,
    numstations,
    chords,
    twists,
    airfoils,
    rpm;
    radialstations=nothing,
    tipgap=0.0, #non-dimensional relative to blade length
    skews=nothing,
    rakes=nothing,
    reynolds=nothing,
    cascade=nothing,
    mach=nothing,
)

    #Set up radial stations if not defined by user.
    if radialstations != nothing
        @assert numstations == length(radialstations)
    else
        #use linear spacing. user can define custom spacing otherwise
        radialstations = collect(range(0.0, 1.0 - tipgap; length=nstations))
    end

    #if none of the following are set, set them to zeros
    if skews == nothing
        skews = [0.0 for i in 1:numstations]
    end

    if rakes == nothing
        rakes = [0.0 for i in 1:numstations]
    end

    if reynolds == nothing
        reynolds = [0.0 for i in 1:numstations]
    end

    if cascade == nothing
        cascade = [0.0 for i in 1:numstations]
    end

    if mach == nothing
        mach = [0.0 for i in 1:numstations]
    end

    #Define Rotor Object
    return Rotor(
        xlocation,
        numblades,
        radialstations,
        tipgap,
        chords,
        twists,
        skews,
        rakes,
        airfoils,
        reynolds,
        cascade,
        mach,
        rpm,
    )
end

"""
    initialize_blade(DuctSplines, Rotor)

Initilialize needed blade information for various calculations during the solution process.

**Arguments:**
 - `DuctSplines::DuctTAPE.DuctSplines` : DuctSplines object containing splines for duct wall and hub
 - `Rotor::DuctTAPE.Rotor` : Rotor object for which to define blade information
"""
function initialize_blade(DuctSplines, Rotor)

    #unpack splines for convenience
    wallspline = DuctSplines.wallinnerspline
    hubspline = DuctSPlines.hubspline

    #find r-position of wall and hub at rotor location
    rhub = hubspline(Rotor.xlocation)
    rtip = wallspline(Rotor.xlocation)

    #get dimensional radial station locations
    #use linear transform to go from non-dim to dimensional range
    rdim =
        rhub .+
        (rhub - rtip) / (Rotor.radialstations[2] - Rotor.radialstations[1]) .*
        (Rotor.radialstations .- Rotor.radialstations[1])

    #calculate swept area
    sweptannulus = pi * (rtip^2 - rhub^2)
    sweptarea = pi * rtip^2

    return Blade(rhub, rtip, rdim, sweptannulus, sweptarea)
end

"""
need to set up source and vortex lines representing starting points of rotor wakes.
have this function called by a more global function that holds all the source and vortex lines.
"""
function initialize_rotor_wake() end

##################################
##### ----- FOR SOLVER ----- #####
##################################

"""
"""
function blade_section_gamma(W, chord, cl)
    return 0.5 * W * chord * cl #eqn 75 in dfdc docs
end

