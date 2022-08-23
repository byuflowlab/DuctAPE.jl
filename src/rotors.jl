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
 - `chords::Array{Float}` : array of chord lengths at radial stations defining rotor blade
 - `twists::Array{Float}` : array of twist values (in degrees) at radial stations defining rotor blade
 - `skews::Array{Float}` : array of skew values (similar to sweep) at radial stations defining rotor blade
 - `rakes::Array{Float}` : array of rake values (similar to dihedral) at radial stations defining rotor blade
 - `airfoils::Array{Airfoil}` : array of airfoil data objects at radial stations defining rotor blade
 - `reynolds::Array{Float}` : array of reynolds numbers at radial stations defining rotor blade
 - `solidity:Array{Float}` : array of rotor solidity at radial stations defining rotor blade
 - `machs::Array{Float}` : array of machs numbers at radial stations defining rotor blade (currently unsupported).
"""
struct Rotor{TF,TI,TR,TG,TC,TT,TSk,TRa,TAF,TRe,TSo,TM,TRpm}
    xlocation::TF
    numblades::TI
    radialstations::TR
    tipgap::TG
    chords::TC
    twists::TT
    skews::TSk
    rakes::TRa
    airfoils::TAF
    reynolds::TRe
    solidities::TSo
    machs::TM #TODO: need to add compressibility corrections to solver before this can be used.
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
    solidities=nothing,
    machs=nothing,
)

    #Set up radial stations if not defined by user.
    if radialstations != nothing
        #make sure that the radial stations that are defined meet the requirements to not break stuff...
        @assert numstations == length(radialstations) &&
            ismonotonic(radialstations) &&
            all(>=(0), radialstations)
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

    if solidities == nothing
        solidities = [0.0 for i in 1:numstations]
    end

    if machs == nothing
        machs = [0.0 for i in 1:numstations]
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
        solidities,
        machs,
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
    hubspline = DuctSplines.hubspline

    #find r-position of wall and hub at rotor location
    rhub = max(0.0, hubspline(Rotor.xlocation))
    rtip = wallspline(Rotor.xlocation)

    #get dimensional radial station locations
    #use linear transform to go from non-dim to dimensional range
    rdim =
        rhub .+
        (rtip - rhub) / (Rotor.radialstations[end] - Rotor.radialstations[1]) .*
        (Rotor.radialstations .- Rotor.radialstations[1])

    # check that the rotor radial positions are all positive and that the array is monotonic to avoid issues including hanging in wake grid generation
    @assert ismonotonic(rdim) && all(>=(0), rdim)

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

"""
since wake relaxes and is not aligned with aft rotor radial stations, need to reinterpolate rotor based on updated radial stations
"""
function reinterpolate_rotor!(wakegrid, rotor, rotoridx)

    #rename things for convenience
    gridxs = wakegrid.x_grid_points
    gridrs = wakegrid.r_grid_points
    nr = length(rotor.radialstations)
    new_rad_stash = gridxs[rotoridx, :]

    ## -- Calculate New Section Properties -- ##

    # update chords
    if rotor.chords != nothing
        new_chords = FLOWMath.Akima(rotor.radialstations, rotor.chords)
        for i in 1:nr
            rotor.chords[i] = new_chords(new_rad_stash[i])
        end
    end

    # update twists
    if rotor.twists != nothing
        new_twists = FLOWMath.Akima(rotor.radialstations, rotor.twists)
        for i in 1:nr
            rotor.twists[i] = new_twists(new_rad_stash[i])
        end
    end

    # updates skews
    if rotor.twists != nothing
        new_skews = FLOWMath.Akima(rotor.radialstations, rotor.skews)
        for i in 1:nr
            rotor.skews[i] = new_skews(new_rad_stash[i])
        end
    end

    # update rakes
    if rotor.rakes != nothing
        new_rakes = FLOWMath.Akima(rotor.radialstations, rotor.rakes)
        for i in 1:nr
            rotor.rakes[i] = new_rakes(new_rad_stash[i])
        end
    end

    # update reynolds
    if rotor.reynolds != nothing
        new_reynolds = FLOWMath.Akima(rotor.radialstations, rotor.reynolds)
        for i in 1:nr
            rotor.reynolds[i] = new_reynolds(new_rad_stash[i])
        end
    end

    # update solidity
    if rotor.solidities != nothing
        new_solidities = FLOWMath.Akima(rotor.radialstations, rotor.solidities)
        for i in 1:nr
            rotor.solidities[i] = new_solidities(new_rad_stash[i])
        end
    end

    # update machs
    if rotor.machs != nothing
        new_machs = FLOWMath.Akima(rotor.radialstations, rotor.machs)
        for i in 1:nr
            rotor.machs[i] = new_mach(new_rad_stash[i])
        end
    end

    # update radial stations at the end after using both old and new for splines
    for i in 1:nr
        rotor.radialstations[i] = new_rad_stash[i]
    end

    return nothing
end

##################################
##### ----- FOR SOLVER ----- #####
##################################

"""
"""
function blade_section_gamma(W, chord, cl)
    return 0.5 * W * chord * cl #eqn 75 in dfdc docs
end

