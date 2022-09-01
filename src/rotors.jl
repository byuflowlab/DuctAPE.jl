#=
Types and Functions related to rotors

Authors: Judd Mehr,
=#

###############################
##### ----- EXPORTS ----- #####
###############################

## -- TYPES

export RotorGeometry, BladeDimensions

## -- FUNCTIONS

export initialize_rotor_geometry, initialize_blade_dimensions

#######################################
##### ----- COMPOSITE TYPES ----- #####
#######################################

"""
    RotorGeometry{TF, TI, TA, TC, TT, TAF, TR, TM}

**Fields:**
 - `xlocation::Float` : x location of rotor plane, non-dimensional based on duct chord (max TE location - min LE location of hub/wall)
 - `numblades::Int` : number of rotor blades
 - `radialstations::Array{Float}` : array of radial stations defining rotor blade, non-dimensional with hub=0, tip=1
 - `tipgap::Float` : gap between blade tip and duct wall (not implemented yet)
 - `chords::Array{Float}` : array of chord lengths at radial stations defining rotor blade, non-dimensional based on blade tip radius
 - `twists::Array{Float}` : array of twist values (in degrees) at radial stations defining rotor blade
 - `skews::Array{Float}` : array of skew values (similar to sweep) at radial stations defining rotor blade, non-dimensional based on rotor tip radius. (note: this is for reference only, the solver can't use this information)
 - `rakes::Array{Float}` : array of rake values (similar to dihedral) at radial stations defining rotor blade, non-dimensional based on rotor tip radius. (note: this is for reference only right now. it may be implemented into the grid initialization functions later.)
 - `airfoils::Array{Airfoil}` : array of airfoil data objects at radial stations defining rotor blade
 - `solidity:Array{Float}` : array of rotor solidity at radial stations defining rotor blade, chord/distance between blade sections
 - `rpm::Float` : RPM of rotor
"""
struct RotorGeometry{TF,TI,TR,TG,TC,TT,TSk,TRa,TAF,TSo,TRpm}
    xlocation::TF
    numblades::TI
    radialstations::TR
    tipgap::TG
    chords::TC
    twists::TT
    skews::TSk
    rakes::TRa
    airfoils::TAF
    solidities::TSo
    rpm::TRpm
end

"""
    BladeDimensions{TF, TA}

**Fields:**
 - `rhub::Float` : hub radius (dimensional)
 - `rtip::Float` : tip radius (dimensional)
 - `rdim::Array{Float}` : array of dimensional radial stations
 - `cdim::Array{Float}` : array of dimensional chords
 - `tdim::Array{Float}` : array of twists (already dimensional in rotorgeometry)
 - `sweptannulus::Float` : area of blade swept annulus
 - `sweptarea::Float` : area of blade tip swept disk
"""
struct BladeDimensions{TF,TR,TC,TT}
    hubr::TF
    tipr::TF
    rdim::TR
    cdim::TC
    tdim::TT
    sweptannulus::TF
    sweptarea::TF
end

"""
    BladeAero{TRe,TMa,TCl,TCd,TM,TG,TW,TVa,TVt}

Blade aerodynamic values.
TODO: unsure where to put these, or if they even need to be put anywhere.

**Fields:**
 - `reynolds::Array{Float}` : local section reynolds numbers
 - `mach::Array{Float}` : local section mach numbers
 - `cl::Array{Float}` : local section coefficients of lift
 - `cd::Array{Float}` : local section coefficients of drag
 - `cm::Array{Float}` : local section coefficients of moment
 - `Gamma::Array{Float}` : local section circulations
 - `W::Array{Float}` : local inflow velocity
 - `vax::Array{Float}` : local axial velocity
 - `vtan::Array{Float}` : local tangential velocity
"""
struct BladeAero{TRe,TMa,TCl,TCd,TM,TG,TW,TVa,TVt}
    reynolds::TRe
    mach::TMa
    cl::TCl
    cd::TCd
    cm::TM
    Gamma::TG
    W::TW
    vax::TVa
    vtan::TVt
end

"""
    RotorVelocities{TA}

**Fields:**
 - `induced_axial_velocities::Array{Float}` : local section induced axial velocities
 - `induced_radial_velocities::Array{Float}` : local section induced radial velocities
 - `induced_tangential_velocities::Array{Float}` : local section induced tangential (circumferential) velocities
 - `absolute_axial_velocities::Array{Float}` : local section absolute axial velocities
 - `absolute_radial_velocities::Array{Float}` : local section absolute radial velocities
 - `absolute_tangential_velocities::Array{Float}` : local section absolute tangential velocities
 - `relative_axial_velocities::Array{Float}` : local section relative axial velocities
 - `relative_radial_velocities::Array{Float}` : local section relative radial velocities
 - `relative_tangential_velocities::Array{Float}` : local section relateive tangential velocities
 - `circulation::Array{Float}` : local section circulations
"""
struct RotorVelocities{TA}
    induced_axial_velocities::TA
    induced_radial_velocities::TA
    induced_tangential_velocities::TA
    absolute_axial_velocities::TA
    absolute_radial_velocities::TA
    absolute_tangential_velocities::TA
    relative_axial_velocities::TA
    relative_radial_velocities::TA
    relative_tangential_velocities::TA
    circulation::TA
end

###############################################
##### ----- GEOMETRY INITIALIZATION ----- #####
###############################################

"""
    function initialize_rotor_geometry(
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
        solidities=nothing,
    )

Initialize non-dimensional rotor geometry and other relative information.

if rotor rake is present, need to redo parts of grid initialization (and check that the unused portions of the code work now) to account for different x-locations of start of wake.  NOTE: not sure if rake will work with this solver. Need to think about that more before implementing.

**Arguments:**
 - `xlocation::Float` : x location of rotor relative to duct chord
 - `numblades::Float` : number of rotor blades
 - `numstations::Float` : number of radial stations (the length of the below sectional properties)
 - `chords::Array{Float}` : array of section chord lengths (relative to blade tip radius)
 - `twists::Array{Float}` : array of section twists (90 degrees is aligned with the axial direction) in degrees
 - `airfoils::Array{AFType}` : airfoil objects at each section
 - `rpm::Float` : RPM of rotor

**Keyword Arguments:**
 - `radialstations::Array{Float}` : non-dimensional radial stations (0 = hub radius, 1 = tip radius)
 - `tipgap::Float` : non-dimensional relative to blade length (Not implemented yet)
 - `skews::Array{Float}` : sectional skew distance relative to tip radius (not implemented yet)
 - `rakes::Array{Float}` : sectional rake distances relative to tip radius (not implemented yet)
 - `solidities::Array{Float} or Nothing` : sectional rotor solidities, if using.

**Returns:**
 - `rotorgeometry::DuctTAPE.RotorGeometry` : Rotor geometry object
"""
function initialize_rotor_geometry(
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
    solidities=nothing,
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

    #if skews aren't set, set them to zeros
    if skews == nothing
        skews = [0.0 for i in 1:numstations]
    end

    #if rakes aren't set, set them to zeros
    if rakes == nothing
        rakes = [0.0 for i in 1:numstations]
    end

    # if solidities is nothing, set to -1's to flag airfoil functions later.
    if solidities == nothing
        solidities = [-1.0 for i in 1:numstations]
    end

    #Define Rotor Object
    return RotorGeometry(
        xlocation,
        numblades,
        radialstations,
        tipgap,
        chords,
        twists,
        skews,
        rakes,
        airfoils,
        solidities,
        rpm,
    )
end

"""
    reinterpolate_rotor!(wakegrid, rotor, rotoridx)

Since the wake grid relaxes and is not aligned with aft rotor radial stations, this function reinterpolates rotor data based on updated radial stations.

(The `rotor` inputs is the only one updated by this function.)

**Arguments:**
 - `wakegrid::DuctTAPE.WakeGridGeometry` : wake grid geometry object
 - `rotor::DuctTAPE.RotorGeometry` : the rotor geometry to update
 - `rotoridx::Int` : index in the x direction for where the rotor lies on the wake grid
"""
function reinterpolate_rotor!(wakegrid, rotor, rotoridx)

    #rename things for convenience
    gridxs = wakegrid.x_grid_points
    gridrs = wakegrid.r_grid_points
    nr = length(rotor.radialstations)
    new_rad_stash = gridxs[rotoridx, :]

    ## -- Calculate New Section Properties -- ##

    #TODO: none of this will be nothing at this point (probably), can probably simplify and clean up this function.

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
    if rotor.skews != nothing
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

    # update solidity
    if rotor.solidities != nothing
        new_solidities = FLOWMath.Akima(rotor.radialstations, rotor.solidities)
        for i in 1:nr
            rotor.solidities[i] = new_solidities(new_rad_stash[i])
        end
    end

    # update radial stations at the end after using both old and new for splines
    for i in 1:nr
        rotor.radialstations[i] = new_rad_stash[i]
    end

    return nothing
end

"""
    initialize_blade_dimensions(ductsplines, Rotor)

Initilialize needed blade information for various calculations during the solution process.

**Arguments:**
 - `ductsplines::DuctTAPE.ductsplines` : ductsplines object containing splines for duct wall and hub
 - `Rotor::DuctTAPE.Rotor` : Rotor object for which to define blade information
"""
function initialize_blade_dimensions(ductgeometry, ductsplines, Rotor)

    #unpack splines for convenience
    wallspline = ductsplines.wallinnerspline
    hubspline = ductsplines.hubspline

    #find r-position of wall and hub at rotor location
    rhub = max(0.0, hubspline(Rotor.xlocation * ductgeometry.chord))
    rtip = wallspline(Rotor.xlocation * ductgeometry.chord)

    #get dimensional radial station locations
    #use linear transform to go from non-dim to dimensional range
    rdim = lintran(
        rhub, rtip, Rotor.radialstations[1], Rotor.radialstations[end], Rotor.radialstations
    )

    # dimensional chords
    cdim = Rotor.chords .* rtip

    #twists are unchanged
    tdim = Rotor.twists

    # check that the rotor radial positions are all positive and that the array is monotonic to avoid issues including hanging in wake grid generation
    @assert ismonotonic(rdim) && all(>=(0), rdim)

    #calculate swept area
    sweptannulus = pi * (rtip^2 - rhub^2)
    sweptarea = pi * rtip^2

    return BladeDimensions(rhub, rtip, rdim, cdim, tdim, sweptannulus, sweptarea)
end

"""
    set_rotor_velocities(
        vax, vrad, vtan, vinf, omega, radialstations, wwa=0.0, wwt=0.0, vfac=1.0
    )

Define RotorVelocities object.

**Arguments:**
 - `vax::Array{Float}` : axial velocities
 - `vrad::Array{Float}` : radial velocities
 - `vtan::Array{Float}` : tangential velocities
 - `vinf::Array{Float}` : freestream velocity
 - `omega::Float` : rotation rate in rad/s
 - `wwa::Array{Float} : user defined additional axial velocity (unused)
 - `wwt::Array{Float} : user defined additional tangential velocity (unused)
 - `vfac::Float` : user defined velocity factor (unused)

**Returns:**
 - `rotorvelocities::DuctTAPE.RotorVelocities` : induced, absolute, and relative velocities along rotor blades.
"""
function set_rotor_velocities(
    vax, vrad, vtan, vinf, omega, radialstations, wwa=0.0, wwt=0.0, vfac=1.0
)

    #notes: in dfdc, the following indicies are defined as
    # 1 = axial
    # 2 = radial
    # 3 = tangential

    # get number of stations for convenience
    numstations = length(radialstations)

    for i in 1:numstations

        ## -- Absolute Velocities

        #tangential
        vtanabs[i] = vtan[i] + wwt

        #axial
        vaxabs[i] = (vax[i] + vinf + wwa) * vfac #TODO what is a BB entry? seems vfac is defined there.

        #radial
        vradabs[i] = vrad[i]

        #flow angle?
        vma = sqrt(vax[i]^2 + vrad[i]^2)
        vva = sqrt(vma^2 + vtan[i]^2)
        if vtan[i] != 0.0
            phi_absolute = atan(vma, -vtan[i])
        else
            phi_absolute = 0.5 * pi
        end

        ## -- Relative Velocities

        #tangential
        vtanrel[i] = vtanabs[i] - omega * radialstations[i]

        #axial
        vaxrel[i] = vaxabs[i]

        #radial
        vradrel[i] = vradabs[i]

        #flow angle?
        vmr = sqrt(vaxrel[i]^2 + vradrel[i]^2)
        vvr = sqrt(vmr^2 + vtanrel[i]^2)
        if vtanrel[i] != 0.0
            phi_relative = atan2(vmr, -vtanrel[i])
        else
            phi_relative = 0.5 * pi
        end
    end

    return RotorVelocities(
        vax, vrad, vtan, vaxabs, vradabs, vtanabs, vaxrel, vradrel, vtanrel
    )
end

##################################
##### ----- FOR SOLVER ----- #####
##################################

"""
    blade_section_gamma(W, chord, cl)

Calculate section circulation.

**Arguments:**
 - `W::Float` : Inflow velocity
 - `chord::Float` : section chord length
 - `cl::Float` : section lift coefficient

**Returns:**
 - `gamma::Float` : section circulation
"""
function blade_section_gamma(W, chord, cl)
    return 0.5 * W * chord * cl #eqn 75 in dfdc docs
end

