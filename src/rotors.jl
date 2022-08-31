#=
Types and Functions related to rotors

Authors: Judd Mehr,
=#

#############################
##### ----- TYPES ----- #####
#############################

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
end

######################################
##### ----- INITIALIZATION ----- #####
######################################

"""

also get dimensional rotor blade section locations. These should inform the wake grid generation function. NOTE: definitely need to start wake generation with rotor locations as well as LE/TE location of hub/duct.  Grid should line up with each of those.

if rotor rake is present, need to redo parts of grid initialization (and check that the unused portions of the code work now) to account for different x-locations of start of wake.  NOTE: not sure if rake will work with this solver. Need to think about that more before implementing.
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

    #if none of the following are set, set them to zeros
    if skews == nothing
        skews = [0.0 for i in 1:numstations]
    end

    if rakes == nothing
        rakes = [0.0 for i in 1:numstations]
    end

    if solidities == nothing
        solidities = [0.0 for i in 1:numstations]
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

#FOR initialize_blade_aero: in re-write, can probably use ccblade like this.
### -- SET UP CCBLADE CALL
## ccblade rotor object from blade dimensions
#ccb_rotor = ccb.Rotor(blade[i].rhub, blade[i].rtip, nbld)

## ccblade blade section object from blade data
#ccb_sections =
#    ccb.Section.(blade[i].rdim, blade[i].cdim, blade[i].tdim, rotors[i].airfoils)

## ccblade operating point object using calculated average axial velocity and induced tangential velocity.
## TODO: for some reason, the tangential velocity varies across the blade, but the axial velocity is taken as an average.  FIGURE OUT HOW TO USE CCBLADE OUT.U AND OUT.V TO MAKE THIS BETTER.
#ccb_op = ccb.OperatingPoint(
#    [average_axial_velocity for i in 1:length(induced_tangential_velocity)],
#    induced_tangential_velocity,
#    freestream.rho;
#    mu=freestream.mu,
#    asound=freestream.asound,
#)

## solve ccblade problem
#ccb_out = ccb.solve.(Ref(ccb_rotor), ccb_section, ccb_op)

##get thrust
#thrust, _ = ccb.thrusttorque(ccb_rotor, ccb_sections, ccb_out)

### - Apply CCBlade Outputs

## set B*circulation values on rotor sections
#b_circ_rotor[i, :] = 0.5 .* ccb_out.cl .* ccb_out.W .* blade[i].cdim .* nbld

"""
pretty much just use CCBlade and then do some stuff with the outputs...
"""
function initialize_blade_aero(rotors, blades, wakegrid, freestream; niter=10, rlx=0.5)

    ## -- INITIALIZE VARIABLES -- ##

    # average axial velocity
    average_axial_velocity = freestream.vinf

    # B*Gamma at grid points
    b_gamma_grid = [0.0 for i in length(wakegrid[:, 1]), j in length(wakegrid[1, :])]

    # B*gamma at rotor stations
    b_circ_rotor = [0.0 for i in length(rotors), j in length(rotors[i].radialstations)]

    # DeltaH values at rotor stations
    delta_enthalpy_grid = [0.0 for i in length(wakegrid[:, 1]), j in length(wakegrid[1, :])]

    #Induced velocities at rotor stations
    rotor_velocities = Array{RotorVelocities}(undef,length(rotors))

    for i in 1:numrotors

        #rename for convenience
        nbld = rotors[i].numblades
        yrc = rotorpanels[i].panel_centers
        yrp = blade.rdim
        omege = get_omega(rotors[i].rpm)

        # set up input velocities
        induced_axial_velocity = average_axial_velocity - freestream.vinf

        #iterate to find induced axial velocity
        for iter in 1:niter

            #initialize total thrust for current rotor
            total_thrust = 0.0

            for r in 1:(nr - 1)
                # xi = yrc[r][2] / blade[i].rtip #not needed due to custom airfoil data inputs

                if i == 1
                    #if first rotor, nothing has induced a tangential velocity yet
                    induced_tangential_velocity = 0.0
                else
                    # if not the first rotor, then previous rotors have induced tangential velocities. Find the value one grid station in front of the current rotor.
                    induced_tangential_velocity =
                        b_gamma_grid[wakegrid.rotoridxs[i] - 1, :] ./
                        (2.0 * pi .* getindex.(rotorpanels[i].panel_centers, 2))
                    #TODO: it seems that the b*Gamma values here are not aligned with the panel centers, since we have the grid centers along the grid lines, not between them at the rotor section centers... THIS WILL LIKELY LEAD TO INDEXING ERRORS SOMEWHERE
                end #if first rotor

                # Sum up velocity components to get totals in axial and tangential directions
                # TODO: for some reason the axial velocity is not sectional, but averaged. Not sure why that works when the sectional tangential comonents are used. Perhaps since this isn't applying the duct influence at this point, we can't know the sectional axial components.
                average_axial_velocity = freestream.vinf + induced_axial_velocity
                total_section_tangential_velocity =
                    induced_tangential_velocity - yrc[r] * omega

                # calculate inflow velocity
                W = sqrt(total_section_tangential_velocity^2 + average_axial_velocity^2)

                #get inflow angle
                phi = atan(average_axial_velocity, -total_section_tangential_velocity)

                #get angle of attack
                alpha = blade.twist[r] - phi

                #calculate reynolds number
                reynolds = blade.cdim[r] * abs(W) * freestream.rho / freestream.mu

                #caluclate mach number
                mach = W / freestream.asound

                #don't actually need these since using custom airfoil data inputs, also the solidity (section_sigma) is already in the rotor objects (TODO: though will need to be added to the re-interpolate rotor function )
                # section_sigma = nbld * blade.cdim[r] / (2.0 * pi * yrc[r])
                # section_stagr = 0.5*pi - blade.twist[r] #this is only used for some sort of fit from a book for applying a correction for cascades. It won't be needed if cascade data is used in the first place.

                # get airfoil data using CCBlade or similar functions.
                # NOTE: this does not include any corrections at this point. including dfdc corrections for solidity, prandtl-glauert, etc.  Assumes any airfoil data used contains any solidity/mach dependencies required.
                cl, cd = get_clcd(
                    rotors[i].airfoils[r], alpha, reynolds, mach, rotor[i].solidities[r]
                )

                # update rotor section circulation
                b_circ_new = 0.5 * cl * W * blade.cdim[r] * nblds
                b_circ_old = b_circ_rotor[i, r]
                b_circ_change = b_circ_new - b_circ_old
                b_circ_rotor[i, r] += rlx * b_circ_change #under relaxation

                # blade section length
                dr = yrp[r + 1] - yrp[r]

                #total thrust
                total_thrust -=
                    b_circ_rotor[i, r] *
                    freestream.rho *
                    total_section_tangential_velocity *
                    dr

                # wwa and wwt seem to be additional inflow values based on some sort of velocity input file for dfdc. It's probably safe to keep them as zero for now and set them to zero in functions where they are used (probably default inputs)
                # wwa, wwt = uvinfl(yrc[r])

                # set rotor slipstream velocities
                rotor_induced_axial_velocites[i, r] = induced_axial_velocity
                rotor_induced_radial_velocites[i, r] = 0.0
                rotor_induced_tangential_velocities[i, r] =
                    total_section_tangential_velocity +
                    yrc[r] * omega +
                    b_circ_rotor[i, r] / (2.0 * pi * yrc[r])
            end #for rotor panel centers

            #TODO: add for non-zero tipgap, assumes you added another section at the "end" of the rotor where the gap is. basically you'll need a panel that contributes nothing between the rotor tip and the duct wall.
            #if rotors[i].tipgap != 0.0
            #    b_circ_rotor[i,end] = 0.0
            #end

            #update induced and average axial velocities using momentum theory
            vhsq = thrust / (freestream.rho * blade[i].sweptarea) #0.5V^2 (v half square, where V is related to the induced axial velocity) from thrust equation: T = 0.5*rho*A*(vout^2 - vin^2)
            if i == 1
                #if at the first rotor, then vinf will be main component
                induced_axial_velocity =
                    -0.5 * freestream.vinf + sqrt((0.5 * freestream.vinf)^2 + vhsq)
            else
                # if not at first rotor, first rotor will have added to vinf
                average_axial_velocity = freestream.vinf + induced_axial_velocity
                induced_axial_velocity +=
                    -0.5 * average_axial_velocity +
                    sqrt((0.5 * average_axial_velocity)^2 + vhsq)
            end
        end #for iterations

        #update average_axial_velocity so that next rotor will see correct induced velocity.
        #TODO: this seems overly complicated for a julia implementation.  Probably could go back and simplify/clean things up since scoping rules are different than Fortran.
        average_axial_velocity = induced_axial_velocity + freestream.vinf

        # update grid circulation values
        update_grid_circulation!(
            b_gamma_grid,
            b_circ_rotor[i, :],
            delta_enthalpy_grid,
            omega,
            wakegrid.rotoridxs[i],
        )

        # calculate absolute and relative velocity components on rotor sections
        rotor_velocities[i] = set_rotor_velocities(
            rotor_induced_axial_velocites,
            rotor_induced_radial_velocites,
            rotor_induced_tangential_velocities,
            freestream.vinf,
            omega,
            blade.radialstations,
        )
    end #for numrotors

    #return all the initialized rotor aero arrays
    return GridAero(b_gamma_grid, b_circ_rotor, delta_enthalpy_grid), rotor_velocities
end

"""
"""
function update_grid_circulation!(
    b_gamma_grid, b_circ_rotor, delta_enthalpy_grid, omega, rotoridx
)

    #rename for convenience
    nx, nr = size(b_gamma_grid) #TODO: this might be backwards...

    #loop through radial locations
    for i in 1:(nr - 1)

        #circulation jump across rotor disk
        b_gamma_grid_change = b_circ_rotor[i]

        #enthalpy jump across rotor disk
        delta_enthalpy_grid_change = omega * b_circ_rotor[i] / (2.0 * pi)

        #convect changes downstream from current rotor index
        for j in rotoridx:(nx - 1)
            b_gamma_grid[i, j] += b_gamma_grid_change
            delta_enthalpy_grid += delta_enthalpy_grid_change
        end
    end

    return nothing
end

"""
See line 301 in rotoper.f
"""
function set_rotor_velocities(
    vax, vrad, vtan, vinf, omega, radialstations, wwa=0.0, wwt=0.0, vfac=1.0
)

    #notes:
    # 1 = axial
    # 2 = radial
    # 3 = tangential
    #

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

##################################
##### ----- FOR SOLVER ----- #####
##################################

"""
"""
function blade_section_gamma(W, chord, cl)
    return 0.5 * W * chord * cl #eqn 75 in dfdc docs
end

