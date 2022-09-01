#=
Types and Functions Related to the Ducted Rotor System

Authors: Judd Mehr,
=#

###############################
##### ----- EXPORTS ----- #####
###############################

## -- TYPES

export OperatingConditions

## -- FUNCTIONS

export initialize_system_aerodynamics

#######################################
##### ----- COMPOSITE TYPES ----- #####
#######################################

"""
    OperatingConditions{TVI,TVR,TF}

**Fields:**
 - `vinf::Array{Float}` : Array of freestream velocities
 - `vref::Array{Float}` : Array of reference velocities
 - `rho::Float` : air density value
 - `vso::Float` : speed of sound value
 - `mu::Float` : air viscosity value
"""
struct OperatingConditions{TVI,TVR,TF}
    vinf::TVI
    vref::TVR
    rho::TF
    vso::TF
    mu::TF
    # boundarylayer::TB
end

"""
    SystemAero{TG,TH,TS,TC,TR}

**Fields:**
 - `b_gamma_grid::Matrix{Float}` : B*Γ values on wake grid
 - `delta_enthalpy_grid::Matrix{Float}` : ΔH (enthalpy) values on wake grid
 - `delta_entropy_grid::Matrix{Float}` : ΔS (entropy) values on wake grid
 - `b_circ_rotors::Array{Array{Float}}` : B*γ values at rotor blades (one array per rotor in increasing order of x location)
 - `rotor_source_strengths::Array{Array{Float}}` : σ (source strength) values at rotor blades (one array per rotor in increasing order of x location)
 - `control_point_velocities::Matrix{Float}` : velocities at control points
"""
struct SystemAero{TG,TH,TS,TC,TR,TV}
    b_gamma_grid::TG
    delta_enthalpy_grid::TH
    delta_entropy_grid::TS
    b_circ_rotors::TC
    rotor_source_strengths::TR
    control_point_velocities::TV
end

####################################
##### ----- AERODYNAMICS ----- #####
####################################

#FOR initialize_system_aerodynamics: in re-write, can probably use ccblade like this.
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
#b_circ_rotors[i, :] = 0.5 .* ccb_out.cl .* ccb_out.W .* blade[i].cdim .* nbld

"""
    initialize_system_aerodynamics(
        rotors, blades, wakegrid, freestream; niter=10, rlx=0.5
    )

Initialize system aerodynamics for rotors and wakes.

**Arguments:**
 - `rotors::Array{DuctTAPE.RotorGeometry}` : array of rotor geometries
 - `blades::Array{DuctTAPE.BladeDimensions}` : array of dimensional blade geometries
 - `wakegrid::DuctTAPE.WakeGridGeometry` : wake grid geometry
 - `freestream::DuctTAPE.OperatingConditions` : freestream information

**Returns:**
 - `systemaero::DuctTAPE.SystemAero` : aerodynamic values for rotor sections and wake grid
 - `rotorvelocities::DuctTAPE.RotorVelocities` : velocities along rotor blades
"""
function initialize_system_aerodynamics(
    rotors, blades, wakegrid, freestream; niter=10, rlx=0.5
)

    ## -- INITIALIZE VARIABLES -- ##

    # average axial velocity
    average_axial_velocity = freestream.vinf

    # B*Gamma at grid points
    b_gamma_grid = [0.0 for i in length(wakegrid[:, 1]), j in length(wakegrid[1, :])]

    # DeltaH values across wake
    delta_enthalpy_grid = [0.0 for i in length(wakegrid[:, 1]), j in length(wakegrid[1, :])]

    # DeltaS values across wake
    delta_entropy_grid = [0.0 for i in length(wakegrid[:, 1]), j in length(wakegrid[1, :])]

    # B*gamma at rotor stations
    b_circ_rotors = [0.0 for i in length(rotors), j in length(rotors[i].radialstations)]

    # rotor_source_strengths at rotor stations
    rotor_source_strengths = [
        0.0 for i in length(rotors), j in length(rotors[i].radialstations)
    ]

    #Induced velocities at rotor stations
    rotor_velocities = Array{RotorVelocities}(undef, length(rotors))

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
                b_circ_old = b_circ_rotors[i, r]
                b_circ_change = b_circ_new - b_circ_old
                b_circ_rotors[i, r] += rlx * b_circ_change #under relaxation

                # blade section length
                dr = yrp[r + 1] - yrp[r]

                #total thrust
                total_thrust -=
                    b_circ_rotors[i, r] *
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
                    b_circ_rotors[i, r] / (2.0 * pi * yrc[r])
            end #for rotor panel centers

            #TODO: add for non-zero tipgap, assumes you added another section at the "end" of the rotor where the gap is. basically you'll need a panel that contributes nothing between the rotor tip and the duct wall.
            #if rotors[i].tipgap != 0.0
            #    b_circ_rotors[i,end] = 0.0
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

        # initialize/update grid circulation and enthalpy values based on rotor circulation
        update_grid_circulation!(
            b_gamma_grid,
            b_circ_rotors[i, :],
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
    #TODO: maybe this should be panelaero rather than gridaero...
    return SystemAero(
        b_gamma_grid,
        delta_enthalpy_grid,
        delta_entropy_grid,
        [b_circ_rotors[i, :] for i in numrotors],
        rotor_source_strengths,
    ),
    rotor_velocities
end

"""
    set_grid_aero!(
        b_gamma_grid,
        delta_enthalpy_grid,
        delta_entropy_grid,
        b_circ_rotors,
        rotor_source_strengths,
        control_point_velocities,
        omegas,
        rotoridxs,
    )

sets grid aero data from rotor disk jump data

**Arguments:**
 - `b_gamma_grid::Matrix{Float}` : B*Γ values on wake grid
 - `delta_enthalpy_grid::Matrix{Float}` : ΔH (enthalpy) values on wake grid
 - `delta_entropy_grid::Matrix{Float}` : ΔS (entropy) values on wake grid
 - `b_circ_rotors::Array{Array{Float}}` : B*γ values at rotor blades (one array per rotor in increasing order of x location)
 - `rotor_source_strengths::Array{Array{Float}}` : σ (source strength) values at rotor blades (one array per rotor in increasing order of x location)
 - `control_point_velocities::Matrix{Float}` : velocities at control points
 - `omegas::Array{Float}` : rotation rate (rad/s) of the rotors
 - `rotoridxs::Array{Int}` : x-indicies of rotor locations

## Other Dispatches:

    set_grid_aero!(grid_aerodynamics, omegas, rotoridxs)

Same, but inputs are mostly contained in grid_aerodynamics object.

**Unique Arguments:**
- `grid_aerodynamics::DuctTAPE.SystemAero` : circulations, source strengths, velocities, etc. of system grid.


    set_grid_aero!(b_gamma_grid, b_circ_rotor, delta_enthalpy_grid, omega, rotoridx)

Doesn't set entropy on grid (used in initializing rotor aerodynamics).
Also only does one rotor at a time (so inputs are for only a single rotor).
"""
function set_grid_aero!(
    b_gamma_grid,
    delta_enthalpy_grid,
    delta_entropy_grid,
    b_circ_rotors,
    rotor_source_strengths,
    control_point_velocities,
    omegas,
    rotoridxs,
)

    #rename for convenience
    nx, nr = size(b_gamma_grid)

    for r in 1:length(rotors)
        for j in 1:nr

            #circulation jump across rotor disk
            b_gamma_grid_change = b_circ_rotors[r, j]

            #enthalpy jump across rotor disk
            delta_enthalpy_grid_change = omegas[r] * b_circ_rotors[r, j] / (2.0 * pi)

            #entropy jump across rotor disk (from rotor drag sources)
            #TODO: the index in question here is the index of the control point on the rotor at the station, j, we are at.
            meridional_velocity = sqrt(
                control_point_velocities[WHAT_INDEX, 1]^2 +
                control_point_velocities[WHAT_INDEX, 2]^2,
            )
            #TODO: the index in question here will be the indexs of the panel edges on either side of the rotor control point (these should be the rotor_source_panels since we're talking about source strengths here)
            sigma_avg =
                0.5 *
                (rotor_source_strengths[WHAT_INDEX] + rotor_source_strengths[WHAT_INDEX])
            delta_entropy_grid_change = meridional_velocity * sigma_avg

            #convect changes downstream from current rotor index
            for i in rotoridxs[r]:(nx - 1)
                b_gamma_grid[i, j] += b_gamma_grid_change
                delta_enthalpy_grid[i, j] += delta_enthalpy_grid_change
                delta_entropy_grid[i, j] += delta_entropy_grid_change
            end
        end
    end
    return nothing
end

function set_grid_aero!(grid_aerodynamics, omegas, rotoridxs)

    #rename for convenience
    b_gamma_grid = grid_aerodynamics.b_gamma_grid
    b_circ_rotors = grid_aerodynamics.b_circ_rotors
    delta_enthalpy_grid = grid_aerodynamics.delta_enthalpy_grid
    delta_entropy_grid = grid_aerodynamics.delta_entropy_grid
    control_point_velocities = grid_aerodynamics.control_point_velocities
    rotor_source_strengths = grid_aerodynamics.rotor_source_strengths

    set_grid_aero!(
        b_gamma_grid,
        b_circ_rotors,
        delta_enthalpy_grid,
        delta_entropy_grid,
        control_point_velocities,
        rotor_source_strengths,
        omegas,
        rotoridxs,
    )

    return nothing
end

function set_grid_aero!(b_gamma_grid, delta_enthalpy_grid, b_circ_rotor, omega, rotoridx)

    #rename for convenience
    nx, nr = size(b_gamma_grid)

    for j in 1:nr

        #circulation jump across rotor disk
        b_gamma_grid_change = b_circ_rotor[j]

        #enthalpy jump across rotor disk
        delta_enthalpy_grid_change = omega * b_circ_rotor[j] / (2.0 * pi)

        #convect changes downstream from current rotor index
        for i in rotoridx:(nx - 1)
            b_gamma_grid[i, j] += b_gamma_grid_change
            delta_enthalpy_grid[i, j] += delta_enthalpy_grid_change
        end
    end

    return nothing
end
