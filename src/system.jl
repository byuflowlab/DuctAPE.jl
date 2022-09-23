#=
Types and Functions Related to the Ducted Rotor System

Authors: Judd Mehr,
=#

####################################
##### - AERODYNAMICS - #####
####################################

#FOR initialize_system_aerodynamics: in re-write, can probably use ccblade like this.
###  SET UP CCBLADE CALL
## ccblade rotor object from blades dimensions
#ccb_rotor = ccb.Rotor(blade[i].rhub, blade[i].rtip, nbld)

## ccblade blade section object from blade data
#ccb_sections =
#    ccb.Section.(blade[i].rdim, blade[i].cdim, blade[i].tdim, rotors[i].airfoils)

## ccblade operating point object using calculated average axial velocity and induced tangential velocity.
## TODO: for some reason, the tangential velocity varies across the blade, but the axial velocity is taken as an average.  FIGURE OUT HOW TO USE CCBLADE OUT.U AND OUT.V TO MAKE THIS BETTER.
#ccb_op = ccb.OperatingPoint(
#    [Vm_avg for i in 1:length(induced_tangential_velocity)],
#    induced_tangential_velocity,
#    freestream.rho;
#    mu=freestream.mu,
#    asound=freestream.vso,
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
 - `freestream::DuctTAPE.Freestream` : freestream information

**Returns:**
 - `systemaero::DuctTAPE.SystemAero` : aerodynamic values for rotor sections and wake grid
 - `rotorvelocities::DuctTAPE.RotorVelocities` : velocities along rotor blades
 - `Vm_avg::Float` : initial guess for average axial velocity in the wakes

# Notes:
This function starts by initializing the rotor operating point from a constant freestream.
An iteration using momentum theory is performed to find the self-induced tangential velocities on the rotor.
`set_grid_flow_data` is called at the end to initialize the circulation and enthalpy values on the wake grid (the entropy values are set to zero for now and are initialized later).
`set_rotor_velocities` is also called at the end of this function to save the various rotor velocities.
"""
function initialize_system_aerodynamics(
    rotors, blades, wakegrid, rotorpanels, freestream; niter=10, rlx=0.5
)

    ##  INITIALIZE VARIABLES  ##

    # average meridional velocity
    #TODO: this would likely be a better guess for initialization (see system_initializaiton function for initializing vmavg) if it was a vector such that the effects of each rotor covered only the sections applied to them.  Right now, the final velocity (after the last rotor) is applied everywhere at initialization.
    Vm_avg = freestream.vinf

    #rename for convenience
    numrotors = length(rotors)

    # B*Gamma at grid points
    b_gamma_grid = [
        0.0 for i in 1:length(wakegrid.x_grid_points[:, 1]),
        j in 1:length(wakegrid.x_grid_points[1, :])
    ]

    # DeltaH values across wake
    delta_enthalpy_grid = [
        0.0 for i in 1:length(wakegrid.x_grid_points[:, 1]),
        j in 1:length(wakegrid.x_grid_points[1, :])
    ]

    # DeltaS values across wake
    delta_entropy_grid = [
        0.0 for i in 1:length(wakegrid.x_grid_points[:, 1]),
        j in 1:length(wakegrid.x_grid_points[1, :])
    ]

    #control point velocities
    #TODO: is this just wake grid values?
    control_point_velocities = [
        0.0 for i in 1:(length(wakegrid.x_grid_points[:, 1]) - 1),
        j in 1:(length(wakegrid.x_grid_points[1, :]) - 1)
    ]

    # B*gamma at rotor stations
    b_circ_rotors = [0.0 for i in 1:numrotors, j in 1:length(rotors[1].radialstations)]

    # rotor_source_strengths at rotor stations
    rotor_source_strengths = [
        0.0 for i in 1:numrotors, j in 1:length(rotors[1].radialstations)-1
    ]

    #Induced velocities at rotor stations
    rotor_velocities = Array{RotorVelocities}(undef, length(rotors))

    for i in 1:numrotors

        #rename for convenience
        nbld = rotors[i].numblades
        yrc = getindex.(rotorpanels[i].panel_centers, 2)
        yrp = blades[i].rdim
        omega = get_omega(rotors[i].rpm)
        nr = length(rotors[i].radialstations)

        # set up input velocities
        induced_axial_velocity = Vm_avg - freestream.vinf

        #initialize output velocities
        rotor_induced_axial_velocites = [0.0 for r in 1:nr]
        rotor_induced_radial_velocites = [0.0 for r in 1:nr]
        rotor_induced_tangential_velocities = [0.0 for r in 1:nr]

        #iterate to find induced axial velocity
        for iter in 1:niter

            #initialize total thrust for current rotor
            total_thrust = 0.0

            for r in 1:(nr - 1)
                # xi = yrc[r] / blade[i].rtip #not needed due to custom airfoil data inputs

                if i == 1
                    #if first rotor, nothing has induced a tangential velocity yet
                    induced_tangential_velocity = 0.0
                else
                    # if not the first rotor, then previous rotors have induced tangential velocities. Find the value one grid station in front of the current rotor.
                    # comes from: v_θ = Γ/2πr, eqn 48 (or really 70 and 71) in dfdc docs
                    induced_tangential_velocity =
                        b_gamma_grid[wakegrid.rotoridxs[i] - 1, r] / (2.0 * pi .* yrc[r])
                end #if first rotor

                # Sum up velocity components to get totals in axial and tangential directions
                # TODO: for some reason the axial velocity is not sectional, but averaged. Not sure why that works when the sectional tangential comonents are used. Perhaps since this isn't applying the duct influence at this point, we can't know the sectional axial components.
                Vm_avg = freestream.vinf + induced_axial_velocity
                total_section_tangential_velocity =
                    induced_tangential_velocity - yrc[r] * omega

                # calculate inflow velocity
                W = sqrt(total_section_tangential_velocity^2 + Vm_avg^2)

                #get inflow angle
                phi = atan(Vm_avg, -total_section_tangential_velocity)

                #get angle of attack
                alpha = blades[i].tdim[r] - phi

                #calculate reynolds number
                reynolds = blades[i].cdim[r] * abs(W) * freestream.rho / freestream.mu

                #caluclate mach number
                mach = W / freestream.vso

                #don't actually need these since using custom airfoil data inputs, also the solidity (section_sigma) is already in the rotor objects (TODO: though will need to be added to the re-interpolate rotor function )
                # section_sigma = nbld * blades.cdim[r] / (2.0 * pi * yrc[r])
                # section_stagr = 0.5*pi - blades.twist[r] #this is only used for some sort of fit from a book for applying a correction for cascades. It won't be needed if cascade data is used in the first place.

                # get airfoil data using CCBlade or similar functions.
                # NOTE: this does not include any corrections at this point. including dfdc corrections for solidity, prandtl-glauert, etc.  Assumes any airfoil data used contains any solidity/mach dependencies required.
                cl, cd = get_clcd(
                    rotors[i].airfoils[r], alpha, reynolds, mach, rotors[i].solidities[r]
                )

                # update rotor section circulation
                b_circ_new = 0.5 * cl * W * blades[i].cdim[r] * nbld #eqn 1.82?
                b_circ_old = b_circ_rotors[i, r]
                b_circ_change = b_circ_new - b_circ_old
                b_circ_rotors[i, r] += rlx * b_circ_change #under relaxation

                # blades section length
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
                rotor_induced_axial_velocites[r] = induced_axial_velocity
                rotor_induced_radial_velocites[r] = 0.0
                rotor_induced_tangential_velocities[r] =
                    total_section_tangential_velocity +
                    yrc[r] * omega +
                    b_circ_rotors[i, r] / (2.0 * pi * yrc[r])
            end #for rotor panel centers

            #TODO: add for non-zero tipgap, assumes you added another section at the "end" of the rotor where the gap is. basically you'll need a panel that contributes nothing between the rotor tip and the duct wall.
            #if rotors[i].tipgap != 0.0
            #    b_circ_rotors[i,end] = 0.0
            #end

            #update induced and average axial velocities using momentum theory
            vhsq = total_thrust / (freestream.rho * blades[i].sweptarea) #0.5V^2 (v half square, where V is related to the induced axial velocity) from thrust equation: T = 0.5*rho*A*(vout^2 - vin^2)
            if i == 1
                #if at the first rotor, then vinf will be main component
                #NOTE: this is not in the dfdc documentation. Looks to simply be a first guess to set up the system before the solver can be called.
                induced_axial_velocity =
                    -0.5 * freestream.vinf + sqrt((0.5 * freestream.vinf)^2 + vhsq)
            else
                # if not at first rotor, first rotor will have added to vinf
                Vm_avg = freestream.vinf + induced_axial_velocity

                induced_axial_velocity += -0.5 * Vm_avg + sqrt((0.5 * Vm_avg)^2 + vhsq)
            end
        end #for iterations

        #update Vm_avg so that next rotor will see correct induced velocity.
        #TODO: this seems overly complicated for a julia implementation.  Probably could go back and simplify/clean things up since scoping rules are different than Fortran.
        Vm_avg = induced_axial_velocity + freestream.vinf

        # initialize/update grid circulation and enthalpy values based on rotor circulation
        set_grid_flow_data!(
            b_gamma_grid,
            delta_enthalpy_grid,
            b_circ_rotors[i, :],
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
            blades[i].rdim,
        )
    end #for numrotors

    #return all the initialized rotor aero arrays
    #TODO: maybe this should be panelaero rather than gridaero...
    return b_gamma_grid,
        delta_enthalpy_grid,
        delta_entropy_grid,
        [b_circ_rotors[i, :] for i in 1:numrotors],
        rotor_source_strengths,
        control_point_velocities,
    rotor_velocities,
    Vm_avg
end

"""
    set_grid_flow_data!(
        b_gamma_grid,
        delta_enthalpy_grid,
        delta_entropy_grid,
        b_circ_rotors,
        rotor_source_strengths,
        control_point_velocities,
        omegas,
        rotoridxs,
    )

Sets grid aero data from rotor disk jump data

Aero data includes the three "tilde" values in the DFDC theory: circulation (Gamma), enthalpy (H), and entropy (S), which are all defined based on rotor values, then carried downstream at each wake station, summing up effects from each rotor.

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

    set_grid_flow_data!(grid_aerodynamics, omegas, rotoridxs)

Same, but inputs are mostly contained in grid_aerodynamics object.

**Unique Arguments:**
- `grid_aerodynamics::DuctTAPE.SystemAero` : circulations, source strengths, velocities, etc. of system grid.


    set_grid_flow_data!(b_gamma_grid, b_circ_rotor, delta_enthalpy_grid, omega, rotoridx)

Doesn't set entropy on grid (used in initializing rotor aerodynamics).
Also only does one rotor at a time (so inputs are for only a single rotor).


# TODOs:
- probably actually want to have a initialize_grid_aero function for initialization, then an update function for later use.
That way, things can be broken up a bit better so that not everything has to be done in the initialize_system_aerodynamics function (or maybe that can run everything, but the current content can be broken out into smaller functions for testing).
"""
function set_grid_flow_data!(
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
            #from eqn 1.28
            delta_enthalpy_grid_change = omegas[r] * b_circ_rotors[r, j] / (2.0 * pi)

            #entropy jump across rotor disk (from rotor drag sources)
            meridional_velocity = sqrt(
                control_point_velocities[rotoridxs[r], j][1]^2 +
                control_point_velocities[rotoridxs[r], j][2]^2,
            )
            #TODO: the index in question here will be the indexs of the panel edges on either side of the rotor control point (these should be the rotor_source_panels since we're talking about source strengths here)
            sigma_avg =
                0.5 *
                (rotor_source_strengths[WHAT_INDEX] + rotor_source_strengths[WHAT_INDEX])

            #TODO: this doesn't seem to match either entropy expression in dfdc docs.
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

function set_grid_flow_data!(grid_aerodynamics, omegas, rotoridxs)

    #rename for convenience
    b_gamma_grid = grid_aerodynamics.b_gamma_grid
    b_circ_rotors = grid_aerodynamics.b_circ_rotors
    delta_enthalpy_grid = grid_aerodynamics.delta_enthalpy_grid
    delta_entropy_grid = grid_aerodynamics.delta_entropy_grid
    control_point_velocities = grid_aerodynamics.control_point_velocities
    rotor_source_strengths = grid_aerodynamics.rotor_source_strengths

    set_grid_flow_data!(
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

function set_grid_flow_data!(
    b_gamma_grid, delta_enthalpy_grid, b_circ_rotor, omega, rotoridx
)
    #This function is analogous to ROTBG2GRD in dfdc inigrid.f line 434

    #rename for convenience
    nx, nr = size(b_gamma_grid)

    for j in 1:nr

        #circulation jump across rotor disk
        b_gamma_grid_change = b_circ_rotor[j]

        #enthalpy jump across rotor disk
        #from eqn 1.28
        delta_enthalpy_grid_change = omega * b_circ_rotor[j] / (2.0 * pi)

        #convect changes downstream from current rotor index
        for i in rotoridx:(nx - 1)
            b_gamma_grid[i, j] += b_gamma_grid_change
            delta_enthalpy_grid[i, j] += delta_enthalpy_grid_change
        end
    end

    return nothing
end

"""
TODO: this would be much better if in the system aero initialization function, the average axial velocity variable was saved before/after each rotor. therefore before rotor 1, everythign would be vinf, then after rotor 1, you'd have vinf + the induced velcity, then for subsequent rotors you add the additional induced velocities.
There's also not really any need to do things in a loop here like in dfdc, everything can just be defined directly using the rotoridxs indices, etc.
"""
function initialize_vm_average(initial_guesses, wakegrid, vinf)

    #get indices for convenience
    nx, nr = size(wakegrid.x_grid_points)

    vm_average = [0.0 for i in 1:nx, j in 1:nr]

    for i in 1:nx
        for j in 1:nr
            if j == 1 || j == nr
                vm_average[i, j] = 0.5 * (initial_guesses + vinf)
            else
                vm_average[i, j] = initial_guesses
            end
        end
    end

    return vm_average
end

"""
gamma_theta values at panel centers
uses eqn 45 in dfdc theory
"""
function calculate_gamma_theta(system_aero, wakegrid, vm_average)

    #rename for convenience
    nx_grid, nr_grid = size(wakegrid.x_grid_points)
    b_gamma_grid = system_aero.b_gamma_grid
    delta_enthalpy_grid = system_aero.delta_enthalpy_grid

    # initialize gamma_theta
    gamma_theta = [0.0 for i in 1:nx_grid, j in 1:nr_grid]

    #set gamma_theta on wakes between walls
    for r in 2:(numradialstations - 1) #loop through internal rotor radial stations (between walls)
        for x in 1:nx_grid #loop through wake panels (assume wake starts at first rotor

            #H is ΔH2 - ΔH1, (where ΔH is constant on rotor panels, i.e., between vortex wakes)
            H = (delta_enthalpy_grid[x, r], delta_enthalpy_grid[x, r - 1])

            #similar for Gamma
            Gamma = (b_gamma_grid[x, r], b_gamma_grid[x, r - 1])

            #use the r at the location we care about (radial station where wake is generated from)
            r = wakegridgeometry.r_grid_points[x, r]

            #similar for vm_average
            vmavg = vm_average[x, r]

            #get gamma_i on vortex wake panel edge
            gamma_theta[x, r] = calc_gamma_i(H, Gamma, r, vmavg)
        end
    end

    #set gamma_theta on hub trailing edge wake
    for x in (wakegrid.hubTEidx):nx_grid

        #no enthalpy jump across wall
        #TODO: why is this not (0.0, delta_enthalpy_grid[x,1])?
        H = (0.0, 0.0)

        #no Gamma inside wall
        Gamma = (0.0, b_gamma_grid[x, 1])

        #r is wall r
        r = wakegridgeometry.r_grid_points[x, 1]

        #similar for vm_average
        vmavg = vm_average[x, 1]

        #get gamma_i on hub TE wake
        gamma_theta[x, 1] = calc_gamma_i(H, Gamma, r, vmavg)
    end

    #Taper gamma_theta along hub wall from full value at start of TE wake to start of rotor wake
    #TODO: is this really necessary?  What does it do?
    for x in 1:(wakegrid.hubTEidx - 1)
        taper = 1.0 - (wakegrid.hubTEidx - (x - 1)) / wakegrid.hubTEidx
        gamma_theta[x, 1] = gamma_theta[hubTEidx, 1] * taper
    end

    #set gamma_theta on duct trailing edge wake
    for x in (wakegrid.hubTEidx):nx_grid

        #no enthalpy jump across wall
        H = (delta_enthalpy_grid[x, nr_grid - 1], 0.0)

        #no Gamma inside wall
        Gamma = (b_gamma_grid[x, nr_grid - 1], 0.0)

        #r is wall r
        #TODO: why is it nr-1 then? should it not be nr?
        r = wakegridgeometry.r_grid_points[x, nr_grid - 1]

        #similar for vm_average
        vmavg = vm_average[x, nr_grid - 1]

        #get gamma_i on duct wall TE wake
        gamma_theta[x, end] = calc_gamma_i(H, Gamma, r, vmavg)
    end

    #Taper gamma_theta along hub wall from full value at start of TE wake to start of rotor wake
    #TODO: again, why is this done?
    for x in 1:(wakegrid.wallTEidx - 1)
        taper = 1.0 - (wakegrid.wallTEidx - (x - 1)) / wakegrid.wallTEidx
        gamma_theta[x, end] = gamma_theta[wallTEidx, end] * taper
    end

    #return gamma_theta
    return gamma_theta
end

"""
equation 45 (or 61) in dfdc theory
"""
function calc_gamma_i(H, Gamma, r, vmavg)
    if vmag == 0.0
        #don't divide by zero...
        return 0.0
    else
        return ((H[2] - H[1]) - 0.5 * (Gamma[2]^2 - Gamma[1]^2) / (2.0 * pi * r)^2) / vmag
    end
end
