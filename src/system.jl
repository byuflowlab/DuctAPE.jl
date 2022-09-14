#=
Types and Functions Related to the Ducted Rotor System

Authors: Judd Mehr,
=#

###############################
##### - EXPORTS - #####
###############################

##  TYPES

export Freestream

##  FUNCTIONS

export initialize_system_aerodynamics

#######################################
##### - COMPOSITE TYPES - #####
#######################################

"""
    Freestream{TVI,TVR,TF}

**Fields:**
 - `vinf::Float` : Freestream velocities
 - `vref::Float` : Reference velocities
 - `rho::Float` : Air density value
 - `vso::Float` : Speed of sound value
 - `mu::Float` : Air dynamic viscosity value
"""
struct Freestream{TVI,TVR,TF}
    vinf::TVI
    vref::TVR
    rho::TF
    vso::TF
    mu::TF
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
`set_grid_aero` is called at the end to initialize the circulation and enthalpy values on the wake grid (the entropy values are set to zero for now and are initialized later).
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
        0.0 for i in 1:numrotors, j in 1:length(rotors[1].radialstations)
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

                induced_axial_velocity +=
                    -0.5 * Vm_avg +
                    sqrt((0.5 * Vm_avg)^2 + vhsq)
            end
        end #for iterations

        #update Vm_avg so that next rotor will see correct induced velocity.
        #TODO: this seems overly complicated for a julia implementation.  Probably could go back and simplify/clean things up since scoping rules are different than Fortran.
        Vm_avg = induced_axial_velocity + freestream.vinf

        # initialize/update grid circulation and enthalpy values based on rotor circulation
        set_grid_aero!(
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
    return SystemAero(
        b_gamma_grid,
        delta_enthalpy_grid,
        delta_entropy_grid,
        [b_circ_rotors[i, :] for i in 1:numrotors],
        rotor_source_strengths,
        control_point_velocities,
    ),
    rotor_velocities,
    Vm_avg
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

    set_grid_aero!(grid_aerodynamics, omegas, rotoridxs)

Same, but inputs are mostly contained in grid_aerodynamics object.

**Unique Arguments:**
- `grid_aerodynamics::DuctTAPE.SystemAero` : circulations, source strengths, velocities, etc. of system grid.


    set_grid_aero!(b_gamma_grid, b_circ_rotor, delta_enthalpy_grid, omega, rotoridx)

Doesn't set entropy on grid (used in initializing rotor aerodynamics).
Also only does one rotor at a time (so inputs are for only a single rotor).


# TODOs:
- probably actually want to have a initialize_grid_aero function for initialization, then an update function for later use.
That way, things can be broken up a bit better so that not everything has to be done in the initialize_system_aerodynamics function (or maybe that can run everything, but the current content can be broken out into smaller functions for testing).
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
        r = wakegridgeometry.r_grid_points[x, nr_grid - 1]

        #similar for vm_average
        vmavg = vm_average[x, nr_grid - 1]

        #get gamma_i on duct wall TE wake
        gamma_theta[x, end] = calc_gamma_i(H, Gamma, r, vmavg)
    end

    #Taper gamma_theta along hub wall from full value at start of TE wake to start of rotor wake
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

"""
see qaic.f
call qaic1 for each control point
call qaic1 again for something if closed body? (not sure what this is doing)

in qaic1, call panaic
"""
function compute_controlpoint_velocities(
    control_points,
    panel_edges,
    wall_vorticities,
    source_strengths,
    wake_vorticities,
    vinf,
    alpha,
)

    #zero out control point velocities

    #Get velocities induced by wall (duct and hub) panels and rotor source panels (call panaic stuff)

    #Add velocities induced by wake vortex panels (line 215 in qaic.f)

    #do trailing edge panel calculations and such (see 280ish in qaic.f) (also consider moving to a different place)

    #Add freestream to velocity field

    return control_point_velocities,
    velocity_wall_jacobian, velocity_source_jacobian,
    velocity_wake_jacobian
end

"""
see panaic line 31 in aic.f
"""
function calculate_wall_and_blade_induced_velocity!(control_point_velocities, gamma_wall, gamma_hub, sigma_blades)

    # call lamp function

    # calculate velocity

    return nothing
end

#TODO: YOU ARE HERE: get the below functions fully into julia (get ins and outs figured out along with what they all do...
"""
see lampc in lamp.f
"""
function lampshade_influences(; romtol=1e-6)

    #TODO: consider replace romberg stuff with gauss-legendre quadrature
    #TODO: note that there is a singularity that needs to be delt with, so maybe wait to update integration method until later.

    # lampshade meridional length^2
    DELSQ = (X1 - X2)^2 + (R1 - R2)^2

    # reference length for convergence tolerance
    REFL = 0.5 * (R1 + R2)

    # field point is assumed to be at midpoint
    XF = 0.5 * (X1 + X2)
    RF = 0.5 * (R1 + R2)

    # evaluate integrals on increasingly fine grids
    #     (start with two intervals to avoid landing right on the midpoint)
    for IROM in 1:NROMX
        NT = 2^IROM

        UG1I[IROM] = 0.0
        VG1I[IROM] = 0.0
        UG2I[IROM] = 0.0
        VG2I[IROM] = 0.0
        US1I[IROM] = 0.0
        VS1I[IROM] = 0.0
        US2I[IROM] = 0.0
        VS2I[IROM] = 0.0

        # visit the midpoints of each of the NT intervals
        for IT in 1:NT
            T = (IT - 0.5) / NT
            TB = 1.0 - T

            DT = 1.0 / NT

            XT = X1 * TB + X2 * T
            RT = R1 * TB + R2 * T

            #TODO: write ring function find out returns
            RING!(XT, RT, XF, RF, UGT, VGT, UST, VST)

            # singular parts of velocities in the limit  XT,RT -> XF,RF
            DSQ = (XT - XF)^2 + (RT - RF)^2
            UGA =
                (RT - RF) / (4.0 * pi * DSQ) -
                0.5 * log(DSQ / (64.0 * RF^2)) / (8.0 * pi * RF)
            VGA = -(XT - XF) / (4.0 * pi * DSQ)
            USA = -(XT - XF) / (4.0 * pi * DSQ)
            VSA = -(RT - RF) / (4.0 * pi * DSQ) - 0.5 * log(DSQ / RF^2) / (8.0 * pi * RF)

            # accumulate integrals, with singular parts (at t=0.5) removed
            UG1I[IROM] = UG1I[IROM] + DT * (UGT * TB - UGA)
            VG1I[IROM] = VG1I[IROM] + DT * (VGT * TB - VGA)
            UG2I[IROM] = UG2I[IROM] + DT * (UGT * T - UGA)
            VG2I[IROM] = VG2I[IROM] + DT * (VGT * T - VGA)
            US1I[IROM] = US1I[IROM] + DT * (UST * TB - USA)
            VS1I[IROM] = VS1I[IROM] + DT * (VST * TB - VSA)
            US2I[IROM] = US2I[IROM] + DT * (UST * T - USA)
            VS2I[IROM] = VS2I[IROM] + DT * (VST * T - VSA)
        end
        # Romberg sequence using all previous grid results
        for KROM in IROM:-1:2
            # weight needed to cancel lowest-order error terms in KROM level
            W = 2.0^(2 * (IROM - KROM + 1))
            WM1 = W - 1.0

            # put Richardson extrapolation for KROM level into KROM-1 level
            UG1I[KROM - 1] = (W * UG1I[KROM] - UG1I[KROM - 1]) / WM1
            VG1I[KROM - 1] = (W * VG1I[KROM] - VG1I[KROM - 1]) / WM1
            UG2I[KROM - 1] = (W * UG2I[KROM] - UG2I[KROM - 1]) / WM1
            VG2I[KROM - 1] = (W * VG2I[KROM] - VG2I[KROM - 1]) / WM1
            US1I[KROM - 1] = (W * US1I[KROM] - US1I[KROM - 1]) / WM1
            VS1I[KROM - 1] = (W * VS1I[KROM] - VS1I[KROM - 1]) / WM1
            US2I[KROM - 1] = (W * US2I[KROM] - US2I[KROM - 1]) / WM1
            VS2I[KROM - 1] = (W * VS2I[KROM] - VS2I[KROM - 1]) / WM1
        end

        if IROM > 1
            #- compare the best-current and best-previous integrals
            ERRUG1 = UG1I[1] - UG1I[2]
            ERRVG1 = VG1I[1] - VG1I[2]
            ERRUG2 = UG2I[1] - UG2I[2]
            ERRVG2 = VG2I[1] - VG2I[2]
            ERRUS1 = US1I[1] - US1I[2]
            ERRVS1 = VS1I[1] - VS1I[2]
            ERRUS2 = US2I[1] - US2I[2]
            ERRVS2 = VS2I[1] - VS2I[2]

            ERR = max(
                abs(ERRUG1),
                abs(ERRVG1),
                abs(ERRUG2),
                abs(ERRVG2),
                abs(ERRUS1),
                abs(ERRVS1),
                abs(ERRUS2),
                abs(ERRVS2),
            )

            if ERR * REFL < ROMTOL
                convergence = false
            end
        end #if
    end #for

    DELSQ = (X1 - X2)^2 + (R1 - R2)^2
    DELS = sqrt(DELSQ)

    # analytically-integrated singular parts which were removed
    UGAI = (1.0 + log(16.0 * RF / DELS)) / (4.0 * pi * RF)
    VGAI = 0.0
    USAI = 0.0
    VSAI = (1.0 + log(2.0 * RF / DELS)) / (4.0 * pi * RF)

    # return final results, with removed parts added back on
    UG1 = (UG1I[1] + UGAI * 0.5) * DELS
    VG1 = (VG1I[1] + VGAI * 0.5) * DELS
    UG2 = (UG2I[1] + UGAI * 0.5) * DELS
    VG2 = (VG2I[1] + VGAI * 0.5) * DELS
    US1 = (US1I[1] + USAI * 0.5) * DELS
    VS1 = (VS1I[1] + VSAI * 0.5) * DELS
    US2 = (US2I[1] + USAI * 0.5) * DELS
    VS2 = (VS2I[1] + VSAI * 0.5) * DELS

    return UG, US
end

"""


     Computes the velocities induced by a vortex/source ring
     located at XV with radius RV with unit circulation
     and unit source/length density.

     Adapted from routines provided by J. Kerwin.

  Input:
     XV,RV  x,r of the ring
     XF,RF  x,r of the field point

  Output:
     UX,UR  velocity at XF,RF for unit circulation (+ counterclockwise)
     SX,SR  velocity at XF,RF for unit source/perimeter

"""
function RING!(XV, RV, XF, RF, UX, UR, SX, SR)
    if RV <= 0.0
        # zero-radius ring
        UX = 0.0
        UR = 0.0
        SX = 0.0
        SR = 0.0
        return nothing
    end

    # this fails if R=1 and X=0  (on the ring itself)
    R = RF / RV
    X = (XV - XF) / RV

    if R == 1.0 && X == 0.0
        UX = 0.0
        UR = 0.0
        SX = 0.0
        SR = 0.0
        return nothing
    end

    if RF == 0.0
        # Control Point on the axis
        UX = 1.0 / sqrt(1.0 + X^2)^3 / (2.0 * RV)
        UR = 0.0
        SX = -X / sqrt(1.0 + X^2)^3 / (2.0 * RV)
        SR = 0.0

    else
        # Control Point not on X-axis
        XRP = X^2 + (1.0 + R)^2
        XRM = X^2 + (1.0 - R)^2

        SRP = sqrt(XRP)

        AK = XRM / XRP
        ELLEK!(AK, ELE, ELK)

        F = 2.0 / XRM

        UX = (1.0 / SRP) * (ELK - ELE * (1.0 + F * (R - 1.0))) / (2.0 * PI * RV)
        UR = (X / (SRP * R)) * (ELK - ELE * (1.0 + F * R)) / (2.0 * PI * RV)

        SX = (X / SRP) * (-ELE * F) / (2.0 * PI * RV)
        SR = (1.0 / (SRP * R)) * (ELK - ELE * (1.0 + F * (R - R * R))) / (2.0 * PI * RV)

        #streamfunction due to vortex
        #       PSI = ((1.0 - 2.0*XRP*R)*ELK - ELE)*RV / (2.0*PI*sqrt(XRP))
    end

    return nothing
end

"""

ELLIPTIC FUNCTIONS ROUTINE

Adapted from routines provided by J. Kerwin.

Input
     AK     elliptic-integral argument

Output
     ELE    complete elliptic integral of the second kind
     ELK    complete elliptic integral of the first  kind

"""
function ELLEK!(AK, ELE, ELK)
    ALK = -log(AK)

    ELE = 1.00000000000
    +(0.44325141463 + (0.06260601220 + (0.04757383546 + 0.01736506451 * AK) * AK) * AK) * AK
    +(
        (0.24998368310 + (0.09200180037 + (0.04069697526 + 0.00526449639 * AK) * AK) * AK) *
        AK,
    ) * ALK

    ELK = 1.38629436112
    +(0.09666344259 + (0.03590092383 + (0.03742563713 + 0.01451196212 * AK) * AK) * AK) * AK
    +(
        0.50000000000 +
        (0.12498593597 + (0.06880248576 + (0.03328355346 + 0.00441787012 * AK) * AK) * AK) *
        AK,
    ) * ALK

    return nothing
end
