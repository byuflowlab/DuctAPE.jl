#=

Various Post-processing functions

=#

"""

Post-processes the converged rotor and wake strengths and outputs various performance metrics.

# Arguments:
- `Gamr::Matrix{Float}` :
- `sigr::Matrix{Float}` :
- `gamw::vector{Float}` :
- `inputs::NamedTuple` :

# Keyword Arguments:
- `write_outputs::Bool=false` :
- `outfile::String="dfdc_data.jl"` :
- `checkoutfileexists::Bool=false` :
- `tuple_name::String="ductape_data"` :

# Returns:
- `outputs::NamedTuple` : Named tuple of output values including:
  - `Gamr::Matrix{Float}` : Circulation strengths.
  - `gamw::Vector{Float}` : Wake panel node strengths.
  - `sigr::Matrix{Float}` : Rotor source panel node strengths.
  - `gamb::Vector{Float}` : Body panel node strengths.
  - `intermediate_values::NamedTuple` : Values used in intermediate calculations.
  - `reference_values::NamedTuple` : Reference Values.
  - `bodies::NamedTuple` : Outputs associated with the duct and center body.
  - `body_wakes::NamedTuple` : Outputs associated with the duct and centerbody wakes.
  - `rotors::NamedTuple` : Outputs associated with the rotors.
  - `totals::NamedTuple` : Outputs associated with the system as a whole.

See extended help for details of what the inputs and outputs objects contain.

# Extended help

### Input Details:

  - **`inputs`**

### Output Details:

  - **`intermediate_values`**
      - `vz_rotor::Matrix{Float}` :
      - `vr_rotor::Matrix{Float}` :
      - `vtheta_rotor::Matrix{Float}` :
      - `Wz_rotor::Matrix{Float}` :
      - `Wtheta_rotor::Matrix{Float}` :
      - `Wm_rotor::Matrix{Float}` :
      - `Wmag_rotor::Matrix{Float}` :
      - `vz_rotor_b::Matrix{Float}` :
      - `vr_rotor_b::Matrix{Float}` :
      - `vz_rotor_w::Matrix{Float}` :
      - `vr_rotor_w::Matrix{Float}` :
      - `vz_rotor_r::Matrix{Float}` :
      - `vr_rotor_r::Matrix{Float}` :
      - `vz_wake::Matrix{Float}` :
      - `vr_wake::Matrix{Float}` :
      - `Wz_wake::Matrix{Float}` :
      - `Wm_wake::Matrix{Float}` :
      - `vz_wake_b::Matrix{Float}` :
      - `vr_wake_b::Matrix{Float}` :
      - `vz_wake_w::Matrix{Float}` :
      - `vr_wake_w::Matrix{Float}` :
      - `vz_wake_r::Matrix{Float}` :
      - `vr_wake_r::Matrix{Float}` :
      - `Gamma_tilde::Matrix{Float}` :
      - `H_tilde::Matrix{Float}` :
      - `phi::Matrix{Float}` :
      - `alpha::Matrix{Float}` :
      - `cl::Matrix{Float}` :
      - `cd::Matrix{Float}` :
  - **`reference_values`**
    - `Vinf::Float` :
    - `Vref::Float` :
    - `Rref::Float` :
    - `rhoinf::Float` :
    - `muinf::Float` :
    - `asound::Float` :
    - `Omega::Float` :
    - `B::Float` :
  - **`bodies`**
    - `thrust::Float` :
    - `cp_casing::Vector{Float}` :
    - `zpts_casing::Vector{Float}` :
    - `cp_nacelle::Vector{Float}` :
    - `zpts_nacelle::Vector{Float}` :
    - `cp_centerbody::Vector{Float}` :
    - `zpts_centerbody::Vector{Float}` :
    - `Vtan::Vector{Float}` :
    - `vtan_casing::Vector{Float}` :
    - `vtan_nacelle::Vector{Float}` :
    - `vtan_centerbody::Vector{Float}` :
    - `vtan_vinf::Vector{Float}` :
    - `vtan_body::Vector{Float}` :
    - `vtan_casing_b::Vector{Float}` :
    - `vtan_nacelle_b::Vector{Float}` :
    - `vtan_centerbody_b::Vector{Float}` :
    - `vtan_wake::Vector{Float}` :
    - `vtan_casing_w::Vector{Float}` :
    - `vtan_nacelle_w::Vector{Float}` :
    - `vtan_centerbody_w::Vector{Float}` :
    - `vtan_rotors::Vector{Float}` :
    - `vtan_casing_r::Vector{Float}` :
    - `vtan_nacelle_r::Vector{Float}` :
    - `vtan_centerbody_r::Vector{Float}` :
  - **`body_wakes`**
    - `vtan_ductwake::Vector{Float}` :
    - `cp_ductwake::Vector{Float}` :
    - `zpts_ductwake::Vector{Float}` :
    - `vtan_hubwake::Vector{Float}` :
    - `cp_hubwake::Vector{Float}` :
    - `zpts_hubwake::Vector{Float}` :
  - **`rotors`**
    - `efficiency::Vector{Float}` :
    - `inviscid_thrust::Vector{Float}` :
    - `inviscid_thrust_dist::Matrix{Float}` :
    - `viscous_thrust::Vector{Float}` :
    - `viscous_thrust_dist::Matrix{Float}` :
    - `thrust::Vector{Float}` :
    - `CT::Vector{Float}` :
    - `inviscid_torque::Vector{Float}` :
    - `inviscid_torque_dist::Matrix{Float}` :
    - `viscous_torque::Vector{Float}` :
    - `viscous_torque_dist::Matrix{Float}` :
    - `torque::Vector{Float}` :
    - `CQ::Vector{Float}` :
    - `inviscid_power::Vector{Float}` :
    - `inviscid_power_dist::Matrix{Float}` :
    - `viscous_power::Vector{Float}` :
    - `viscous_power_dist::Matrix{Float}` :
    - `power::Vector{Float}` :
    - `CP::Vector{Float}` :
    - `cl::Matrix{Float}` :
    - `cd::Matrix{Float}` :
    - `alpha::Matrix{Float}` :
    - `phi::Matrix{Float}` :
    - `blade_normal_force_per_unit_span::Matrix{Float}` :
    - `blade_tangential_force_per_unit_span::Matrix{Float}` :
  - **`totals`**
    - `thrust::Float` :
    - `torque::Float` :
    - `power::Float` :
    - `CT::Float` :
    - `CQ::Float` :
    - `CP::Float` :
    - `total_efficiency::Float` :
    - `induced_efficiency::Float` :
    - `ideal_efficiency::Float` :
"""
function post_process(
    Gamr,
    sigr,
    gamw,
    inputs;
    write_outputs=false,
    outfile="dfdc_data.jl",
    checkoutfileexists=false,
    tuple_name="ductape_data",
)

    # - things contained in iv tuple
    # vz_rotor,
    # vr_rotor,
    # vtheta_rotor,
    # Wz_rotor,
    # Wtheta_rotor,
    # Wm_rotor,
    # Wmag_rotor,
    # vzfrombody,
    # vr_rotor_b,
    # vz_rotor_w,
    # vr_rotor_w,
    # vz_rotor_r,
    # vr_rotor_r,
    # Gamma_tilde,
    # H_tilde,
    # phi,
    # alpha,
    # cl,
    # cd,

    iv = get_intermediate_values(Gamr, sigr, gamw, inputs) #intermediate values

    # get problem dimensions
    nr, nrotor = size(Gamr)
    nw = nr + 1

    # - Extract convenient input fields - #
    Vinf = inputs.freestream.Vinf
    Vref = inputs.reference_parameters.Vref
    Rref = inputs.reference_parameters.Rref
    rhoinf = inputs.freestream.rhoinf
    muinf = inputs.freestream.muinf
    asound = inputs.freestream.asound
    didr = inputs.ductidsaftofrotors
    hidr = inputs.hubidsaftofrotors
    Omega = inputs.blade_elements.Omega
    B = inputs.blade_elements.B
    chord = reshape(reduce(vcat, inputs.blade_elements.chords), (nr, nrotor))
    twist = reshape(reduce(vcat, inputs.blade_elements.twists), (nr, nrotor))
    #stagger is twist angle but from axis
    stagger = 0.5 * pi .- twist
    solidity = reshape(reduce(vcat, inputs.blade_elements.solidity), (nr, nrotor))
    afparamsin = reshape(reduce(vcat, inputs.blade_elements.inner_airfoil), (nr, nrotor))
    afparamsout = reshape(reduce(vcat, inputs.blade_elements.outer_airfoil), (nr, nrotor))
    innerfrac = reshape(reduce(vcat, inputs.blade_elements.inner_fraction), (nr, nrotor))
    rpc = inputs.rotor_panel_centers
    rpe = inputs.rotor_panel_edges
    Rhub = rpe[1, :]
    Rtip = rpe[end, :]
    rpl = reshape(
        reduce(vcat, (p -> p.influence_length).(inputs.rotor_source_panels)), (nr, nrotor)
    )

    ref = (; Vinf, Vref, Rref, rhoinf, muinf, asound, Omega, B)

    # - Rotor Thrust - #
    # inviscid thrust
    rotor_inviscid_thrust, rotor_inviscid_thrust_dist = inviscid_rotor_trust(
        iv.Wtheta_rotor, iv.Gamma_tilde, rpl, rhoinf
    )

    # viscous thrust
    rotor_viscous_thrust, rotor_viscous_thrust_dist = viscous_rotor_thrust(
        iv.Wz_rotor, iv.Wmag_rotor, B, chord, rpl, iv.cd, rhoinf
    )

    # total thrust
    rotor_thrust = rotor_inviscid_thrust .+ rotor_viscous_thrust

    # - Rotor Torque - #
    # inviscid torque
    rotor_inviscid_torque, rotor_inviscid_torque_dist = inviscid_rotor_torque(
        iv.Wz_rotor, rpc, rpl, iv.Gamma_tilde, rhoinf
    )

    # viscous torque
    rotor_viscous_torque, rotor_viscous_torque_dist = viscous_rotor_torque(
        iv.Wtheta_rotor, iv.Wmag_rotor, B, chord, rpc, rpl, iv.cd, rhoinf
    )

    # - Rotor Power - #
    # inviscid power
    rotor_inviscid_power = inviscid_rotor_power(rotor_inviscid_torque', Omega)

    rotor_inviscid_power_dist = similar(rotor_inviscid_torque_dist) .= 0.0
    for ir in 1:nrotor
        rotor_inviscid_power_dist[:, ir] = inviscid_rotor_power(
            rotor_inviscid_torque_dist[:, ir], Omega[ir]
        )
    end

    # viscous power
    rotor_viscous_power = viscous_rotor_power(rotor_viscous_torque', Omega)
    rotor_viscous_power_dist = similar(rotor_viscous_torque_dist) .= 0.0
    for ir in 1:nrotor
        rotor_inviscid_power_dist[:, ir] = viscous_rotor_power(
            rotor_inviscid_torque_dist[:, ir], Omega[ir]
        )
    end

    ## -- Surface Velocity on Bodies -- ##
    vtan_tuple = get_body_tangential_velocities(inputs, gamw, sigr)

    (;
        # Total and components
        Vtan,
        vtan_vinf, # tangential velocity from freestream
        vtan_body,
        body_jump_term,
        vtan_wake,
        vtan_rotors,
        # Splits:
        # duct inner surface (casing)
        vtan_casing,
        vtan_casing_b,
        casing_jump_term,
        vtan_casing_w,
        vtan_casing_r,
        # duct outer surface (nacelle)
        vtan_nacelle,
        vtan_nacelle_b,
        nacelle_jump_term,
        vtan_nacelle_w,
        vtan_nacelle_r,
        # center body
        vtan_centerbody,
        vtan_centerbody_b,
        centerbody_jump_term,
        vtan_centerbody_w,
        vtan_centerbody_r,
    ) = vtan_tuple

    ## -- Pressure on Bodies -- ##
    _, cp_casing, cp_nacelle, cp_centerbody, _, zpts_casing, zpts_nacelle, zpts_centerbody = get_body_cps(
        Vtan,
        Gamr,
        sigr,
        iv.Wm_rotor,
        Vinf,
        Vref,
        B,
        Omega,
        didr,
        hidr,
        inputs.body_vortex_panels,
        inputs.isduct,
    )

    ## -- Pressure on Body Wakes -- ##
    cp_ductwake, vtan_ductwake = get_bodywake_cps(
        Gamr,
        inputs.vz_dww,
        inputs.vr_dww,
        gamw,
        inputs.vz_dwr,
        inputs.vr_dwr,
        sigr,
        inputs.vz_dwb,
        inputs.vr_dwb,
        inputs.gamb[1:(inputs.body_vortex_panels.totnode)],
        inputs.duct_wake_panels,
        iv.Wm_rotor,
        Omega,
        B,
        Vinf,
        Vref;
        body="duct",
    )

    cp_hubwake, vtan_hubwake = get_bodywake_cps(
        Gamr,
        inputs.vz_hww,
        inputs.vr_hww,
        gamw,
        inputs.vz_hwr,
        inputs.vr_hwr,
        sigr,
        inputs.vz_hwb,
        inputs.vr_hwb,
        inputs.gamb[1:(inputs.body_vortex_panels.totnode)],
        inputs.hub_wake_panels,
        iv.Wm_rotor,
        Omega,
        B,
        Vinf,
        Vref;
        body="hub",
    )

    ## -- Duct Outputs -- ##
    # - Put duct pressures together - #
    duct_cp = [cp_casing; cp_nacelle]

    # - Calculate Thrust from Bodies - #

    body_thrust, _ = forces_from_pressure(
        duct_cp, inputs.body_vortex_panels; rhoinf=rhoinf, Vref=Vref
    )

    ## -- Total Outputs -- ##

    # - Total Thrust - #
    rotor_thrust = rotor_inviscid_thrust .+ rotor_viscous_thrust
    total_thrust = sum([rotor_inviscid_thrust'; rotor_viscous_thrust'; body_thrust])

    # - Total Torque - #
    rotor_torque = rotor_inviscid_torque .+ rotor_viscous_torque
    total_torque = sum([rotor_inviscid_torque; rotor_viscous_torque])

    # - Total Power - #
    rotor_power = rotor_inviscid_power .+ rotor_viscous_power
    total_power = sum([rotor_inviscid_power; rotor_viscous_power])

    # - Total Efficiency - #
    rotor_efficiency = get_total_efficiency(rotor_thrust, rotor_power, Vinf)
    total_efficiency = get_total_efficiency(total_thrust, total_power, Vinf)

    # - Induced Efficiency - #
    induced_efficiency = [
        get_induced_efficiency(
            rotor_inviscid_thrust[ir], body_thrust, rotor_inviscid_power[ir], Vinf
        ) for ir in 1:nrotor
    ]

    # - Ideal Efficiency - #
    ideal_efficiency = get_ideal_efficiency(total_thrust, rhoinf, Vinf, Rref)

    # - Blade Loading - #
    blade_normal_force_per_unit_span, blade_tangential_force_per_unit_span = get_blade_loads(
        iv.Wmag_rotor, iv.phi, iv.cl, iv.cd, chord, rhoinf
    )

    # - Thrust and Torque Coefficients - #
    rotor_CT, rotor_CQ, rotor_CP = tqpcoeff(rotor_thrust, rotor_torque, rhoinf, Omega, Rref)
    total_CT, total_CQ, total_CP = tqpcoeff(total_thrust, total_torque, rhoinf, Omega, Rref)

    ## -- Assemble Output Tuple -- ##

    #TODO: standardize the output naming conventions
    outs = (;
        # - States - #
        Gamr,
        gamw,
        sigr,
        gamb=inputs.gamb,
        # - Intermediate Values - #
        intermediate_values=iv,
        # - Reference Values - #
        reference_values=ref,
        # - Body Values - #
        bodies=(;
            # body thrust
            thrust=body_thrust,
            # surface pressures
            cp_casing,
            zpts_casing,
            cp_nacelle,
            zpts_nacelle,
            cp_centerbody,
            zpts_centerbody,
            #individual body velocity contributions
            Vtan,
            vtan_casing,
            vtan_nacelle,
            vtan_centerbody,
            vtan_vinf,
            vtan_body,
            vtan_casing_b,
            vtan_nacelle_b,
            vtan_centerbody_b,
            vtan_wake,
            vtan_casing_w,
            vtan_nacelle_w,
            vtan_centerbody_w,
            vtan_rotors,
            vtan_casing_r,
            vtan_nacelle_r,
            vtan_centerbody_r,
        ),
        # - Body Wake Values - #
        body_wakes=(;
            # surface velocities and pressures
            vtan_ductwake,
            cp_ductwake,
            zpts_ductwake=inputs.duct_wake_panels.controlpoint[1, :],
            vtan_hubwake,
            cp_hubwake,
            zpts_hubwake=inputs.hub_wake_panels.controlpoint[1, :],
        ),
        # - Rotor Values - #
        rotors=(;
            efficiency=rotor_efficiency,
            inviscid_thrust=rotor_inviscid_thrust,
            inviscid_thrust_dist=rotor_inviscid_thrust_dist,
            viscous_thrust=rotor_viscous_thrust,
            viscous_thrust_dist=rotor_viscous_thrust_dist,
            thrust=rotor_thrust,
            CT = rotor_CT,
            # rotor torque
            inviscid_torque=rotor_inviscid_torque,
            inviscid_torque_dist=rotor_inviscid_torque_dist,
            viscous_torque=rotor_viscous_torque,
            viscous_torque_dist=rotor_viscous_torque_dist,
            torque=rotor_torque,
            CQ = rotor_CQ,
            # rotor power
            inviscid_power=rotor_inviscid_power,
            inviscid_power_dist=rotor_inviscid_power_dist,
            viscous_power=rotor_viscous_power,
            viscous_power_dist=rotor_viscous_power_dist,
            power=rotor_power,
            CP = rotor_CP,
            # - Blade Element Values - #
            cl=iv.cl,
            cd=iv.cd,
            alpha=iv.alpha,
            phi=iv.phi,
            blade_normal_force_per_unit_span,
            blade_tangential_force_per_unit_span,
        ),
        # - Total Values - #
        totals=(;
            thrust=total_thrust,
            torque=total_torque,
            power=total_power,
            CT=total_CT[1],
            CQ=total_CQ[1],
            CP=total_CP[1],
            total_efficiency=total_efficiency[1],
            induced_efficiency,
            ideal_efficiency,
        ),
    )

    if write_outputs
        write_data(
            outs, outfile; tuple_name=tuple_name, checkoutfileexists=checkoutfileexists
        )
    end

    return outs
end

######################################################################
#                                                                    #
#                        Velocity Functions                          #
#                                                                    #
######################################################################

"""
"""
function get_body_tangential_velocities(inputs, gamw, sigr)
    # - rename for convenience - #
    (; AICt, AICtr, AICtw, body_vortex_panels, gamb) = inputs
    nrotor = size(sigr, 2)

    # - Influence from Freestream - #
    Vs = inputs.freestream.Vinf * [1.0; 0.0] # axisymmetric, so no radial component
    Vsmat = repeat(Vs; outer=(1, size(inputs.body_vortex_panels.controlpoint, 2))) # need velocity on each panel
    vtan_vinf = [dot(v, t) for (v, t) in zip(eachcol(Vsmat), eachcol(body_vortex_panels.tangent))]

    # initialize total
    Vtan = copy(vtan_vinf)

    # - Velocity Contributions from body - #

    # Panel-on-panel influence
    vtan_body = AICt * gamb[1:size(AICt, 2)]
    Vtan .+= vtan_body

    # Jump Term
    jumpduct = (gamb[1:(body_vortex_panels.npanel[1])] + gamb[2:(body_vortex_panels.nnode[1])]) / 2
    jumphub =
        (
            gamb[(body_vortex_panels.nnode[1]):(body_vortex_panels.totpanel)] +
            gamb[(body_vortex_panels.nnode[1] + 1):(body_vortex_panels.totpanel + 1)]
        ) / 2.0
    body_jump_term = [jumpduct; jumphub]
    Vtan .-= body_jump_term / 2.0

    # - Velocity Contributions from wake - #
    vtan_wake = AICtw * gamw
    Vtan .+= vtan_wake

    # - Velocity Contributions from rotors - #
    vtan_rotors = similar(vtan_wake,(length(Vtan),nrotor))
    for jrotor in 1:nrotor
        @views vtan_rotors[:,jrotor] = AICtr[:,:,jrotor] * sigr[:,jrotor]
        Vtan .+= vtan_rotors[:,jrotor]
    end

    # - Split Velocities associates with inner and outer duct and hub - #
    # total tangential velocities
    vtan_casing, vtan_nacelle, vtan_centerbody, _, _, _ = split_bodies(Vtan, body_vortex_panels)

    # tangential velocities due to body
    vtan_casing_b, vtan_nacelle_b, vtan_centerbody_b, _, _, _ = split_bodies(
        vtan_body, body_vortex_panels
    )

    # jump terms
    casing_jump_term, nacelle_jump_term, centerbody_jump_term, _, _, _ = split_bodies(body_jump_term, body_vortex_panels)

    # tangential velocities due to wake
    vtan_casing_w, vtan_nacelle_w, vtan_centerbody_w, _, _, _ = split_bodies(vtan_wake, body_vortex_panels)

    # tangential velocities due to rotors
    vtan_casing_r, vtan_nacelle_r, vtan_centerbody_r, _, _, _ = split_bodies(vtan_rotors, body_vortex_panels)

    return (;
        # Total and components
        Vtan,
        vtan_vinf, # tangential velocity from freestream
        vtan_body,
        body_jump_term,
        vtan_wake,
        vtan_rotors,
        # Splits:
        # duct inner surface (casing)
        vtan_casing,
        vtan_casing_b,
        casing_jump_term,
        vtan_casing_w,
        vtan_casing_r,
        # duct outer surface (nacelle)
        vtan_nacelle,
        vtan_nacelle_b,
        nacelle_jump_term,
        vtan_nacelle_w,
        vtan_nacelle_r,
        # center body
        vtan_centerbody,
        vtan_centerbody_b,
        centerbody_jump_term,
        vtan_centerbody_w,
        vtan_centerbody_r,
    )
end

"""
Calculate tangential velocity for a given net circulation and radial location
"""
function calculate_vtheta(Gamma_tilde, r)
    T = promote_type(eltype(Gamma_tilde), eltype(r))
    vtheta = zeros(eltype(T), length(r))

    for (i, (gti, ri)) in enumerate(zip(Gamma_tilde, r))
        if isapprox(ri, 0.0)
            vtheta[i] = 0.0
        else
            vtheta[i] = gti ./ (2.0 * pi * ri)
        end
    end

    return vtheta
end

"""
Calculate the induced velocities on one of the body wakes (unit velocity inputs determine which one)
"""
function calculate_induced_velocities_on_bodywake(
    vz_w, vr_w, gamw, vz_r, vr_r, sigr, vz_b, vr_b, gamb, Vinf
)

    # problem dimensions
    nrotor = size(sigr, 2) # number of rotors
    np = size(vz_b, 1) # number of panels in bodywake

    # initialize outputs
    vz = Vinf * ones(eltype(gamw), np) # axial induced velocity
    vr = zeros(eltype(gamw), np) # radial induced velocity

    # add body induced velocities
    @views vz[:] .+= vz_b * gamb
    @views vr[:] .+= vr_b * gamb

    # add wake induced velocities
    @views vz[:] .+= vz_w * gamw
    @views vr[:] .+= vr_w * gamw

    # add rotor induced velocities
    for jrotor in 1:nrotor
        @views vz[:] .+= vz_r[jrotor] * sigr[:, jrotor]
        @views vr[:] .+= vr_r[jrotor] * sigr[:, jrotor]
    end

    # return raw induced velocities
    return vz, vr
end

######################################################################
#                                                                    #
#                        Pressure Functions                          #
#                                                                    #
######################################################################

"""
Calculate steady pressure coefficient
"""
function steady_cp(vs, vinf, vref)
    return (vinf^2 .- vs .^ 2) / vref^2
end

"""
only used in post-process for cp.
expression not in dfdc theory, comes from source code.
"""
function calculate_entropy_jumps(sigr, Wm_rotor)
    # average sigr's
    sigr_avg = similar(Wm_rotor, size(sigr,1)-1, size(sigr,2)) .=0
    for (i,s) in enumerate(eachcol(sigr))
        sigr_avg[:,i] = (s[2:end]+s[1:end-1])/2.0
    end

    # multiply by Wm_rotor's
    return sigr_avg .* Wm_rotor
end


"""
Calculate change in pressure coefficient aft of rotor, due to rotor
"""
function delta_cp(deltaH, deltaS, Vtheta, Vref)
    if isapprox(Vref, 0.0)
        return 0.0
    else
        return (2.0 * (deltaH - deltaS) .- Vtheta .^ 2) / Vref^2
    end
end

"""
Calculate net circulation and enthalpy and entropy disk jumps
"""
function calculate_rotor_jumps(Gamr, Omega, B, sigr, Wm_rotor)

    # - Calculate net circulations - #
    Gamma_tilde = calculate_net_circulation(Gamr, B)

    # - Calculate enthalpy disk jump - #
    Htilde = calculate_enthalpy_jumps(Gamr, Omega, B)

    # - Calculate entropy disk jump - #
    Stilde = calculate_entropy_jumps(sigr, Wm_rotor)

    return Gamma_tilde, Htilde, Stilde
end

"""
Calculate change in pressure coefficient due to rotors specifically on the body panels aft of the rotors
"""
function calculate_body_delta_cp!(
    cp, Gamr, sigr, Wm_rotor, Vref, Omega, B, body_vortex_panels, didr, hidr
)

    ## -- Calculate change in pressure coefficient -- ##

    Gamma_tilde, Htilde, Stilde = calculate_rotor_jumps(Gamr, Omega, B, sigr, Wm_rotor)

    nrotor = size(Gamr, 2)

    for ir in 1:nrotor

        # - Get the tangential velocities on the bodies - #
        v_theta_duct = calculate_vtheta(
            Gamma_tilde[end, ir], body_vortex_panels.controlpoint[2, didr[ir]]
        )
        v_theta_hub = calculate_vtheta(
            Gamma_tilde[1, ir], body_vortex_panels.controlpoint[2, hidr[ir]]
        )

        # assemble change in cp due to enthalpy and entropy behind rotor(s)
        cp[didr[ir]] .+= delta_cp(Htilde[end, ir], Stilde[end, ir], v_theta_duct, Vref)
        cp[hidr[ir]] .+= delta_cp(Htilde[1, ir], Stilde[1, ir], v_theta_hub, Vref)
    end

    return nothing
end

"""
Calculate change in pressure coefficient due to rotors specifically on the body wakes
"""
function calculate_bodywake_delta_cp(Gamr, sigr, Wm_rotor, Vref, Omega, B, r; body="duct")

    ## -- Calculate change in pressure coefficient -- ##

    Gamma_tilde, Htilde, Stilde = calculate_rotor_jumps(Gamr, Omega, B, sigr, Wm_rotor)

    # - Get the tangential velocities on the bodies - #
    if body == "duct"
        gt = Gamma_tilde[end]
        ht = Htilde[end]
        st = Stilde[end]
    else
        gt = Gamma_tilde[1, end]
        ht = Htilde[1, end]
        st = Stilde[1, end]
    end

    v_theta_wake = calculate_vtheta(gt, r)

    # assemble change in cp due to enthalpy and entropy behind rotor(s)
    deltacp = delta_cp(ht, st, v_theta_wake, Vref)

    return deltacp
end

"""
calculate pressure coefficient distribution on duct/hub walls
formulation taken from DFDC source code. TODO: derive where the expressions came from.

"""
function get_body_cps(
    Vs, Gamr, sigr, Wm_rotor, Vinf, Vref, B, Omega, didr, hidr, body_vortex_panels, isduct
)

    # - Calculate standard pressure coefficient expression - #
    cp = steady_cp(Vs, Vinf, Vref)

    # - add the change in Cp on the walls due to enthalpy, entropy, and vtheta - #
    calculate_body_delta_cp!(
        cp, Gamr, sigr, Wm_rotor, Vref, Omega, B, body_vortex_panels, didr, hidr
    )

    # - Split body strengths into inner/outer duct and hub - #
    cpdi, cpdo, cph, xdi, xdo, xh = split_bodies(cp, body_vortex_panels; duct=isduct)

    return cp, cpdi, cpdo, cph, body_vortex_panels.controlpoint[:, 1], xdi, xdo, xh
end

"""
Calculate the pressure coefficient distributions on one of the body wakes
"""
function get_bodywake_cps(
    Gamr,
    vz_w,
    vr_w,
    gamw,
    vz_r,
    vr_r,
    sigr,
    vz_b,
    vr_b,
    gamb,
    panels,
    Wm_rotor,
    Omega,
    B,
    Vinf,
    Vref;
    body="duct",
)

    # - Get "surface" velocities - #

    # get induced velocities
    vz_bodywake, vr_bodywake = calculate_induced_velocities_on_bodywake(
        vz_w, vr_w, gamw, vz_r, vr_r, sigr, vz_b, vr_b, gamb, Vinf
    )

    # get "surface" velocities
    Vmat = [vz_bodywake vr_bodywake]
    vtan = [dot(v, t) for (v, t) in zip(eachrow(Vmat), panels.tangent)]

    # - Get steady pressure coefficients - #
    cp_steady = steady_cp(vtan, Vinf, Vref)

    # - Get delta cp - #
    deltacp = calculate_bodywake_delta_cp(
        Gamr, sigr, Wm_rotor, Vref, Omega, B, panels.controlpoint[2, :]; body=body
    )

    return cp_steady .+ deltacp, vtan
end

"""
Calculate dimensional and non-dimensional axial force on a single body
"""
function forces_from_pressure(cps, panels; rhoinf=1.225, Vref=1.0)

    # - rename for convenience - #
    #just want x-component of normals since it's axisymmetric
    ns = panels.normal[1,:]
    #radial positions
    rs = panels.controlpoint[2, :]
    #panel lengths
    ds = panels.influence_length
    # dimensions
    np = length(cps)

    # - initialize - #
    cfx = 0.0 # axial force coefficient (all others are zero for axisymmetric case)
    # - rectangular integration due to constant panel strengths. - #
    for i in 1:np
        cfx += cps[i] * ns[i] * ds[i] * 2.0 * pi * rs[i]
    end

    #dimensionalize
    q = 0.5 * rhoinf * Vref^2

    #note, thrust is in negative x-direction
    return cfx * q, cfx
end

######################################################################
#                                                                    #
#                       Rotor Aero Performance                       #
#                                                                    #
######################################################################

function inviscid_rotor_trust(Wtheta, Gamma_tilde, rotor_panel_length, rhoinf)

    # problem dimensions
    nr, nrotor = size(Gamma_tilde)

    # initialize
    dTi = similar(Gamma_tilde) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            # section thrust
            dTi[ir, irotor] =
                -rhoinf *
                Gamma_tilde[ir, irotor] *
                Wtheta[ir, irotor] *
                rotor_panel_length[ir, irotor]
        end
    end

    #sum the section thrust
    Tinv = sum(dTi; dims=1)

    return Tinv, dTi
end

function viscous_rotor_thrust(
    Wz_rotor, Wmag_rotor, B, chord, rotor_panel_length, cd, rhoinf
)

    # get dimensions
    nr, nrotor = size(Wz_rotor)

    #initialize
    dTv = similar(Wz_rotor) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            hrwc = 0.5 * rhoinf * Wmag_rotor[ir, irotor] * chord[ir, irotor]
            bdr = B[irotor] * rotor_panel_length[ir, irotor]
            dTv[ir, irotor] = -hrwc * cd[ir, irotor] * Wz_rotor[ir, irotor] * bdr
        end
    end

    Tvisc = sum(dTv; dims=1)

    return Tvisc, dTv
end

function inviscid_rotor_torque(
    Wz_rotor, rotor_panel_center, rotor_panel_length, Gamma_tilde, rhoinf
)

    # dimensions
    nr, nrotor = size(Gamma_tilde)

    # initialize
    dQi = similar(Gamma_tilde) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            rdr = rotor_panel_center[ir, irotor] * rotor_panel_length[ir, irotor]
            dQi[ir, irotor] = rhoinf * Gamma_tilde[ir, irotor] * Wz_rotor[ir, irotor] * rdr
        end
    end

    Qinv = sum(dQi; dims=1)

    return Qinv, dQi
end

function viscous_rotor_torque(
    Wtheta_rotor, Wmag_rotor, B, chord, rotor_panel_center, rotor_panel_length, cd, rhoinf
)

    # dimensions
    nr, nrotor = size(Wtheta_rotor)

    # initialize
    dQv = similar(chord) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            hrwc = 0.5 * rhoinf * Wmag_rotor[ir, irotor] * chord[ir, irotor]
            brdr =
                B[irotor] * rotor_panel_center[ir, irotor] * rotor_panel_length[ir, irotor]
            dQv[ir, irotor] = -hrwc * cd[ir, irotor] * Wtheta_rotor[ir, irotor] * brdr
        end
    end

    Qvisc = sum(dQv; dims=1)

    return Qvisc, dQv
end

function inviscid_rotor_power(Qinv, Omega)
    return Qinv .* Omega
end

function viscous_rotor_power(Qvisc, Omega)
    return Qvisc .* Omega
end

function get_total_efficiency(total_thrust, total_power, Vinf)
    eta = zeros(eltype(total_thrust), length(total_thrust))
    for i in 1:length(total_thrust)
        if Vinf <= 0.0 || total_power[i] <= 0.0 || total_thrust[i] <= 0.0
            #do nothing
        else
            eta[i] = total_thrust[i] * Vinf / total_power[i]
        end
    end
    return eta
end

function get_induced_efficiency(Tinv, Tduct, Pinv, Vinf)
    if Vinf <= 0.0 || Pinv <= 0.0
        return 0.0
    else
        return Vinf * (Tinv .+ Tduct) ./ Pinv
    end
end

function get_ideal_efficiency(total_thrust, rhoinf, Vinf, Rref)
    if Vinf != 0.0
        TC = total_thrust / (0.5 * rhoinf * Vinf^2 * pi * Rref^2)
        return 2.0 / (1.0 + sqrt(max(TC, -1.0) + 1.0))
    else
        return 0.0
    end
end

function tqpcoeff(thrust, torque, rhoinf, Omega, Rref)
    T = promote_type(eltype(thrust), eltype(torque), eltype(Omega))
    CT = zeros(T, length(Omega))
    CQ = zeros(T, length(Omega))
    CP = zeros(T, length(Omega))

    for (i, o) in enumerate(Omega)
        if isapprox(o, 0.0)
            CT[i] = CQ[i] = CP[i] = 0.0
        else
            # reference diameter
            D = 2.0 * Rref

            # rototion in rev per second
            n = o / (2.0 * pi)

            # thrust coefficient
            CT[i] = thrust[i] / (rhoinf * n^2 * D^4)

            # torque coefficient
            CQ[i] = torque[i] / (rhoinf * n^2 * D^5)

            # power coefficient
            CP[i] = CQ[i] * o
        end
    end

    return CT, CQ, CP
end

function get_blade_aero(
    Wmag_rotor,
    Wm_rotor,
    Wtheta_rotor,
    solidity,
    stagger,
    chord,
    twist,
    afparamsin,
    afparamsout,
    innerfrac,
    rhoinf,
    muinf,
    asound,
)

    # dimensions
    nr, nrotor = size(Wmag_rotor)

    # initialize
    cl = similar(Wmag_rotor) .= 0.0
    cd = similar(Wmag_rotor) .= 0.0
    phi = similar(Wmag_rotor) .= 0.0
    alpha = similar(Wmag_rotor) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            reynolds = chord[ir, irotor] * abs(Wmag_rotor[ir, irotor]) * rhoinf / muinf

            phi[ir, irotor], alpha[ir, irotor] = calculate_inflow_angles(
                Wm_rotor[ir, irotor], Wtheta_rotor[ir, irotor], twist[ir, irotor]
            )

            #get inner values
            clin, cdin, _ = c4b.dfdceval(
                Wmag_rotor[ir, irotor],
                reynolds,
                solidity[ir, irotor],
                stagger[ir, irotor],
                alpha[ir, irotor],
                afparamsin[ir, irotor],
                asound,
            )
            # get outer values
            clout, cdout, _ = c4b.dfdceval(
                Wmag_rotor[ir, irotor],
                reynolds,
                solidity[ir, irotor],
                stagger[ir, irotor],
                alpha[ir, irotor],
                afparamsout[ir, irotor],
                asound,
            )

            # interpolate inner and outer values
            cl[ir, irotor] = fm.linear([0.0; 1.0], [clin, clout], innerfrac[ir, irotor])
            cd[ir, irotor] = fm.linear([0.0; 1.0], [cdin, cdout], innerfrac[ir, irotor])
        end
    end

    return cl, cd, phi, alpha
end

function get_blade_loads(Wmag_rotor, phi, cl, cd, chords, rhoinf)#, Rhub, Rtip, rpc,rpl ,B, Omega)

    # dimensions
    nr, nrotor = size(Wmag_rotor)

    # initialize
    Np = similar(Wmag_rotor) .= 0.0
    Tp = similar(Wmag_rotor) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            # rename for convenience
            cphi = cos(phi[ir, irotor])
            sphi = sin(phi[ir, irotor])

            # resolve lift and drag into normal and tangential coefficients
            cn = cl[ir, irotor] * cphi - cd[ir, irotor] * sphi
            ct = cl[ir, irotor] * sphi + cd[ir, irotor] * cphi

            # get the normal and tangential loads per unit length N' and T'
            Np[ir, irotor] =
                cn * 0.5 * rhoinf * Wmag_rotor[ir, irotor]^2 * chords[ir, irotor]
            Tp[ir, irotor] =
                ct * 0.5 * rhoinf * Wmag_rotor[ir, irotor]^2 * chords[ir, irotor]
        end
    end

    # Npfull = [zeros(nrotor)'; Np; zeros(nrotor)']
    # Tpfull = [zeros(nrotor)'; Tp; zeros(nrotor)']

    # TODO: consider comparing this with DFDC versions for them

    ## -- Integrate Loads to get Thrust and Torque
    # add hub/tip for complete integration.  loads go to zero at hub/tip.
    # rfull = [Rhub; rpc; Rtip]

    # thrust and torqe distributions
    # thrust = Npfull
    # torque = Tpfull .* rfull

    # integrate Thrust and Torque (trapezoidal)
    # T = B * fm.trapz(rfull, thrust)
    # Q = B * fm.trapz(rfull, torque)
    # - Actually use rectangle rather than trapezoid integration
    # T = B * sum(rpl.*Np)
    # Q = B * sum(rpl.* Tp.*rpc)
    # P = Q * Omega

    return Np, Tp
end

######################################################################
#                                                                    #
#                         Intermediate Values                        #
#                                                                    #
######################################################################
"""
"""
function get_intermediate_values(Gamr, sigr, gamw, inputs)

    # - Extract commonly used items from precomputed inputs - #
    blade_elements = inputs.blade_elements
    rpc = inputs.rotor_panel_centers
    Vinf = inputs.freestream.Vinf
    freestream = inputs.freestream

    calculate_body_vortex_strengths!(
        inputs.gamb,
        inputs.A_bb,
        inputs.b_bf,
        gamw,
        inputs.A_bw,
        sigr,
        inputs.A_br,
        inputs.RHS;
        post=false,
    )

    vz_rotor, vr_rotor, vtheta_rotor, vz_rotor_b, vr_rotor_b, vz_rotor_w, vr_rotor_w, vz_rotor_r, vr_rotor_r = calculate_induced_velocities_on_rotors(
        blade_elements,
        Gamr,
        inputs.vz_rw,
        inputs.vr_rw,
        gamw,
        inputs.vz_rr,
        inputs.vr_rr,
        sigr,
        inputs.vz_rb,
        inputs.vr_rb,
        inputs.gamb[1:(inputs.body_vortex_panels.totnode)];
        post=true,
    )

    Wz_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor = reframe_rotor_velocities(
        vz_rotor, vr_rotor, vtheta_rotor, Vinf, inputs.blade_elements.Omega, rpc
    )

    Gamma_tilde = calculate_net_circulation(Gamr, blade_elements.B)

    H_tilde = calculate_enthalpy_jumps(Gamr, blade_elements.Omega, blade_elements.B)

    cl, cd, phi, alpha = update_coeffs(
        blade_elements,
        Wm_rotor,
        Wtheta_rotor,
        Wmag_rotor,
        freestream;
        post=true,
        verbose=false,
    )

    vz_wake, vr_wake, vz_wake_b, vr_wake_b, vz_wake_r, vr_wake_r, vz_wake_w, vr_wake_w = calculate_induced_velocities_on_wakes(
        inputs.vz_ww,
        inputs.vr_ww,
        gamw,
        inputs.vz_wr,
        inputs.vr_wr,
        sigr,
        inputs.vz_wb,
        inputs.vr_wb,
        inputs.gamb[1:(inputs.body_vortex_panels.totnode)];
        post=true,
    )

    Wz_wake, Wm_wake = reframe_wake_velocities(vz_wake, vr_wake, Vinf; post=true)

    return (;
        # velocities at rotor
        vz_rotor,
        vr_rotor,
        vtheta_rotor,
        Wz_rotor,
        Wtheta_rotor,
        Wm_rotor,
        Wmag_rotor,
        # raw induced velocities on rotor
        vz_rotor_b,
        vr_rotor_b,
        vz_rotor_w,
        vr_rotor_w,
        vz_rotor_r,
        vr_rotor_r,
        # velocities at wake
        vz_wake,
        vr_wake,
        Wz_wake,
        Wm_wake,
        # raw induced velocities on wake
        vz_wake_b,
        vr_wake_b,
        vz_wake_w,
        vr_wake_w,
        vz_wake_r,
        vr_wake_r,
        # tilde jumps on blades
        Gamma_tilde,
        H_tilde,
        # blade element angles
        phi,
        alpha,
        # blade element forces
        cl,
        cd,
    )
end

######################################################################
#                                                                    #
#                        Post-Post-Processing                        #
#                                                                    #
######################################################################

"""
probe_poses : matrix of x,r coordinates of locations at which to probe velocity field
"""
function probe_velocity_field(probe_poses, inputs, states; debug=false)

    # - Types - #
    TF = promote_type(eltype(probe_poses), eltype(states))

    # - rename for convenience - #
    (; rotor_source_panels, wake_vortex_panels, body_vortex_panels, blade_elements) =
        inputs

    # get number of rotor blades
    num_blades = blade_elements.B

    # - extract states - #
    gamb, gamw, Gamr, sigr = extract_state_variables(states, inputs)

    # - dimensions - #
    nv = size(probe_poses, 1)
    nr, nrotor = size(Gamr)

    # - initialize - #
    Vzr = zeros(TF, nv, 2)
    Vb = zeros(TF, nv, 2)
    Vr = zeros(TF, nv, 2)
    Vw = zeros(TF, nv, 2)
    vtheta = zeros(TF, nv)

    ##### ----- Axial and Radial Velocity ----- #####
    vfromdoubletpanels!(Vb, probe_poses, body_vortex_panels.nodes, gamb)
    for i in 1:length(rotor_source_panels)
        vfromsourcepanels!(
            Vr,
            probe_poses,
            rotor_source_panels[i].controlpoint,
            rotor_source_panels[i].len,
            sigr[:, i],
        )
    end
    vfromvortexpanels!(
        Vw, probe_poses, wake_vortex_panels.controlpoint, wake_vortex_panels.len, gamw
    )

    Vzr .+= Vb .+ Vr .+ Vw

    ###### ----- Tangential Velocity ----- #####
    # TODO: there is something here that is making things jump around.
    # TODO: thinking about it more, probably can't try and make it continuous
    # TODO: comment all of this out and just find the Gamma_tilde for the associated blade element and use the probe radial position.  I don't think there's a way to properly smooth this out.

    # reshape the wake panel control points into the wake sheets
    nsheets = size(Gamr, 1) + 1
    nwakex = Int(wake_vortex_panels.totpanel / nsheets)
    wakecpx = reshape(wake_vortex_panels.controlpoint[:, 1], (nwakex, nsheets))'
    wakecpr = reshape(wake_vortex_panels.controlpoint[:, 2], (nwakex, nsheets))'
    rotorzloc = inputs.blade_elements.rotorzloc

    # get B*Circulation on each rotor at each wake shedding location
    if size(Gamr, 1) > 1
        Gambar = [
            Gamr[1, :]'
            (Gamr[2:end, :] .+ Gamr[1:(end - 1), :]) ./ 2
            Gamr[end, :]'
        ]
    else
        Gambar = [
            Gamr[1, :]
            Gamr[end, :]
        ]
    end
    BGambar = reshape(Gambar, (nsheets, nrotor)) .* num_blades'

    # Get Gamma_tilde on rotors
    Gamma_tilde_rotor = similar(BGambar) .= 0.0
    for irotor in 1:nrotor
        Gamma_tilde_rotor[:, irotor] .+= 0.5 * BGambar[:, irotor]
        Gamma_tilde_rotor[:, (irotor + 1):end] .+= BGambar[:, irotor]
    end

    # Get Gamma_tilde at each wake control point
    Gamma_tilde = cumsum(BGambar; dims=2)
    Gamma_tilde_grid = similar(wakecpx) .= 0.0
    xids = [searchsortedfirst(wakecpx[1, :], rotorzloc[i]) for i in 1:nrotor]
    for i in 1:nrotor
        if i == nrotor
            Gamma_tilde_grid[:, xids[i]:end] .= Gamma_tilde[:, i]
        else
            Gamma_tilde_grid[:, xids[i]:(xids[i + 1] - 1)] .= Gamma_tilde[:, i]
        end
    end

    for (ip, probe) in enumerate(eachrow(probe_poses))

        # - Find wake stations just in front, and just behind - #
        wzid1 = findlast(x -> x <= probe[1], wakecpx[1, :])
        wzid2 = findfirst(x -> x >= probe[1], wakecpx[1, :])

        # check if outside the wake or on the edge
        if isnothing(wzid1)
            if probe[1] < rotorzloc[1]
                #outside of wake
                vtheta[ip] = 0.0
                continue
            else
                wzid1 = wzid2
            end

        elseif isnothing(wzid2)
            if probe[1] > wakecpx[end, end]
                #outside of wake
                vtheta[ip] = 0.0
                continue
            else
                wzid2 = wzid1
            end
        end

        # - Find the rotor just in front - #
        xrid = findlast(x -> x <= probe[1], rotorzloc)

        # - Find the wake station just below and just above - #
        wrid1 = findlast(x -> x <= probe[2], wakecpr[:, wzid1])
        wrid2 = findfirst(x -> x >= probe[2], wakecpr[:, wzid2])

        # Check if outside the wake, or on the edge
        if isnothing(wrid1)
            if probe[2] < (wakecpr[wzid1, 1] + wakecpr[wzid2, 1]) / 2
                #outside of wake
                vtheta[ip] = 0.0
                continue
            else
                wrid1 = wrid2
            end

        elseif isnothing(wrid2)
            if probe[2] > (wakecpr[end, wzid1] + wakecpr[end, wzid2]) / 2
                #outside of wake
                vtheta[ip] = 0.0
                continue
            else
                wrid2 = wrid1
            end
        end

        # check if on rotor or aligned with wake control points
        if isapprox(probe[1], rotorzloc[xrid])
            # On the a rotor need to use self-induced rotor tangential velocity and interpolate in r only

            #on edges, just use the edge value
            if wrid1 == wrid2
                vtheta[ip] =
                    Gamma_tilde_rotor[wrid2, xrid] ./ (2 * pi .* wakecpr[wrid2, wzid1])
            else
                vthetaself1 = Gamma_tilde_rotor[wrid1, xrid]# ./ (2 * pi * wakecpr[wrid1, wzid1])
                vthetaself2 = Gamma_tilde_rotor[wrid2, xrid]# ./ (2 * pi * wakecpr[wrid2, wzid1])
                vtheta[ip] =
                    fm.linear(
                        [wakecpr[wrid1, wzid1]; wakecpr[wrid2, wzid1]],
                        [vthetaself1; vthetaself2],
                        probe[2],
                    ) ./ (2 * pi * probe[2])
            end

        elseif wzid1 == wzid2
            # inline with wake control points need to interpolate in r only

            #on edges, just use the edge value
            if wrid1 == wrid2
                vtheta[ip] =
                    Gamma_tilde_grid[wrid1, wzid1] ./ (2 * pi * wakecpr[wrid1, wzid1])
            else
                vtheta1 = Gamma_tilde_grid[wrid1, wzid1]# ./ (2 * pi * wakecpr[wrid1, wzid1])
                vtheta2 = Gamma_tilde_grid[wrid2, wzid1]#./ (2 * pi * wakecpr[wrid2, wzid1])
                vtheta[ip] =
                    fm.linear(
                        [wakecpr[wrid1, wzid1]; wakecpr[wrid2, wzid1]],
                        [vtheta1; vtheta2],
                        probe[2],
                    ) ./ (2 * pi * probe[2])
            end

        elseif wrid1 == wrid2
            #wake radius doesn't change, just use first radial point
            vtheta[ip] = Gamma_tilde_grid[wrid1, wzid1] ./ (2 * pi * wakecpr[wrid1, wzid1])

        else

            # use rotor self induced value if it's closer than the nearest wake control point
            if probe[1] < rotorzloc[xrid]
                vtx2r1 = Gamma_tilde_rotor[wrid1, xrid]# ./ (2 * pi * wakecpr[wrid1])
                vtx2x2 = Gamma_tilde_rotor[wrid2, xrid]# ./ (2 * pi * wakecpr[wrid2])
            else
                vtx2r1 = Gamma_tilde_grid[wrid1, wzid2]# ./ (2 * pi * wakecpr[wrid1])
                vtx2x2 = Gamma_tilde_grid[wrid2, wzid2]# ./ (2 * pi * wakecpr[wrid2])
            end

            if probe[1] > rotorzloc[xrid]
                vtx1r1 = Gamma_tilde_rotor[wrid1, xrid]# ./ (2 * pi * wakecpr[wrid1])
                vtx1r2 = Gamma_tilde_rotor[wrid2, xrid]# ./ (2 * pi * wakecpr[wrid2])
            else
                vtx1r1 = Gamma_tilde_grid[wrid1, wzid1]# ./ (2 * pi * wakecpr[wrid1])
                vtx1r2 = Gamma_tilde_grid[wrid2, wzid1]# ./ (2 * pi * wakecpr[wrid2])
            end

            r1 = fm.linear(
                [wakecpx[wrid1, wzid1]; wakecpx[wrid1, wzid2]],
                [wakecpr[wrid1, wzid1]; wakecpr[wrid1, wzid2]],
                probe[1],
            )
            r2 = fm.linear(
                [wakecpx[wrid2, wzid1]; wakecpx[wrid2, wzid2]],
                [wakecpr[wrid2, wzid1]; wakecpr[wrid2, wzid2]],
                probe[1],
            )

            vt1 = fm.linear(
                [wakecpx[wrid1, wzid1]; wakecpx[wrid1, wzid2]], [vtx1r1; vtx2r1], probe[1]
            )
            vt2 = fm.linear(
                [wakecpx[wrid2, wzid1]; wakecpx[wrid2, wzid2]], [vtx1r2; vtx2x2], probe[1]
            )

            vtheta[ip] = fm.linear([r1; r2], [vt1; vt2], probe[2]) ./ (2 * pi * probe[2])
        end
    end

    if debug
        return Vzr[:, 1],
        Vzr[:, 2], vtheta, Vb[:, 1], Vb[:, 2], Vr[:, 1], Vr[:, 2], Vw[:, 1],
        Vw[:, 2]
    else
        return Vzr[:, 1], Vzr[:, 2], vtheta
    end
end
