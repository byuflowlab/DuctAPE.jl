#=

Various Post-processing functions

=#

function post_process(states, inputs)

    # - things contained in iv tuple
    # mub,
    # gamw,
    # Gamr,
    # sigr,
    # vx_rotor,
    # vr_rotor,
    # vtheta_rotor,
    # Wx_rotor,
    # Wtheta_rotor,
    # Wm_rotor,
    # Wmag_rotor,
    # vxfrombody,
    # vrfrombody,
    # vxfromwake,
    # vrfromwake,
    # vxfromrotor,
    # vrfromrotor,
    # Gamma_tilde,
    # H_tilde,
    # phi,
    # alpha,
    # cl,
    # cd,

    iv = get_intermediate_values(states, inputs) #intermediate values

    # get problem dimensions
    nr, nrotor = size(iv.Gamr)
    nw = nr + 1

    # - Extract convenient input fields - #
    Vinf = inputs.Vinf
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
    rpl = reshape(reduce(vcat, (p -> p.len).(inputs.rotor_source_panels)), (nr, nrotor))

    ref = (; Vinf, Vref, Rref, rhoinf, muinf, asound, Omega, B)

    # get blade element coefficients
    clift, cdrag, inflow_angle, angle_of_attack = get_blade_aero(
        iv.Wmag_rotor,
        iv.Wm_rotor,
        iv.Wtheta_rotor,
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

    # - Rotor Thrust - #
    # inviscid thrust
    rotor_inviscid_thrust, rotor_inviscid_thrust_dist = inviscid_rotor_trust(
        iv.Wtheta_rotor, iv.Gamma_tilde, rpl, rhoinf
    )

    # viscous thrust
    rotor_viscous_thrust, rotor_viscous_thrust_dist = viscous_rotor_thrust(
        iv.Wx_rotor, iv.Wmag_rotor, B, chord, rpl, iv.cd, rhoinf
    )

    # total thrust
    rotor_thrust = rotor_inviscid_thrust .+ rotor_viscous_thrust

    # - Rotor Torque - #
    # inviscid torque
    rotor_inviscid_torque, rotor_inviscid_torque_dist = inviscid_rotor_torque(
        iv.Wx_rotor, rpc, rpl, iv.Gamma_tilde, rhoinf
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
    # body_surface_velocity, duct_inner_vs, duct_outer_vs, hub_vs, vtot_on_body, vsfromvinf, vsfrombodypanels, vsfromTEpanels, vsfromgradmu, vsfromwake, vsfromrotors
    bodyv = get_body_vs(iv.mub, iv.gamw, iv.sigr, inputs)

    (;
        vtot_on_body,
        body_surface_velocity,
        duct_inner_vs,
        duct_outer_vs,
        hub_vs,
        vsfromvinf,
        vsfrombodypanels,
        duct_inner_vs_from_body,
        duct_outer_vs_from_body,
        hub_vs_from_body,
        vsfromTEpanels,
        duct_inner_vs_from_TE,
        duct_outer_vs_from_TE,
        hub_vs_from_TE,
        vsfromgradmu,
        duct_inner_vs_from_gradmu,
        duct_outer_vs_from_gradmu,
        hub_vs_from_gradmu,
        vsfromwake,
        duct_inner_vs_from_wake,
        duct_outer_vs_from_wake,
        hub_vs_from_wake,
        vsfromrotors,
        duct_inner_vs_from_rotor,
        duct_outer_vs_from_rotor,
        hub_vs_from_rotor,
    ) = bodyv

    ## -- Pressure on Bodies -- ##
    body_cp, duct_inner_cp, duct_outer_cp, hub_cp, body_x, duct_inner_x, duct_outer_x, hub_x = get_body_cps(
        body_surface_velocity,
        iv.Gamr,
        iv.sigr,
        iv.Wm_rotor,
        Vinf,
        Vref,
        B,
        Omega,
        didr,
        hidr,
        inputs.body_doublet_panels,
        inputs.isduct,
    )

    ## -- Pressure on Body Wakes -- ##
    ductwake_cp, ductwake_vs = get_bodywake_cps(
        iv.Gamr,
        inputs.vx_dww,
        inputs.vr_dww,
        iv.gamw,
        inputs.vx_dwr,
        inputs.vr_dwr,
        iv.sigr,
        inputs.vx_dwb,
        inputs.vr_dwb,
        iv.mub,
        inputs.vx_dwbte,
        inputs.vr_dwbte,
        (p -> p.idx).(inputs.body_doublet_panels.TEnodes),
        inputs.duct_wake_panels,
        iv.Wm_rotor,
        Omega,
        B,
        Vinf,
        Vref;
        body="duct",
    )

    hubwake_cp, hubwake_vs = get_bodywake_cps(
        iv.Gamr,
        inputs.vx_hww,
        inputs.vr_hww,
        iv.gamw,
        inputs.vx_hwr,
        inputs.vr_hwr,
        iv.sigr,
        inputs.vx_hwb,
        inputs.vr_hwb,
        iv.mub,
        inputs.vx_hwbte,
        inputs.vr_hwbte,
        (p -> p.idx).(inputs.body_doublet_panels.TEnodes),
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
    duct_cp = [duct_inner_cp; duct_outer_cp]

    # - Calculate Thrust from Bodies - #

    body_thrust, _ = forces_from_pressure(
        duct_cp, inputs.body_doublet_panels; rhoinf=rhoinf, Vref=Vref
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
    CT, CQ = tqcoeff(total_thrust, total_torque, rhoinf, Omega, Rref)

    ## -- Assemble Output Tuple -- ##

    return (;
        # - Geometry - #
        Rref,
        rotor_panel_center=rpc,
        chord,
        twist,
        inputs.duct_coordinates,
        inputs.hub_coordinates,
        inputs.body_doublet_panels,
        inputs.rotor_source_panels,
        inputs.wake_vortex_panels,
        # - Intermediate Values - #
        intermediate_values=iv,
        # - Reference Values - #
        reference_values=ref,
        # - States - #
        mub=iv.mub,
        gamw=iv.gamw,
        Gamr=iv.Gamr,
        sigr=iv.sigr,
        # - Body Values - #
        # body thrust
        body_thrust,
        # surface velocities and pressures
        # duct_inner_vs,
        duct_inner_cp,
        duct_inner_x,
        # duct_outer_vs,
        duct_outer_cp,
        duct_outer_x,
        # hub_vs,
        hub_cp,
        hub_x,
        #individual body velocity contributions
        vtot_on_body,
        body_surface_velocity,
        duct_inner_vs,
        duct_outer_vs,
        hub_vs,
        vsfromvinf,
        vsfrombodypanels,
        duct_inner_vs_from_body,
        duct_outer_vs_from_body,
        hub_vs_from_body,
        vsfromTEpanels,
        duct_inner_vs_from_TE,
        duct_outer_vs_from_TE,
        hub_vs_from_TE,
        vsfromgradmu,
        duct_inner_vs_from_gradmu,
        duct_outer_vs_from_gradmu,
        hub_vs_from_gradmu,
        vsfromwake,
        duct_inner_vs_from_wake,
        duct_outer_vs_from_wake,
        hub_vs_from_wake,
        vsfromrotors,
        duct_inner_vs_from_rotor,
        duct_outer_vs_from_rotor,
        hub_vs_from_rotor,
        # body_surface_velocity,
        # vtot_on_body,
        # vsfromvinf,
        # vsfrombodypanels,
        # vsfromTEpanels,
        # vsfromgradmu,
        # vsfromwake,
        # vsfromrotors,
        # - Body Wake Values - #
        # surface velocities and pressures
        ductwake_vs,
        ductwake_cp,
        ductwake_x=inputs.duct_wake_panels.controlpoint[:, 1],
        hubwake_vs,
        hubwake_cp,
        hubwake_x=inputs.hub_wake_panels.controlpoint[:, 1],
        # - Rotor Values - #
        rotor_efficiency,
        rotor_inviscid_thrust,
        rotor_inviscid_thrust_dist,
        rotor_viscous_thrust,
        rotor_viscous_thrust_dist,
        rotor_thrust,
        CT,
        # rotor torque
        rotor_inviscid_torque,
        rotor_inviscid_torque_dist,
        rotor_viscous_torque,
        rotor_viscous_torque_dist,
        rotor_torque,
        CQ,
        # rotor power
        rotor_inviscid_power,
        rotor_inviscid_power_dist,
        rotor_viscous_power,
        rotor_viscous_power_dist,
        rotor_power,
        # - Blade Values - #
        clift,
        cdrag,
        angle_of_attack,
        inflow_angle,
        blade_normal_force_per_unit_span,
        blade_tangential_force_per_unit_span,
        # - Total Values - #
        total_thrust,
        total_torque,
        total_power,
        total_efficiency=total_efficiency[1],
        induced_efficiency,
        ideal_efficiency,
    )
end

######################################################################
#                                                                    #
#                        Velocity Functions                          #
#                                                                    #
######################################################################

function get_body_vs(mub, gamw, sigr, inputs)
    # - rename for convenience - #
    (; body_doublet_panels, A_br, A_bw) = inputs
    nrotor = size(sigr, 2)

    # - Influence from Freestream - #
    Vinf = inputs.Vinf * [1.0 0.0] # axisymmetric, so no radial component
    vsfromvinf = repeat(Vinf, body_doublet_panels.totpanel) # need velocity on each panel

    ## -- Velocity Contributions from body -- ##

    # - Body-induced Surface Velocity - #
    vsfrombodypanels = vfromdoubletpanels(
        body_doublet_panels.controlpoint, body_doublet_panels.nodes, mub
    )

    # - "Wake"-induced Surface Velocity - #
    vsfromTEpanels = vfromTE(
        body_doublet_panels.controlpoint, body_doublet_panels.TEnodes, mub
    )

    # - ∇μ/2 surface velocity - #
    vsfromgradmu = vfromgradmutry3b(body_doublet_panels, mub)

    ## -- Velocity Contributions from rotors and wake -- ##
    # - wake velocity components - #
    vsfromwake = vfromvortexpanels(
        inputs.body_doublet_panels.controlpoint,
        inputs.wake_vortex_panels.controlpoint,
        inputs.wake_vortex_panels.len,
        gamw,
    )

    # - rotor velocity components - #
    vsfromrotors = similar(vsfromwake) .= 0.0
    for jrotor in 1:nrotor
        vfromsourcepanels!(
            vsfromrotors,
            inputs.body_doublet_panels.controlpoint,
            inputs.rotor_source_panels[jrotor].controlpoint,
            inputs.rotor_source_panels[jrotor].len,
            sigr[:, jrotor],
        )
    end

    ## -- Total Velocity -- ##
    #theoretically, vtot_on_body dot nhat should be zero and vtot_on_body dot that should be norm(vtot_on_body)
    vtot_on_body =
        vsfromvinf .+ vsfrombodypanels .+ vsfromTEpanels .+ vsfromgradmu .+ vsfromwake .+
        vsfromrotors

    # - Get magnitude and split - #
    vs_body = norm.(eachrow(vtot_on_body))

    vsdi, vsdo, vsh, _, _, _ = split_bodies(vs_body, body_doublet_panels; duct=true)
    vbdi, vbdo, vbh, _, _, _ = split_bodies(
        vsfrombodypanels, body_doublet_panels; duct=true
    )
    vtedi, vtedo, vteh, _, _, _ = split_bodies(
        vsfromTEpanels, body_doublet_panels; duct=true
    )
    vgmdi, vgmdo, vgmh, _, _, _ = split_bodies(vsfromgradmu, body_doublet_panels; duct=true)
    vwdi, vwdo, vwh, _, _, _ = split_bodies(vsfromwake, body_doublet_panels; duct=true)
    vrdi, vrdo, vrh, _, _, _ = split_bodies(vsfromrotors, body_doublet_panels; duct=true)

    return (;
        vtot_on_body=vtot_on_body,
        body_surface_velocity=vs_body,
        duct_inner_vs=vsdi,
        duct_outer_vs=vsdo,
        hub_vs=vsh,
        vsfromvinf=vsfromvinf,
        vsfrombodypanels=vsfrombodypanels,
        duct_inner_vs_from_body=vbdi,
        duct_outer_vs_from_body=vbdo,
        hub_vs_from_body=vbh,
        vsfromTEpanels=vsfromTEpanels,
        duct_inner_vs_from_TE=vtedi,
        duct_outer_vs_from_TE=vtedo,
        hub_vs_from_TE=vteh,
        vsfromgradmu=vsfromgradmu,
        duct_inner_vs_from_gradmu=vgmdi,
        duct_outer_vs_from_gradmu=vgmdo,
        hub_vs_from_gradmu=vgmh,
        vsfromwake=vsfromwake,
        duct_inner_vs_from_wake=vwdi,
        duct_outer_vs_from_wake=vwdo,
        hub_vs_from_wake=vwh,
        vsfromrotors=vsfromrotors,
        duct_inner_vs_from_rotor=vrdi,
        duct_outer_vs_from_rotor=vrdo,
        hub_vs_from_rotor=vrh,
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
    vx_w, vr_w, gamw, vx_r, vr_r, sigr, vx_b, vr_b, mub, vx_bte, vr_bte, TEidxs, Vinf
)

    # problem dimensions
    nrotor = size(sigr, 2) # number of rotors
    np = size(vx_b, 1) # number of panels in bodywake

    # initialize outputs
    vx = Vinf * ones(eltype(gamw), np) # axial induced velocity
    vr = zeros(eltype(gamw), np) # radial induced velocity

    # add body induced velocities
    @views vx[:] .+= vx_b * mub
    @views vr[:] .+= vr_b * mub

    # add body TE induced velocities
    @views vx[:] .+= vx_bte * mub[TEidxs]
    @views vr[:] .+= vr_bte * mub[TEidxs]

    # add wake induced velocities
    @views vx[:] .+= vx_w * gamw
    @views vr[:] .+= vr_w * gamw

    # add rotor induced velocities
    for jrotor in 1:nrotor
        @views vx[:] .+= vx_r[jrotor] * sigr[:, jrotor]
        @views vr[:] .+= vr_r[jrotor] * sigr[:, jrotor]
    end

    # return raw induced velocities
    return vx, vr
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
function calculate_rotor_jumps(Gamr, Omega, B, sigr, Vm_rotor)

    # - Calculate net circulations - #
    Gamma_tilde = calculate_net_circulation(Gamr, B)

    # - Calculate enthalpy disk jump - #
    Htilde = calculate_enthalpy_jumps(Gamr, Omega, B)

    # - Calculate entropy disk jump - #
    Stilde = calculate_entropy_jumps(sigr, Vm_rotor)

    return Gamma_tilde, Htilde, Stilde
end

"""
Calculate change in pressure coefficient due to rotors specifically on the body panels aft of the rotors
"""
function calculate_body_delta_cp!(
    cp, Gamr, sigr, Vm_rotor, Vref, Omega, B, body_doublet_panels, didr, hidr
)

    ## -- Calculate change in pressure coefficient -- ##

    Gamma_tilde, Htilde, Stilde = calculate_rotor_jumps(Gamr, Omega, B, sigr, Vm_rotor)

    nrotor = size(Gamr, 2)

    for ir in 1:nrotor

        # - Get the tangential velocities on the bodies - #
        v_theta_duct = calculate_vtheta(
            Gamma_tilde[end, ir], body_doublet_panels.controlpoint[didr[ir], 2]
        )
        v_theta_hub = calculate_vtheta(
            Gamma_tilde[1, ir], body_doublet_panels.controlpoint[hidr[ir], 2]
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
function calculate_bodywake_delta_cp(Gamr, sigr, Vm_rotor, Vref, Omega, B, r; body="duct")

    ## -- Calculate change in pressure coefficient -- ##

    Gamma_tilde, Htilde, Stilde = calculate_rotor_jumps(Gamr, Omega, B, sigr, Vm_rotor)

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
    Vs, Gamr, sigr, Vm_rotor, Vinf, Vref, B, Omega, didr, hidr, body_doublet_panels, isduct
)

    # - Calculate standard pressure coefficient expression - #
    cp = steady_cp(Vs, Vinf, Vref)

    # - add the change in Cp on the walls due to enthalpy, entropy, and vtheta - #
    calculate_body_delta_cp!(
        cp, Gamr, sigr, Vm_rotor, Vref, Omega, B, body_doublet_panels, didr, hidr
    )

    # - Split body strengths into inner/outer duct and hub - #
    cpdi, cpdo, cph, xdi, xdo, xh = split_bodies(cp, body_doublet_panels; duct=isduct)

    return cp, cpdi, cpdo, cph, body_doublet_panels.controlpoint[:, 1], xdi, xdo, xh
end

"""
Calculate the pressure coefficient distributions on one of the body wakes
"""
function get_bodywake_cps(
    Gamr,
    vx_w,
    vr_w,
    gamw,
    vx_r,
    vr_r,
    sigr,
    vx_b,
    vr_b,
    mub,
    vx_bte,
    vr_bte,
    TEidxs,
    panels,
    Vm_rotor,
    Omega,
    B,
    Vinf,
    Vref;
    body="duct",
)

    # - Get "surface" velocities - #

    # get induced velocities
    vx_bodywake, vr_bodywake = calculate_induced_velocities_on_bodywake(
        vx_w, vr_w, gamw, vx_r, vr_r, sigr, vx_b, vr_b, mub, vx_bte, vr_bte, TEidxs, Vinf
    )

    # get "surface" velocities
    Vmat = [vx_bodywake vr_bodywake]
    vs = [dot(v, t) for (v, t) in zip(eachrow(Vmat), panels.tangent)]

    # - Get steady pressure coefficients - #
    cp_steady = steady_cp(vs, Vinf, Vref)

    # - Get delta cp - #
    deltacp = calculate_bodywake_delta_cp(
        Gamr, sigr, Vm_rotor, Vref, Omega, B, panels.controlpoint[:, 2]; body=body
    )

    return cp_steady .+ deltacp, vs
end

"""
Calculate dimensional and non-dimensional axial force on a single body
"""
function forces_from_pressure(cps, panels; rhoinf=1.225, Vref=1.0)

    # - rename for convenience - #
    #just want x-component of normals since it's axisymmetric
    ns = panels.normal[:, 1]
    #radial positions
    rs = panels.controlpoint[:, 2]
    #panel lengths
    ds = panels.len
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
    Wx_rotor, Wmag_rotor, B, chord, rotor_panel_length, cd, rhoinf
)

    # get dimensions
    nr, nrotor = size(Wx_rotor)

    #initialize
    dTv = similar(Wx_rotor) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            hrwc = 0.5 * rhoinf * Wmag_rotor[ir, irotor] * chord[ir, irotor]
            bdr = B[irotor] * rotor_panel_length[ir, irotor]
            dTv[ir, irotor] = -hrwc * cd[ir, irotor] * Wx_rotor[ir, irotor] * bdr
        end
    end

    Tvisc = sum(dTv; dims=1)

    return Tvisc, dTv
end

function inviscid_rotor_torque(
    Wx_rotor, rotor_panel_center, rotor_panel_length, Gamma_tilde, rhoinf
)

    # dimensions
    nr, nrotor = size(Gamma_tilde)

    # initialize
    dQi = similar(Gamma_tilde) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            rdr = rotor_panel_center[ir, irotor] * rotor_panel_length[ir, irotor]
            dQi[ir, irotor] = rhoinf * Gamma_tilde[ir, irotor] * Wx_rotor[ir, irotor] * rdr
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
        if Vinf == 0.0 || total_power[i] <= 0.0 || total_thrust[i] <= 0.0
            #do nothing
        else
            eta[i] = total_thrust[i] * Vinf / total_power[i]
        end
    end
    return eta
end

function get_induced_efficiency(Tinv, Tduct, Pinv, Vinf)
    if Vinf == 0.0 || Pinv == 0.0
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

function tqcoeff(thrust, torque, rhoinf, Omega, Rref)
    T = promote_type(eltype(thrust), eltype(torque), eltype(Omega))
    CT = zeros(T, length(Omega))
    CQ = zeros(T, length(Omega))

    for (i, o) in enumerate(Omega)
        if isapprox(o, 0.0)
            CT[i] = CQ[i] = 0.0
        else
            # reference diameter
            D = 2.0 * Rref

            # rototion in rev per second
            n = o / (2.0 * pi)

            # thrust coefficient
            CT[i] = thrust / (rhoinf * n^2 * D^4)

            # torque coefficient
            CQ[i] = torque / (rhoinf * n^2 * D^5)
        end
    end

    return CT, CQ
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
            clin, cdin, _ = c3b.dfdceval(
                Wmag_rotor[ir, irotor],
                reynolds,
                solidity[ir, irotor],
                stagger[ir, irotor],
                alpha[ir, irotor],
                afparamsin[ir, irotor],
                asound,
            )
            # get outer values
            clout, cdout, _ = c3b.dfdceval(
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
function get_intermediate_values(states, inputs)

    # - Extract commonly used items from precomputed inputs - #
    blade_elements = inputs.blade_elements
    rpc = inputs.rotor_panel_centers
    Vinf = inputs.Vinf

    # - Extract states - #
    mub, gamw, Gamr, sigr = extract_state_variables(states, inputs)

    _, _, _, vxfrombody, vrfrombody, vxfromwake, vrfromwake, vxfromrotor, vrfromrotor = calculate_induced_velocities_on_rotors(
        blade_elements,
        Gamr,
        inputs.vx_rw,
        inputs.vr_rw,
        gamw,
        inputs.vx_rr,
        inputs.vr_rr,
        sigr,
        inputs.vx_rb,
        inputs.vr_rb,
        mub,
        inputs.vx_rbte,
        inputs.vr_rbte,
        (p -> p.idx).(inputs.body_doublet_panels.TEnodes);
        debug=true,
    )

    vx_rotor, vr_rotor, vtheta_rotor, Wx_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor = calculate_rotor_velocities(
        Gamr, gamw, sigr, mub, inputs
    )

    Gamma_tilde = calculate_net_circulation(Gamr, blade_elements.B)

    H_tilde = calculate_enthalpy_jumps(Gamr, blade_elements.Omega, blade_elements.B)

    _, _, phi, alpha, cl, cd = calculate_gamma_sigma(
        blade_elements, Wm_rotor, Wtheta_rotor, Wmag_rotor, inputs.freestream; debug=true
    )

    return (;
        mub,
        gamw,
        Gamr,
        sigr,
        vx_rotor,
        vr_rotor,
        vtheta_rotor,
        Wx_rotor,
        Wtheta_rotor,
        Wm_rotor,
        Wmag_rotor,
        vxfrombody,
        vrfrombody,
        vxfromwake,
        vrfromwake,
        vxfromrotor,
        vrfromrotor,
        Gamma_tilde,
        H_tilde,
        phi,
        alpha,
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
    (; rotor_source_panels, wake_vortex_panels, body_doublet_panels, blade_elements) =
        inputs

    # get number of rotor blades
    num_blades = blade_elements.B

    # - extract states - #
    mub, gamw, Gamr, sigr = extract_state_variables(states, inputs)

    # - dimensions - #
    nv = size(probe_poses, 1)
    nr, nrotor = size(Gamr)

    # - initialize - #
    Vxr = zeros(TF, nv, 2)
    Vb = zeros(TF, nv, 2)
    Vr = zeros(TF, nv, 2)
    Vw = zeros(TF, nv, 2)
    vtheta = zeros(TF, nv)

    ##### ----- Axial and Radial Velocity ----- #####
    vfromdoubletpanels!(Vb, probe_poses, body_doublet_panels.nodes, mub)
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

    Vxr .+= Vb .+ Vr .+ Vw

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
        wxid1 = findlast(x -> x <= probe[1], wakecpx[1, :])
        wxid2 = findfirst(x -> x >= probe[1], wakecpx[1, :])

        # check if outside the wake or on the edge
        if isnothing(wxid1)
            if probe[1] < rotorzloc[1]
                #outside of wake
                vtheta[ip] = 0.0
                continue
            else
                wxid1 = wxid2
            end

        elseif isnothing(wxid2)
            if probe[1] > wakecpx[end, end]
                #outside of wake
                vtheta[ip] = 0.0
                continue
            else
                wxid2 = wxid1
            end
        end

        # - Find the rotor just in front - #
        xrid = findlast(x -> x <= probe[1], rotorzloc)

        # - Find the wake station just below and just above - #
        wrid1 = findlast(x -> x <= probe[2], wakecpr[:, wxid1])
        wrid2 = findfirst(x -> x >= probe[2], wakecpr[:, wxid2])

        # Check if outside the wake, or on the edge
        if isnothing(wrid1)
            if probe[2] < (wakecpr[wxid1, 1] + wakecpr[wxid2, 1]) / 2
                #outside of wake
                vtheta[ip] = 0.0
                continue
            else
                wrid1 = wrid2
            end

        elseif isnothing(wrid2)
            if probe[2] > (wakecpr[end, wxid1] + wakecpr[end, wxid2]) / 2
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
                    Gamma_tilde_rotor[wrid2, xrid] ./ (2 * pi .* wakecpr[wrid2, wxid1])
            else
                vthetaself1 = Gamma_tilde_rotor[wrid1, xrid]# ./ (2 * pi * wakecpr[wrid1, wxid1])
                vthetaself2 = Gamma_tilde_rotor[wrid2, xrid]# ./ (2 * pi * wakecpr[wrid2, wxid1])
                vtheta[ip] =
                    fm.linear(
                        [wakecpr[wrid1, wxid1]; wakecpr[wrid2, wxid1]],
                        [vthetaself1; vthetaself2],
                        probe[2],
                    ) ./ (2 * pi * probe[2])
            end

        elseif wxid1 == wxid2
            # inline with wake control points need to interpolate in r only

            #on edges, just use the edge value
            if wrid1 == wrid2
                vtheta[ip] =
                    Gamma_tilde_grid[wrid1, wxid1] ./ (2 * pi * wakecpr[wrid1, wxid1])
            else
                vtheta1 = Gamma_tilde_grid[wrid1, wxid1]# ./ (2 * pi * wakecpr[wrid1, wxid1])
                vtheta2 = Gamma_tilde_grid[wrid2, wxid1]#./ (2 * pi * wakecpr[wrid2, wxid1])
                vtheta[ip] =
                    fm.linear(
                        [wakecpr[wrid1, wxid1]; wakecpr[wrid2, wxid1]],
                        [vtheta1; vtheta2],
                        probe[2],
                    ) ./ (2 * pi * probe[2])
            end

        elseif wrid1 == wrid2
            #wake radius doesn't change, just use first radial point
            vtheta[ip] = Gamma_tilde_grid[wrid1, wxid1] ./ (2 * pi * wakecpr[wrid1, wxid1])

        else

            # use rotor self induced value if it's closer than the nearest wake control point
            if probe[1] < rotorzloc[xrid]
                vtx2r1 = Gamma_tilde_rotor[wrid1, xrid]# ./ (2 * pi * wakecpr[wrid1])
                vtx2x2 = Gamma_tilde_rotor[wrid2, xrid]# ./ (2 * pi * wakecpr[wrid2])
            else
                vtx2r1 = Gamma_tilde_grid[wrid1, wxid2]# ./ (2 * pi * wakecpr[wrid1])
                vtx2x2 = Gamma_tilde_grid[wrid2, wxid2]# ./ (2 * pi * wakecpr[wrid2])
            end

            if probe[1] > rotorzloc[xrid]
                vtx1r1 = Gamma_tilde_rotor[wrid1, xrid]# ./ (2 * pi * wakecpr[wrid1])
                vtx1r2 = Gamma_tilde_rotor[wrid2, xrid]# ./ (2 * pi * wakecpr[wrid2])
            else
                vtx1r1 = Gamma_tilde_grid[wrid1, wxid1]# ./ (2 * pi * wakecpr[wrid1])
                vtx1r2 = Gamma_tilde_grid[wrid2, wxid1]# ./ (2 * pi * wakecpr[wrid2])
            end

            r1 = fm.linear(
                [wakecpx[wrid1, wxid1]; wakecpx[wrid1, wxid2]],
                [wakecpr[wrid1, wxid1]; wakecpr[wrid1, wxid2]],
                probe[1],
            )
            r2 = fm.linear(
                [wakecpx[wrid2, wxid1]; wakecpx[wrid2, wxid2]],
                [wakecpr[wrid2, wxid1]; wakecpr[wrid2, wxid2]],
                probe[1],
            )

            vt1 = fm.linear(
                [wakecpx[wrid1, wxid1]; wakecpx[wrid1, wxid2]], [vtx1r1; vtx2r1], probe[1]
            )
            vt2 = fm.linear(
                [wakecpx[wrid2, wxid1]; wakecpx[wrid2, wxid2]], [vtx1r2; vtx2x2], probe[1]
            )

            vtheta[ip] = fm.linear([r1; r2], [vt1; vt2], probe[2]) ./ (2 * pi * probe[2])
        end
    end

    if debug
        return Vxr[:, 1],
        Vxr[:, 2], vtheta, Vb[:, 1], Vb[:, 2], Vr[:, 1], Vr[:, 2], Vw[:, 1],
        Vw[:, 2]
    else
        return Vxr[:, 1], Vxr[:, 2], vtheta
    end
end
