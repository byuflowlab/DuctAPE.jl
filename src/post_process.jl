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
    rotor_total_thrust = rotor_inviscid_thrust .+ rotor_viscous_thrust

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

    ## -- Pressure on Bodies -- ##
    duct_inner_cp, duct_outer_cp, hub_cp, duct_inner_vs, duct_outer_vs, hub_vs, duct_inner_x, duct_outer_x, hub_x = get_body_cps(
        iv.gamb,
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
        inputs.body_doublet_panels.endpointidxs,
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
        inputs.vx_dwbte,
        inputs.vr_dwbte,
        inputs.body_doublet_panels.endpointidxs,
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
    total_thrust = sum([rotor_inviscid_thrust'; rotor_viscous_thrust'; body_thrust])

    # - Total Torque - #
    total_torque = sum([rotor_inviscid_torque; rotor_viscous_torque])

    # - Total Power - #
    total_power = sum([rotor_inviscid_power; rotor_viscous_power])

    # - Total Efficiency - #
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
        duct_inner_vs,
        duct_inner_cp,
        duct_inner_x,
        duct_outer_vs,
        duct_outer_cp,
        duct_outer_x,
        hub_vs,
        hub_cp,
        hub_x,
        # - Body Wake Values - #
        # surface velocities and pressures
        ductwake_vs,
        ductwake_cp,
        ductwake_x=inputs.duct_wake_panels.controlpoint[:, 1],
        hubwake_vs,
        hubwake_cp,
        hubwake_x=inputs.hub_wake_panels.controlpoint[:, 1],
        # - Rotor Values - #
        # rotor thrust
        rotor_total_thrust,
        rotor_inviscid_thrust,
        rotor_inviscid_thrust_dist,
        rotor_viscous_thrust,
        rotor_viscous_thrust_dist,
        CT,
        # rotor torque
        rotor_inviscid_torque,
        rotor_inviscid_torque_dist,
        rotor_viscous_torque,
        rotor_viscous_torque_dist,
        CQ,
        # rotor power
        rotor_inviscid_power,
        rotor_inviscid_power_dist,
        rotor_viscous_power,
        rotor_viscous_power_dist,
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
        total_efficiency,
        induced_efficiency,
        ideal_efficiency,
    )

    return out
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
        cp[hidr[ir]] .+= delta_cp(Htilde[1, ir], Stilde[1, ir], v_theta_hub, Vref)
        cp[didr[ir]] .+= delta_cp(Htilde[end, ir], Stilde[end, ir], v_theta_duct, Vref)
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
    gamb, Gamr, sigr, Vm_rotor, Vinf, Vref, B, Omega, dwi, hwi, body_panels, isduct
)

    # - Split body strengths into inner/outer duct and hub - #
    gamdi, gamdo, gamh, xdi, xdo, xh = split_bodies(gamb, body_panels; duct=isduct)

    # - Calculate standard pressure coefficient expression - #
    cpductinner = steady_cp(gamdi, Vinf, Vref)
    cpductouter = steady_cp(gamdo, Vinf, Vref)
    cphub = steady_cp(gamh, Vinf, Vref)

    # - add the change in Cp on the walls due to enthalpy, entropy, and vtheta - #
    calculate_body_delta_cp!(
        cpductinner, cphub, Gamr, sigr, Vm_rotor, Vref, Omega, B, body_panels, dwi, hwi
    )

    return cpductinner, cpductouter, cphub, gamdi, gamdo, gamh, xdi, xdo, xh
end

"""
Calculate the induced velocities on one of the body wakes (unit velocity inputs determine which one)
"""
function calculate_induced_velocities_on_bodywake(
    vx_w, vr_w, gamw, vx_r, vr_r, sigr, vx_b, vr_b, gamb
)

    # problem dimensions
    _, nrotor = size(sigr) # number of rotors
    nwake, _ = size(gamw) # number of wake sheets
    np, _ = size(vx_b[1])

    # initialize outputs
    vx = zeros(eltype(gamw), np) # axial induced velocity
    vr = zeros(eltype(gamw), np) # radial induced velocity

    # add body induced velocities
    @views vx[:] .+= vx_b[1] * gamb
    @views vr[:] .+= vr_b[1] * gamb

    # add wake induced velocities
    for jwake in 1:nwake
        @views vx[:] .+= vx_w[jwake] * gamw[jwake, :]
        @views vr[:] .+= vr_w[jwake] * gamw[jwake, :]
    end

    # add rotor induced velocities
    for jrotor in 1:nrotor
        @views vx[:] .+= vx_r[jrotor] * sigr[:, jrotor]
        @views vr[:] .+= vr_r[jrotor] * sigr[:, jrotor]
    end

    # return raw induced velocities
    return vx, vr
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
        vx_w, vr_w, gamw, vx_r, vr_r, sigr, vx_b, vr_b, gamb
    )

    # get "surface" velocities
    vs =
        (vx_bodywake .+ Vinf) .* cos.(panels.panel_angle) .+
        vr_bodywake .* sin.(panels.panel_angle)

    # - Get steady pressure coefficients - #
    cp_steady = steady_cp(vs, Vinf, Vref)

    # - Get delta cp - #
    deltacp = calculate_bodywake_delta_cp(
        Gamr, sigr, Vm_rotor, Vref, Omega, B, panels.panel_center[:, 2]; body=body
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
    if Vinf == 0.0 || total_power <= 0.0 || total_thrust <= 0.0
        return 0.0
    else
        return total_thrust * Vinf / total_power
    end
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
            clin, cdin, _ = dfdc_clcdcm(
                Wmag_rotor[ir, irotor],
                reynolds,
                solidity[ir, irotor],
                stagger[ir, irotor],
                alpha[ir, irotor],
                afparamsin[ir, irotor],
                asound,
            )
            # get outer values
            clout, cdout, _ = dfdc_clcdcm(
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
#                      get_intermediate_values Intermediate Values                      #
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
        inputs.body_doublet_panels.endpointidxs;
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

#
#
#
#
#
#
#
#
#
#
#
"""
NEED TO TEST. NO GUARENTEES THIS OR CONSTITUENT FUNCTIONS ARE CORRECT. (in fact, they are wrong...)
"""
function probe_velocity_field(
    field_points,
    Vinf;
    body_strengths=nothing,
    body_doublet_panels=nothing,
    wake_strengths=nothing,
    wake_panels=nothing,
    source_strengths=nothing,
    source_panels=nothing,
)

    # Initialize velocities
    Vfield_x = similar(field_points) .= Vinf
    body_induced_x = similar(field_points) .= 0.0
    wake_induced_x = similar(field_points) .= 0.0
    source_induced_x = similar(field_points) .= 0.0
    Vfield_r = similar(field_points) .= 0.0
    body_induced_r = similar(field_points) .= 0.0
    wake_induced_r = similar(field_points) .= 0.0
    source_induced_r = similar(field_points) .= 0.0

    # - add body induced velocity to total
    if !isnothing(body_strengths) && !isnothing(body_doublet_panels)
        # set up influence coefficients based on field points and body panels
        mesh_fpb = generate_field_mesh(body_doublet_panels, field_points)

        A_fpb = assemble_induced_velocity_matrices_infield(
            mesh_fpb, body_doublet_panels, field_points
        )

        # axial components
        vx_s = A_fpb[1]

        # radial components
        vr_s = A_fpb[2]

        # Mutliply things out to get induced velocities
        body_induced_x = vx_s * body_strengths
        body_induced_r = vx_s * body_strengths

        # Add to total velocity field
        Vfield_x .+= body_induced_x
        Vfield_r .+= body_induced_r
    end

    # - add wake induced velocity to total
    if !isnothing(wake_strengths) && !isnothing(wake_panels)
        # set up influence coefficients based on field points and body panels
        mesh_fpw = [
            generate_field_mesh([wake_panels[j]], field_points) for i in 1:1,
            j in 1:length(wake_panels)
        ]

        A_fpw = [
            assemble_induced_velocity_matrices_infield(
                mesh_fpw[i, j], [wake_panels[j]], field_points
            ) for i in 1:1, j in 1:length(wake_panels)
        ]

        # axial components
        vx_s = [A_fpw[i, j][1] for i in 1:1, j in 1:length(wake_panels)]

        # radial components
        vr_s = [A_fpw[i, j][2] for i in 1:1, j in 1:length(wake_panels)]

        # Mutliply things out to get induced velocities
        wake_induced_x = vx_s * wake_strengths
        wake_induced_r = vx_s * wake_strengths

        # Add to total velocity field
        Vfield_x .+= wake_induced_x
        Vfield_r .+= wake_induced_r
    end

    # - add source induced velocity to total
    if !isnothing(source_strengths) && !isnothing(source_panels)
        # set up influence coefficients based on field points and body panels
        mesh_fps = [
            generate_field_mesh([source_panels[j]], field_points) for i in 1:1,
            j in 1:length(source_panels)
        ]

        A_fps = [
            assemble_induced_velocity_matrices_infield(
                mesh_fps[i, j], [source_panels[j]], field_points
            ) for i in 1:1, j in 1:length(source_panels)
        ]

        # axial components
        vx_s = [A_fps[i, j][1] for i in 1:1, j in 1:length(source_panels)]

        # radial components
        vr_s = [A_fps[i, j][2] for i in 1:1, j in 1:length(source_panels)]

        # Mutliply things out to get induced velocities
        source_induced_x = vx_s * source_strengths
        source_induced_r = vx_s * source_strengths

        # Add to total velocity field
        Vfield_x .+= source_induced_x
        Vfield_r .+= source_induced_r
    end

    return Vfield_x,
    Vfield_r,
    body_induced_x,
    body_induced_r,
    wake_induced_x,
    wake_induced_r,
    source_induced_x,
    source_induced_r
end
