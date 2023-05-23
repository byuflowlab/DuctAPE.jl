#=

Various Post-processing functions

=#

function post_process(states, inputs)

    # extract states
    gamb, gamw, Gamr, sigr = extract_state_variables(states, inputs)

    # get problem dimensions
    nr, nrotor = size(Gamr)
    nw = nr + 1

    # - Extract convenient input fields - #
    Vinf = inputs.Vinf
    Vref = inputs.reference_parameters.Vref
    Rref = inputs.reference_parameters.Rref
    rho = inputs.freestream.rho
    mu = inputs.freestream.mu
    asound = inputs.freestream.asound
    dwi = inputs.ductwakeidx
    hwi = inputs.hubwakeidx
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
    rpl = reshape(
        reduce(vcat, (p -> p.panel_length).(inputs.rotor_source_panels)), (nr, nrotor)
    )

    ## -- Rotor Outputs -- ##
    # - Set Up - #
    # get induced velocities on the rotor
    vx_rotor, vr_rotor, vtheta_rotor = calculate_induced_velocities_on_rotors(
        inputs.blade_elements,
        Gamr,
        inputs.vx_rw,
        inputs.vr_rw,
        gamw,
        inputs.vx_rr,
        inputs.vr_rr,
        sigr,
        inputs.vx_rb,
        inputs.vr_rb,
        gamb,
    )

    # reframe velocities
    Wx_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor = reframe_rotor_velocities(
        vx_rotor, vr_rotor, vtheta_rotor, Vinf, Omega, rpc
    )

    # get net circulation
    Gamma_tilde = calculate_net_circulation(Gamr, B)

    # get blade element coefficients
    clift, cdrag, inflow_angle, angle_of_attack = get_blade_aero(
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
        rho,
        mu,
        asound,
    )

    # - Rotor Thrust - #
    # inviscid thrust
    rotor_inviscid_thrust, rotor_inviscid_thrust_dist = inviscid_rotor_trust(
        Wtheta_rotor, Gamma_tilde, rpl, rho
    )

    # viscous thrust
    rotor_viscous_thrust, rotor_viscous_thrust_dist = viscous_rotor_thrust(
        Wx_rotor, Wmag_rotor, B, chord, rpl, cdrag, rho
    )

    # total thrust
    rotor_total_thrust = rotor_inviscid_thrust .+ rotor_viscous_thrust

    # - Rotor Torque - #
    # inviscid torque
    rotor_inviscid_torque, rotor_inviscid_torque_dist = inviscid_rotor_torque(
        Wx_rotor, rpc, rpl, Gamma_tilde, rho
    )

    # viscous torque
    rotor_viscous_torque, rotor_viscous_torque_dist = viscous_rotor_torque(
        Wtheta_rotor, Wmag_rotor, B, chord, rpc, rpl, cdrag, rho
    )

    # - Rotor Power - #
    # inviscid power
    rotor_inviscid_power = inviscid_rotor_power(rotor_inviscid_torque, Omega)
    rotor_inviscid_power_dist = inviscid_rotor_power(rotor_inviscid_torque_dist, Omega)

    # viscous power
    rotor_viscous_power = viscous_rotor_power(rotor_viscous_torque, Omega)
    rotor_viscous_power_dist = viscous_rotor_power(rotor_viscous_torque_dist, Omega)

    ## -- Pressure on Bodies -- ##
    duct_inner_cp, duct_outer_cp, hub_cp, duct_inner_c, duct_outer_x, hub_x = get_cps(
        gamb,
        gamw,
        Gamr,
        sigr,
        Wm_rotor,
        Vinf,
        Vref,
        B,
        Omega,
        dwi,
        hwi,
        inputs.body_panels,
        inputs.isduct,
    )

    ## -- Duct Outputs -- ##
    # - Put duct pressures together - #
    duct_cp = [duct_inner_cp; duct_outer_cp]

    # - Calculate Thrust from Bodies - #

    duct_thrust, _ = forces_from_pressure(
        duct_cp, inputs.body_panels[1]; rho=rho, Vref=Vref
    )

    hub_thrust, _ = forces_from_pressure(hub_cp, inputs.body_panels[2]; rho=rho, Vref=Vref)

    ## -- Total Outputs -- ##

    # - Total Thrust - #
    total_thrust = sum(
        [rotor_inviscid_thrust; rotor_viscous_thrust; duct_thrust; hub_thrust]
    )

    # - Total Torque - #
    total_torque = sum([rotor_inviscid_torque; rotor_viscous_torque])

    # - Total Power - #
    total_power = sum([rotor_inviscid_power; rotor_viscous_power])

    # - Total Efficiency - #
    total_efficiency = get_total_efficiency(total_thrust, total_power, Vinf)

    # - Induced Efficiency - #
    induced_efficiency = get_induced_efficiency(
        rotor_inviscid_thrust, duct_thrust + hub_thrust, rotor_inviscid_power, Vinf
    )

    # - Ideal Efficiency - #
    ideal_efficiency = get_ideal_efficiency(total_thrust, rho, Vinf, Rref)

    # - Blade Loading - #
    blade_normal_force_per_unit_span, blade_tangential_force_per_unit_span = get_blade_loads(
        Wmag_rotor, inflow_angle, clift, cdrag, chord, rho
    )

    # - Thrust and Torque Coefficients - #
    CT, CQ = tqcoeff(total_thrust, total_torque, rho, Omega[1], Rref)

    ## -- Assemble Output Tuple -- ##

    out = (;
        # - Body Values - #
        # body thrust
        duct_thrust,
        hub_thrust,
        body_thrust=duct_thrust + hub_thrust,
        duct_inner_cp,
        duct_inner_c,
        duct_outer_cp,
        duct_outer_x,
        hub_cp,
        hub_x,
        # - Rotor Values - #
        rotor_panel_center=rpc,
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

function steady_cp(vs, vinf, vref)
    return (vinf^2 .- vs .^ 2) / vref^2
end

function delta_cp(deltaH, deltaS, Vtheta, Vref)
    return (2.0 * (deltaH - deltaS) .- Vtheta .^ 2) / Vref^2
end

"""
move to utils.jl
"""
function split_bodies(vec, panels; duct=true)
    # get type of vector for consistent outputs
    TF = eltype(vec)

    #check if duct is used
    if !duct
        #hub only
        return TF[], TF[], vec, TF[], TF[], panels.panel_center[:, 1]
    else
        # get duct leading edge index. assumes duct comes first in vector
        _, leidx = findmin(panels[1].panel_center[:, 1])
        ndpan = length(panels[1].panel_center[:, 1])

        if length(panels) > 1
            #duct and hub
            return vec[1:leidx],
            vec[(leidx + 1):ndpan],
            vec[(ndpan + 1):end],
            panels[1].panel_center[1:leidx, 1],
            panels[1].panel_center[(leidx + 1):ndpan, 1],
            panels[2].panel_center[:, 1]
        else
            #duct only
            return vec[1:leidx],
            vec[(leidx + 1):ndpan],
            TF[],
            panels[1].panel_center[1:leidx, 1],
            panels[1].panel_center[(leidx + 1):ndpan, 1],
            TF[]
        end
    end

    # shouldn't get to this point...
    return nothing
end

"""
"""
function calculate_delta_cp(
    gamb, gamw, Gamr, sigr, Vm_rotor, Vinf, Omega, B, body_panels, dwi, hwi, Vref
)

    ## -- Calculate change in pressure coefficient -- ##
    # - Get the tangential velocities on the bodies - #
    # TODO; changes needed for multiple rotors
    v_theta_duct, v_theta_hub = vtheta_on_body(
        B[1] * Gamr,
        body_panels[1].panel_center[dwi, 2],
        body_panels[2].panel_center[hwi, 2],
    )

    # - Calculate enthalpy disk jump - #
    Htilde = calculate_enthalpy_jumps(Gamr, Omega, B)

    # - Calculate entropy disk jump - #
    Stilde = calculate_entropy_jumps(sigr, Vm_rotor)

    # assemble change in cp due to enthalpy and entropy behind rotor(s)
    deltacphub = delta_cp(Htilde[1], Stilde[1], v_theta_hub, Vref)
    deltacpduct = delta_cp(Htilde[end], Stilde[end], v_theta_duct, Vref)

    return deltacpduct, deltacphub
end

"""
"""
function vtheta_on_body(BGamr, inner_duct_r, hub_r)
    v_theta_duct = BGamr[end] ./ (2.0 * pi * inner_duct_r)
    v_theta_hub = BGamr[1] ./ (2.0 * pi * hub_r)

    return v_theta_duct, v_theta_hub
end

"""
calculate pressure coefficient distribution on duct/hub walls
formulation taken from DFDC source code. TODO: derive where the expressions came from.
"""
function get_cps(
    gamb, gamw, Gamr, sigr, Vm_rotor, Vinf, Vref, B, Omega, dwi, hwi, body_panels, isduct
)

    # - Split body strengths into inner/outer duct and hub - #
    gamdi, gamdo, gamh, xdi, xdo, xh = split_bodies(
        gamb, body_panels; duct=isduct
    )

    # - Calculate standard pressure coefficient expression - #
    cpductinner = steady_cp(gamdi, Vinf, Vref)
    cpductouter = steady_cp(gamdo, Vinf, Vref)
    cphub = steady_cp(gamh, Vinf, Vref)

    # - Calculate the change in Cp on the walls due to enthalpy, entropy, and vtheta - #
    deltacpduct, deltacphub = calculate_delta_cp(
        gamb, gamw, Gamr, sigr, Vm_rotor, Vinf, Omega, B, body_panels, dwi, hwi, Vref
    )

    # - add raw and adjusted cp values together - #
    cpductinner[dwi] .+= deltacpduct
    cphub[hwi] .+= deltacphub

    return cpductinner, cpductouter, cphub, xdi, xdo, xh
end

"""

Calculate dimensional and non-dimensional axial force on a single body
"""
function forces_from_pressure(cps, panels; rho=1.225, Vref=1.0)

    # - rename for convenience - #
    #just want x-component of normals since it's axisymmetric
    ns = panels.panel_normal[:, 1]
    #radial positions
    rs = panels.panel_center[:, 2]
    #panel lengths
    ds = panels.panel_length
    # dimensions
    np = length(cps)

    # - initialize - #
    cfx = 0.0 # axial force coefficient (all others are zero for axisymmetric case)
    # - rectangular integration due to constant panel strengths. - #
    for i in 1:np
        cfx += cps[i] * ns[i] * ds[i] * 2.0 * pi * rs[i]
    end

    #dimensionalize
    q = 0.5 * rho * Vref^2

    #note, thrust is in negative x-direction
    return cfx * q, cfx
end

function inviscid_rotor_trust(Wtheta, Gamma_tilde, rotor_panel_length, rho)

    # problem dimensions
    nr, nrotor = size(Gamma_tilde)

    # initialize
    dTi = similar(Gamma_tilde) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            # section thrust
            dTi[ir, irotor] =
                -rho *
                Gamma_tilde[ir, irotor] *
                Wtheta[ir, irotor] *
                rotor_panel_length[ir, irotor]
        end
    end

    #sum the section thrust
    Tinv = sum(dTi; dims=1)

    return Tinv, dTi
end

function viscous_rotor_thrust(Wx_rotor, Wmag_rotor, B, chord, rotor_panel_length, cd, rho)

    # get dimensions
    nr, nrotor = size(Wx_rotor)

    #initialize
    dTv = similar(Wx_rotor) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            hrwc = 0.5 * rho * Wmag_rotor[ir, irotor] * chord[ir, irotor]
            bdr = B[irotor] * rotor_panel_length[ir, irotor]
            dTv[ir, irotor] = -hrwc * cd[ir, irotor] * Wx_rotor[ir, irotor] * bdr
        end
    end

    Tvisc = sum(dTv; dims=1)

    return Tvisc, dTv
end

function inviscid_rotor_torque(
    Wx_rotor, rotor_panel_center, rotor_panel_length, Gamma_tilde, rho
)

    # dimensions
    nr, nrotor = size(Gamma_tilde)

    # initialize
    dQi = similar(Gamma_tilde) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            rdr = rotor_panel_center[ir, irotor] * rotor_panel_length[ir, irotor]
            dQi[ir, irotor] = rho * Gamma_tilde[ir, irotor] * Wx_rotor[ir, irotor] * rdr
        end
    end

    Qinv = sum(dQi; dims=1)

    return Qinv, dQi
end

function viscous_rotor_torque(
    Wtheta_rotor, Wmag_rotor, B, chord, rotor_panel_center, rotor_panel_length, cd, rho
)

    # dimensions
    nr, nrotor = size(Wtheta_rotor)

    # initialize
    dQv = similar(chord) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            hrwc = 0.5 * rho * Wmag_rotor[ir, irotor] * chord[ir, irotor]
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
    if Vinf != 0.0
        return total_thrust * Vinf / total_power
    else
        return 0.0
    end
end

function get_induced_efficiency(Tinv, Tduct, Pinv, Vinf)
    if Vinf != 0.0
        return Vinf * (Tinv .+ Tduct) ./ Pinv
    else
        return 0.0
    end
end

function get_ideal_efficiency(total_thrust, rho, Vinf, Rref)
    if Vinf != 0.0
        TC = total_thrust / (0.5 * rho * Vinf^2 * pi * Rref^2)
        return 2.0 / (1.0 + sqrt(max(TC, -1.0) + 1.0))
    else
        return 0.0
    end
end

function tqcoeff(thrust, torque, rho, Omega, Rref)

    # reference diameter
    D = 2.0*Rref

    # rototion in rev per second
    n = Omega /(2.0*pi)

    # thrust coefficient
    CT = thrust/(rho*n^2*D^4)

    # torque coefficient
    CQ = torque/(rho*n^2*D^5)

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
    rho,
    mu,
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
            reynolds = chord[ir, irotor] * abs(Wmag_rotor[ir, irotor]) * rho / mu

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

function get_blade_loads(Wmag_rotor, phi, cl, cd, chords, rho)

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
            Np[ir, irotor] = cn * 0.5 * rho * Wmag_rotor[ir, irotor]^2 * chords[ir, irotor]
            Tp[ir, irotor] = ct * 0.5 * rho * Wmag_rotor[ir, irotor]^2 * chords[ir, irotor]
        end
    end

    Npfull = [zeros(nrotor)'; Np; zeros(nrotor)']
    Tpfull = [zeros(nrotor)'; Tp; zeros(nrotor)']

    #TODO: consider comparing this with DFDC versions for them

    # ## -- Integrate Loads to get Thrust and Torque
    # # add hub/tip for complete integration.  loads go to zero at hub/tip.
    # rfull = [blade_elements.Rhub; blade_elements.rbe; blade_elements.Rtip]

    # # thrust and torqe distributions
    # thrust = Npfull
    # torque = Tpfull .* rfull

    # # integrate Thrust and Torque (trapezoidal)
    # T = blade_elements.B * fm.trapz(rfull, thrust)
    # Q = blade_elements.B * fm.trapz(rfull, torque)
    # P = Q * blade_elements.Omega

    return Npfull, Tpfull
end

"""
"""
function dump(states, inputs)

    # - Extract commonly used items from precomputed inputs - #
    blade_elements = inputs.blade_elements
    rpc = inputs.rotor_panel_centers
    Vinf = inputs.Vinf

    # - Extract states - #
    gamb, gamw, Gamr, sigr = extract_state_variables(states, inputs)

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
        gamb;
        debug=true,
    )

    vx_rotor, vr_rotor, vtheta_rotor, Wx_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor = calculate_rotor_velocities(
        Gamr, gamw, sigr, gamb, inputs
    )

    Gamma_tilde = calculate_net_circulation(Gamr, blade_elements.B)

    H_tilde = calculate_enthalpy_jumps(Gamr, blade_elements.Omega, blade_elements.B)

    _, _, phi, alpha, cl, cd = calculate_gamma_sigma(
        blade_elements, Wm_rotor, Wtheta_rotor, Wmag_rotor, inputs.freestream; debug=true
    )

    return (;
        gamb,
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
    body_panels=nothing,
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
    if !isnothing(body_strengths) && !isnothing(body_panels)
        # set up influence coefficients based on field points and body panels
        mesh_fpb = generate_field_mesh(body_panels, field_points)

        A_fpb = assemble_induced_velocity_matrices_infield(
            mesh_fpb, body_panels, field_points
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
