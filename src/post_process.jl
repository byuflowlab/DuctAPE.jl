#=

Various Post-processing functions

=#

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
function calculate_delta_cp(gamb, gamw, Gamr, sigr, inputs, Vref)
    # - Extract convenient input fields - #
    Vinf = inputs.Vinf
    Omega = inputs.blade_elements.Omega
    B = inputs.blade_elements.B
    dwi = inputs.ductwakeidx[1]
    hwi = inputs.hubwakeidx[1]

    # Fill out wake strengths
    wake_vortex_strengths = fill_out_wake_strengths(
        gamw,
        inputs.rotor_indices_in_wake,
        inputs.num_wake_x_panels;
        ductTE_index=inputs.ductTE_index,
        hubTE_index=inputs.hubTE_index,
        interface="hard",
    )

    ## -- Calculate change in pressure coefficient -- ##
    # - Get the meridional and tangential velocities at the rotor - #
    _, _, _, _, _, Vm, _ = calculate_rotor_velocities(
        Gamr, wake_vortex_strengths, sigr, gamb, inputs
    )

    v_theta_duct, v_theta_hub = vtheta_on_body(
        B[1] * Gamr,
        inputs.body_panels[1].panel_center[dwi, 2],
        inputs.body_panels[2].panel_center[hwi, 2],
    )

    # - Calculate enthalpy disk jump - #
    Htilde = calculate_enthalpy_jumps(Gamr, Omega, B)

    # - Calculate entropy disk jump - #
    Stilde = calculate_entropy_jumps(sigr, Vm)

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
function get_cps(states, inputs, Vref)

    # - Extract convenient input fields - #
    Vinf = inputs.Vinf
    dwi = inputs.ductwakeidx[1]
    hwi = inputs.hubwakeidx[1]

    # - Extract States - #
    gamb, gamw, Gamr, sigr = extract_state_variables(states, inputs)

    # - Split body strengths into inner/outer duct and hub - #
    gamdi, gamdo, gamh, xdi, xdo, xh = split_bodies(
        gamb, inputs.body_panels; duct=inputs.isduct
    )

    # - Calculate standard pressure coefficient expression - #
    cpductinner = steady_cp(gamdi, Vinf, Vref)
    cpductouter = steady_cp(gamdo, Vinf, Vref)
    cphub = steady_cp(gamh, Vinf, Vref)

    # - Calculate the change in Cp on the walls due to enthalpy, entropy, and vtheta - #
    deltacpduct, deltacphub = calculate_delta_cp(gamb, gamw, Gamr, sigr, inputs, Vref)

    # - add raw and adjusted cp values together - #
    cpductinner[dwi] .+= deltacpduct
    cphub[hwi] .+= deltacphub

    return cpductinner, cpductouter, cphub, xdi, xdo, xh
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

    # - Fill out wake strengths - #
    wake_vortex_strengths = fill_out_wake_strengths(
        gamw,
        inputs.rotor_indices_in_wake,
        inputs.num_wake_x_panels;
        ductTE_index=inputs.ductTE_index,
        hubTE_index=inputs.hubTE_index,
        interface="hard",
    )

    _, _, _, vxfrombody, vrfrombody, vxfromwake, vrfromwake, vxfromrotor, vrfromrotor = calculate_induced_velocities_on_rotors(
        blade_elements,
        Gamr,
        inputs.vx_rw,
        inputs.vr_rw,
        wake_vortex_strengths,
        inputs.vx_rr,
        inputs.vr_rr,
        sigr,
        inputs.vx_rb,
        inputs.vr_rb,
        gamb;
        debug=true,
    )

    vx_rotor, vr_rotor, vtheta_rotor, Wx_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor = calculate_rotor_velocities(
        Gamr, wake_vortex_strengths, sigr, gamb, inputs
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
