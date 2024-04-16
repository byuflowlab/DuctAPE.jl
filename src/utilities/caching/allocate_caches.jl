"""
"""
function allocate_body_panel_container(total_length, problem_dimensions)
    (;
        ndn,    # number of duct nodes
        ncbn,   # number of centerbody nodes
    ) = problem_dimensions

    nn = [ndn, ncbn]
    # number of panels to generate for each body
    np = nn .- 1
    # number of bodies
    nb = length(np)
    # total number of panels in system
    tp = sum(np)
    # total number of nodes in system
    tn = tp + nb

    return allocate_panel_container(total_length, nn, np, tn, tp, nb)
end

function allocate_rotor_panel_container(total_length, problem_dimensions)
    (;
        nrotor,     # number of rotors
        nws,    # number of wake sheets (also rotor nodes)
    ) = problem_dimensions

    nn = nws * ones(Int, nrotor)
    # number of panels to generate for each body
    np = nn .- 1
    # number of bodies
    nb = length(np)
    # total number of panels in system
    tp = sum(np)
    # total number of nodes in system
    tn = tp + nb

    return allocate_panel_container(total_length, nn, np, tn, tp, nb)
end

function allocate_wake_panel_container(total_length, problem_dimensions)
    (;
        nws,    # number of wake sheets (also rotor nodes)
        nwsn,   # number of nodes in each wake sheet
    ) = problem_dimensions

    nn = nwsn * ones(Int, nws)
    # number of panels to generate for each body
    np = nn .- 1
    # number of bodies
    nb = length(np)
    # total number of panels in system
    tp = sum(np)
    # total number of nodes in system
    tn = tp + nb

    return allocate_panel_container(total_length, nn, np, tn, tp, nb)
end

"""
"""
function allocate_panel_container(total_length, nn, np, tn, tp, nb)
    s = size(nn)
    l = lfs(s)
    nnode = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # number of panels to generate for each body
    npanel = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (1,)
    l = lfs(s)
    nbodies = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (1,)
    l = lfs(s)
    totpanel = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (1,)
    l = lfs(s)
    totnode = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (2, tn)
    l = lfs(s)
    node = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (2, tp)
    l = lfs(s)
    controlpoint = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    normal = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    tangent = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    nodemap = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (tp)
    l = lfs(s)
    influence_length = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nb, 2, 2)
    l = lfs(s)
    endnodes = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    tenode = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (2, nb)
    l = lfs(s)

    itcontrolpoint = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    itnormal = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    ittangent = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    tenormal = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    tendotn = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    tencrossn = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    endnodeidxs = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    endpanelidxs = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    teadjnodeidxs = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nb,)
    l = lfs(s)
    teinfluence_length = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    prescribednodeidxs = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    return (;
        nnode,
        npanel,
        nbodies,
        totpanel,
        totnode,
        node,
        controlpoint,
        normal,
        tangent,
        nodemap,
        influence_length,
        endnodes,
        tenode,
        itcontrolpoint,
        itnormal,
        ittangent,
        tenormal,
        tendotn,
        tencrossn,
        endnodeidxs,
        endpanelidxs,
        teadjnodeidxs,
        teinfluence_length,
        prescribednodeidxs,
    ),
    total_length
end

"""
"""
function allocate_panel_containers(problem_dimensions, total_length)
    body_vortex_panels, total_length = allocate_body_panel_container(
        total_length, problem_dimensions
    )
    rotor_source_panels, total_length = allocate_rotor_panel_container(
        total_length, problem_dimensions
    )
    wake_vortex_panels, total_length = allocate_wake_panel_container(
        total_length, problem_dimensions
    )

    return (; body_vortex_panels, rotor_source_panels, wake_vortex_panels), total_length
end

"""
"""
function allocate_prepost_container_cache(paneling_constants::PanelingConstants)
    problem_dimensions = get_problem_dimensions(paneling_constants)

    return allocate_prepost_container_cache(problem_dimensions)
end

function allocate_prepost_container_cache(problem_dimensions)

    (;
        nrotor,     # number of rotors
        nwn,    # number of wake nodes
        nwp,    # number of wake panels
        ncp,    # number of casing panels
        ndn,    # number of duct nodes
        ncbn,   # number of centerbody nodes
        nbn,    # number of body nodes
        nbp,    # number of body panels
        nws,    # number of wake sheets (also rotor nodes)
        nbe,    # number of blade elements (also rotor panels)
        nwsn,   # number of nodes in each wake sheet
        nbodies,
    ) = problem_dimensions

    # - initialize - #
    total_length = 0

    ### --- PRE-PROCESSING --- ###

    # - COORDINATES - #
    dcshape = (2, ndn)
    dclength = lfs(dcshape)
    rp_duct_coordinates = (;
        index=(total_length + 1):(total_length + dclength), shape=dcshape
    )
    total_length += dclength

    cbcshape = (2, ncbn)
    cbclength = lfs(cbcshape)
    rp_centerbody_coordinates = (;
        index=(total_length + 1):(total_length + cbclength), shape=cbcshape
    )
    total_length += cbclength

    wgshape = (2, nwsn, nws)
    wglength = lfs(wgshape)
    wake_grid = (; index=(total_length + 1):(total_length + wglength), shape=wgshape)
    total_length += wglength

    rs = (nrotor,)
    rl = lfs(rs)
    rotor_indices_in_wake = (; index=(total_length + 1):(total_length + rl), shape=rs)
    total_length += rl

    # - PANELS - #
    panels, total_length = allocate_panel_containers(problem_dimensions, total_length)

    # - INDUCED VELS - #
    s = (nbp, nbn, 2)
    l = lfs(s)
    v_bb = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbp, nrotor * nws, 2)
    l = lfs(s)
    v_br = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbp, nwn, 2)
    l = lfs(s)
    v_bw = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - LINSYS - #
    s = (nbp, nbn)
    l = lfs(s)
    AICn = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (2, nbn)
    l = lfs(s)
    AICpcp = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbp,)
    l = lfs(s)
    vdnb = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (2,)
    l = lfs(s)
    vdnpcp = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    ##### ----- POST PROCESSING ----- #####

    ### --- ROTOR Post-Processing Cache --- ###

    s = (nrotor,)
    l = lfs(s)
    rotor_inviscid_thrust = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_viscous_thrust = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_thrust = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_inviscid_torque = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_viscous_torque = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_torque = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_inviscid_power = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_viscous_power = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_power = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_CT = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_CQ = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_CP = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_efficiency = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    induced_efficiency= (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l


    s = (nbe, nrotor)
    l = lfs(s)
    rotor_inviscid_thrust_dist = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_viscous_thrust_dist = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_inviscid_torque_dist = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_viscous_torque_dist = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_inviscid_power_dist = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_viscous_power_dist = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    blade_normal_force_per_unit_span = (;
        index=(total_length + 1):(total_length + l), shape=s
    )
    total_length += l

    blade_tangential_force_per_unit_span = (;
        index=(total_length + 1):(total_length + l), shape=s
    )
    total_length += l

    cn = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    ct = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    cphi = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    sphi = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l


    ### --- BODY Post-Processing Cache --- ###
    s = (2, nbp)
    l = lfs(s)
    Vtot_in = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Vtot_out = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Vtot_prejump = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vtot_body = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vtot_jump = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vtot_wake = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vtot_rotors = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbp,)
    l = lfs(s)
    Vtan_in = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Vtan_out = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    cp_in = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l
    cp_out = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (ncp,)
    l = lfs(s)
    vtan_casing_in = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vtan_casing_out = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    casing_zpts = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    cp_casing_in = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l
    cp_casing_out = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (ndn - 1 - ncp,)
    l = lfs(s)
    vtan_nacelle_in = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vtan_nacelle_out = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    nacelle_zpts = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    cp_nacelle_in = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l
    cp_nacelle_out = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (ncbn - 1,)
    l = lfs(s)
    vtan_centerbody_in = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vtan_centerbody_out = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    centerbody_zpts = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    cp_centerbody_in = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    cp_centerbody_out = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (ndn - 1,)
    l = lfs(s)
    duct_jump = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (ncbn - 1,)
    l = lfs(s)
    centerbody_jump = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbp,)
    l = lfs(s)
    body_jump_term = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbodies,)
    l = lfs(s)
    body_thrust = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    body_force_coefficient = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    ### --- TOTALS Post-Processing Cache --- ###
    s = (1,)
    l = lfs(s)
    total_thrust = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    total_torque = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    total_power = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    total_CT = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    total_CQ = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    total_CP = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    total_efficiency = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    ideal_efficiency = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # return tuple of initialized cache and associated dimensions
    return (;
        prepost_container_cache=PreallocationTools.DiffCache(zeros(total_length)),
        prepost_container_cache_dims=(;
            ### --- PRE --- ###
            rp_duct_coordinates,
            rp_centerbody_coordinates,
            wake_grid,
            rotor_indices_in_wake,
            panels,
            ivb=(; v_bb, v_br, v_bw),
            AICn,
            AICpcp,
            vdnb,
            vdnpcp,
            ### --- Post --- ###
            # - ROTOR - #
            rotor_inviscid_thrust,
            rotor_inviscid_thrust_dist,
            rotor_viscous_thrust,
            rotor_viscous_thrust_dist,
            rotor_thrust,
            rotor_inviscid_torque,
            rotor_inviscid_torque_dist,
            rotor_viscous_torque,
            rotor_viscous_torque_dist,
            rotor_torque,
            rotor_inviscid_power,
            rotor_inviscid_power_dist,
            rotor_viscous_power,
            rotor_viscous_power_dist,
            rotor_power,
            rotor_CT,
            rotor_CQ,
            rotor_CP,
            rotor_efficiency,
            induced_efficiency,
            blade_normal_force_per_unit_span,
            blade_tangential_force_per_unit_span,
            blade_loading_intermediate_containers=(; cn, ct, cphi, sphi),
            # - BODY - #
            zpts=(; casing_zpts, nacelle_zpts, centerbody_zpts),
            vtan_tuple=(;
                Vtot_in,
                Vtot_out,
                Vtan_in,
                Vtan_out,
                Vtot_prejump,
                vtot_body,
                duct_jump,
                centerbody_jump,
                body_jump_term,
                vtot_jump,
                vtot_wake,
                vtot_rotors,
                vtan_casing_in,
                vtan_casing_out,
                vtan_nacelle_in,
                vtan_nacelle_out,
                vtan_centerbody_in,
                vtan_centerbody_out,
            ),
            cp_tuple=(;
                cp_in,
                cp_out,
                cp_casing_in,
                cp_casing_out,
                cp_nacelle_in,
                cp_nacelle_out,
                cp_centerbody_in,
                cp_centerbody_out,
            ),
            body_thrust,
            body_force_coefficient,
            # - TOTALS - #
            total_thrust,
            total_torque,
            total_power,
            total_CT,
            total_CQ,
            total_CP,
            total_efficiency,
            ideal_efficiency,
        ),
    )
end

"""
"""
function allocate_solve_parameter_cache(
    solve_type::CSORSolverOptions, paneling_constants::PanelingConstants; fd_chunk_size=12, levels=2
)

    # - Get problem dimensions - #
    problem_dimensions = get_problem_dimensions(paneling_constants)

return allocate_solve_parameter_cache(
    solve_type::CSORSolverOptions, problem_dimensions; fd_chunk_size=fd_chunk_size, levels=levels
)
end

function allocate_solve_parameter_extras(solver_options::SIAMFANLEOptions, input_length, total_length)

    s = (input_length,)
    l = lfs(s)
    resid_cache_vec = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (input_length,)
    l = lfs(s)
    jvp_cache_vec = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (input_length,max(2,solver_options.linear_iteration_limit+1))
    l = lfs(s)
    krylov_cache_vec = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    return total_length, (; resid_cache_vec, krylov_cache_vec, jvp_cache_vec)
end

function allocate_solve_parameter_extras(solver_options, input_length, total_length)

    return total_length, (;)
end

function allocate_solve_parameter_cache(
    solve_type::CSORSolverOptions, problem_dimensions; fd_chunk_size=12, levels=2
)

    (;
        nrotor,     # number of rotors
        nwn,    # number of wake nodes
        nwp,    # number of wake panels
        nbn,    # number of body nodes
        nbp,    # number of body panels
        nws,    # number of wake sheets (also rotor nodes)
        nbe,    # number of blade elements (also rotor panels)
    ) = problem_dimensions

    # - initialize - #
    total_length = 0

    # - Initial Guesses - #
    s = (nbe,nrotor)
    l = lfs(s)
    Gamr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nws,nrotor)
    l = lfs(s)
    sigr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwn,)
    l = lfs(s)
    gamw = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # save state dimensions
    state_dims = (; Gamr, sigr, gamw)

    # - Operating Point - #
    s = (1,)
    l = lfs(s)
    Vinf = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rhoinf = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    muinf = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    asound = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nrotor,)
    l = lfs(s)
    Omega = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - Induced Velocities on Rotors - #
    s = (nrotor * nbe, nbn, 2)
    l = lfs(s)
    v_rb = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nrotor * nbe, nrotor * nws, 2)
    l = lfs(s)
    v_rr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nrotor * nbe, nwn, 2)
    l = lfs(s)
    v_rw = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - Induced Velocities on Wakes - #
    s = (nwp, nbn, 2)
    l = lfs(s)
    v_wb = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwp, nrotor * nws, 2)
    l = lfs(s)
    v_wr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwp, nwn, 2)
    l = lfs(s)
    v_ww = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - Linear System - #

    s = (nbn + 2, nbn + 2)
    l = lfs(s)
    A_bb = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbn + 2,)
    l = lfs(s)
    b_bf = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbp, nwn)
    l = lfs(s)
    A_bw = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (2, nwn)
    l = lfs(s)
    A_pw = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbp, nrotor * nws)
    l = lfs(s)
    A_br = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (2, nrotor * nws)
    l = lfs(s)
    A_pr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - Blade Elements - #
    s = (nrotor,)
    l = lfs(s)
    B = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Rtip = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Rhub = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    fliplift = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbe, nrotor)
    l = lfs(s)
    chords = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    twists = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    stagger = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    solidity = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_panel_centers = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    inner_fraction = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - Wake constants - #
    s = (nwn,)
    l = lfs(s)
    wakeK = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # TODO: WHAT TO DO ABOUT AIRFOILS?
    # NOTE: how to know the size of airfoil objects up front?  if it's a DFDC type, this is pretty simple, but if it's a CCBlade type, then the size could be anything depending on how many data points are used
    # airfoils = allocate_airfoil_cache(aftype, nrotor, nbe)

    return (;
        solve_parameter_cache=PreallocationTools.DiffCache(
            zeros(total_length), fd_chunk_size; levels=levels
        ),
        solve_parameter_cache_dims=(;
            state_dims ,
            Gamr,
            sigr,
            gamw,
            operating_point=(; Vinf, rhoinf, muinf, asound, Omega),
            ivr=(; v_rb, v_rr, v_rw),
            ivw=(; v_wb, v_wr, v_ww),
            linsys=(; A_bb, b_bf, A_bw, A_pw, A_br, A_pr),
            blade_elements=(;
                B,
                fliplift,
                chords,
                twists,
                stagger,
                solidity,
                rotor_panel_centers,
                inner_fraction,
                Rtip,
                Rhub,
            ),
            wakeK,
        ),
    )
end

"""
"""
function allocate_solve_parameter_cache(
    solve_type::TS, paneling_constants::PanelingConstants; fd_chunk_size=12, levels=2
) where {TS<:Union{ExternalSolverOptions, MultiSolverOptions}}

    # - Get problem dimensions - #
    problem_dimensions = get_problem_dimensions(paneling_constants)

    return allocate_solve_parameter_cache(
        solve_type, problem_dimensions; fd_chunk_size=fd_chunk_size, levels=levels
    )
end

"""
"""
function allocate_solve_parameter_cache(
    solve_type::TS, problem_dimensions; fd_chunk_size=12, levels=2
) where {TS<:Union{ExternalSolverOptions, MultiSolverOptions}}
    (;
        nrotor,     # number of rotors
        nwn,    # number of wake nodes
        nwp,    # number of wake panels
        nbn,    # number of body nodes
        nbp,    # number of body paneallocate_solve_parameter_extrasheets (also rotor nodes)
        nbe,    # number of blade elements (also rotor panels)
        nws, #number of wake sheets
    ) = problem_dimensions

    # - initialize - #
    total_length = 0

    # - Initial Guesses - #
    s = (nbe, nrotor)
    l = lfs(s)
    vz_rotor = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vtheta_rotor = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwp,)
    l = lfs(s)
    Cm_wake = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    total_length, SIAMFANLE_cache_vecs =  allocate_solve_parameter_extras(solve_type, total_length, total_length)

    # save state dimensions
    state_dims = (; vz_rotor, vtheta_rotor, Cm_wake)

    # - Operating Point - #
    s = (1,)
    l = lfs(s)
    Vinf = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rhoinf = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    muinf = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    asound = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nrotor,)
    l = lfs(s)
    Omega = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - Induced Velocities on Rotors - #
    s = (nrotor * nbe, nbn, 2)
    l = lfs(s)
    v_rb = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nrotor * nbe, nrotor * nws, 2)
    l = lfs(s)
    v_rr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nrotor * nbe, nwn, 2)
    l = lfs(s)
    v_rw = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - Induced Velocities on Wakes - #
    s = (nwp, nbn, 2)
    l = lfs(s)
    v_wb = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwp, nrotor * nws, 2)
    l = lfs(s)
    v_wr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwp, nwn, 2)
    l = lfs(s)
    v_ww = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - Linear System - #

    s = (nbn + 2, nbn + 2)
    l = lfs(s)
    A_bb = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbn + 2,)
    l = lfs(s)
    b_bf = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbp, nwn)
    l = lfs(s)
    A_bw = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (2, nwn)
    l = lfs(s)
    A_pw = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbp, nrotor * nws)
    l = lfs(s)
    A_br = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (2, nrotor * nws)
    l = lfs(s)
    A_pr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - Blade Elements - #
    s = (nrotor,)
    l = lfs(s)
    B = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Rtip = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Rhub = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    fliplift = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbe, nrotor)
    l = lfs(s)
    chords = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    twists = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    stagger = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    solidity = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_panel_centers = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    inner_fraction = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - Wake constants - #
    s = (nwn,)
    l = lfs(s)
    wakeK = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # TODO: WHAT TO DO ABOUT AIRFOILS?
    # NOTE: how to know the size of airfoil objects up front?  if it's a DFDC type, this is pretty simple, but if it's a CCBlade type, then the size could be anything depending on how many data points are used
    # airfoils = allocate_airfoil_cache(aftype, nrotor, nbe)

    return (;
        solve_parameter_cache=PreallocationTools.DiffCache(
            zeros(total_length), fd_chunk_size; levels=levels
        ),
        solve_parameter_cache_dims=(;
            state_dims,
            vz_rotor,
            vtheta_rotor,
            Cm_wake,
            operating_point=(; Vinf, rhoinf, muinf, asound, Omega),
            ivr=(; v_rb, v_rr, v_rw),
            ivw=(; v_wb, v_wr, v_ww),
            linsys=(; A_bb, b_bf, A_bw, A_pw, A_br, A_pr),
            blade_elements=(;
                B,
                fliplift,
                chords,
                twists,
                stagger,
                solidity,
                rotor_panel_centers,
                inner_fraction,
                Rtip,
                Rhub,
            ),
            wakeK,
            SIAMFANLE_cache_vecs...,
        ),
    )
end

"""
"""
function allocate_solve_container_cache(
    solve_type::CSORSolverOptions, paneling_constants::PanelingConstants; fd_chunk_size=12, levels=1
)
    problem_dimensions = get_problem_dimensions(paneling_constants)

    return allocate_solve_container_cache(solve_type, problem_dimensions; fd_chunk_size=fd_chunk_size, levels=levels)
end

"""
"""
function allocate_solve_container_cache(
    solve_type::CSORSolverOptions, problem_dimensions; fd_chunk_size=12, levels=1
)
    (;
        nrotor,  # number of rotors
        nbodies, # number of bodies
        nwn,     # number of wake nodes
        nwp,     # number of wake panels
        nbn,     # number of body nodes
        nbe,     # number of blade elements (also rotor panels)
        nws,     # number of wake sheets (also rotor panel edges)
    ) = problem_dimensions

    # - initialize - #
    total_length = 0

    # Strengths
    s = (nbn + nbodies,)
    l = lfs(s)
    gamb = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l
    rhs = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # Blade Element Values
    s = (nbe, nrotor)
    l = lfs(s)
    beta1 = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Cz_rotor = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Ctheta_rotor = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Cmag_rotor = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    cl = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    cd = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    alpha = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    reynolds = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    mach = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vz_rotor = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vtheta_rotor = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # Circulation
    Gamma_tilde = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    H_tilde = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbe + 1, nrotor)
    l = lfs(s)
    deltaGamma2 = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    deltaH = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # Wake Velocities
    s = (nwp,)
    l = lfs(s)
    vz_wake = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vr_wake = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # State estimates
    s = (nbe, nrotor)
    l = lfs(s)
    Gamr_est = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nws, nrotor)
    l = lfs(s)
    sigr_est = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwn,)
    l = lfs(s)
    gamw_est = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwn,)
    l = lfs(s)
    Cm_avg = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwp,)
    l = lfs(s)
    Cm_wake = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbe, nrotor)
    l = lfs(s)
    deltaG = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    deltaG_prev= (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwn,)
    l = lfs(s)
    deltag = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    deltag_prev = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # Convergence criteria
    s = (nrotor,)
    l = lfs(s)
    maxBGamr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    maxdeltaBGamr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (1,)
    l = lfs(s)
    maxdeltagamw = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    return (;
        solve_container_cache=PreallocationTools.DiffCache(
            zeros(total_length), fd_chunk_size; levels=levels
        ),
        solve_container_cache_dims=(;
            gamb,
            rhs,
            Gamr_est,
            sigr_est,
            gamw_est,
            beta1,
            alpha,
            reynolds,
            mach,
            vz_rotor,
            vtheta_rotor,
            Cz_rotor,
            Ctheta_rotor,
            Cmag_rotor,
            cl,
            cd,
            vz_wake,
            vr_wake,
            Cm_wake,
            Cm_avg,
            Gamma_tilde,
            H_tilde,
            deltaGamma2,
            deltaH,
            deltaG,
            deltaG_prev,
            deltag,
            deltag_prev,
            maxBGamr,
            maxdeltaBGamr,
            maxdeltagamw,
        ),
    )
end

"""
"""
function allocate_solve_container_cache(
    solve_type::TS, paneling_constants::PanelingConstants; fd_chunk_size=12, levels=1
) where {TS<:Union{ExternalSolverOptions,MultiSolverOptions}}
    problem_dimensions = get_problem_dimensions(paneling_constants)

return allocate_solve_container_cache(
    solve_type, problem_dimensions; fd_chunk_size=fd_chunk_size, levels=levels
)
end

"""
"""
function allocate_solve_container_cache(
    solve_type::TS, problem_dimensions; fd_chunk_size=12, levels=1
) where {TS<:Union{ExternalSolverOptions,MultiSolverOptions}}

    (;
        nrotor, # number of rotors
        nwn,    # number of wake nodes
        nwp,    # number of wake panels
        nbn,    # number of body nodes
        nbe,    # number of blade elements (also rotor panels)
    ) = problem_dimensions

    # - initialize - #
    total_length = 0

    # Strengths
    # TODO: is this always going to be 2?, may want to add nb for nbodies to problem dims and then do +nb
    s = (nbn + 2,)
    l = lfs(s)
    gamb = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l
    rhs = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbe, nrotor)
    l = lfs(s)
    Gamr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbe + 1, nrotor)
    l = lfs(s)
    sigr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwn,)
    l = lfs(s)
    gamw = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # Blade Element Values
    s = (nbe, nrotor)
    l = lfs(s)
    beta1 = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Cz_rotor = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Ctheta_rotor = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Cmag_rotor = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    cl = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    cd = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    alpha = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    reynolds = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    mach = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # Circulation
    Gamma_tilde = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    H_tilde = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbe + 1, nrotor)
    l = lfs(s)
    deltaGamma2 = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    deltaH = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # Wake Velocities
    s = (nwp,)
    l = lfs(s)
    vz_wake = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vr_wake = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # State estimates
    s = (nbe, nrotor)
    l = lfs(s)
    vz_est = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vtheta_est = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwp,)
    l = lfs(s)
    Cm_est = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwn,)
    l = lfs(s)
    Cm_avg = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    return (;
        solve_container_cache=PreallocationTools.DiffCache(
            zeros(total_length), fd_chunk_size; levels=levels
        ),
        solve_container_cache_dims=(;
            gamb,
            rhs,
            Gamr,
            sigr,
            gamw,
            beta1,
            alpha,
            reynolds,
            mach,
            Cz_rotor,
            Ctheta_rotor,
            Cmag_rotor,
            cl,
            cd,
            vz_wake,
            vr_wake,
            vz_est,
            vtheta_est,
            Cm_est,
            Cm_avg,
            Gamma_tilde,
            H_tilde,
            deltaGamma2,
            deltaH,
        ),
    )
end