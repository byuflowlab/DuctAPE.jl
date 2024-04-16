"""
"""
function system_residual(state_variables, sensitivity_parameters, constants)
    resid = similar(state_variables) .= 0
    return system_residual!(resid, state_variables, sensitivity_parameters, constants)
end

"""
"""
function system_residual!(resid, state_variables, sensitivity_parameters, constants)

    # - Extract constants - #
    (;
        # solve options for dispatch
        solver_options,
        # comp,
        airfoils,                   # airfoils
        A_bb_LU,                    # linear system left hand side LU decomposition
        idmaps,                     # book keeping items
        solve_parameter_cache_dims, # dimensions for shaping the view of the parameter cache
        solve_container_cache,      # cache for solve_containers used in solve
        solve_container_cache_dims, # dimensions for shaping the view of the solve cache
    ) = constants

    # - Separate out the state variables - #
    # TODO: test this function
    vz_rotor, vtheta_rotor, Cm_wake = extract_state_variables(
        solver_options, state_variables, solve_parameter_cache_dims.state_dims
    )

    # separate out sensitivity_parameters here as well
    solve_parameter_tuple = withdraw_solve_parameter_cache(
        solver_options, sensitivity_parameters, solve_parameter_cache_dims
    )

    (;
        operating_point, # freestream, Omega, etc.
        ivr,             # induced velocities on rotor panels
        ivw,             # induced velocities on wake panels
        linsys,          # includes AIC's for linear system
        blade_elements,  # includes blade element geometry
        wakeK,           # geometric "constants" for wake node strength calculation
    ) = solve_parameter_tuple

    # - Extract and Reset Cache - #
    # get cache vector of correct types
    solve_container_cache_vector = @views PreallocationTools.get_tmp(
        solve_container_cache,
        promote_type(eltype(state_variables), eltype(sensitivity_parameters))(1.0),
    )
    # zero out contents of solve_containers to avoid any potential contamination issues
    solve_container_cache_vector .= 0
    # TODO: test this function
    solve_containers = withdraw_solve_container_cache(
        solver_options, solve_container_cache_vector, solve_container_cache_dims
    )

    # - Estimate New States - #
    # TODO: test this function
    estimate_states!(
        solve_containers,
        vz_rotor,
        vtheta_rotor,
        Cm_wake,
        operating_point,
        ivr,
        ivw,
        (; linsys..., A_bb_LU),
        (; blade_elements..., airfoils...),
        wakeK,
        idmaps,
    )

    # - Get final Residual Values - #
    return update_system_residual!(
        solver_options,
        resid,
        solve_containers.vz_est,
        vz_rotor,
        solve_containers.vtheta_est,
        vtheta_rotor,
        solve_containers.Cm_est,
        Cm_wake,
        solve_parameter_cache_dims,
    )
end

"""
"""
function update_system_residual!(
    solver_options::Union{NonlinearSolveOptions,NLsolveOptions,SIAMFANLEOptions,MinpackOptions},
    resid,
    vz_est,
    vz_rotor,
    vtheta_est,
    vtheta_rotor,
    Cm_est,
    Cm_wake,
    solve_parameter_cache_dims,
)
    vz_est .-= vz_rotor
    vtheta_est .-= vtheta_rotor
    Cm_est .-= Cm_wake

    resid[solve_parameter_cache_dims.state_dims.vz_rotor.index] .= reshape(vz_est, :)
    resid[solve_parameter_cache_dims.state_dims.vtheta_rotor.index] .= reshape(
        vtheta_est, :,
    )
    resid[solve_parameter_cache_dims.state_dims.Cm_wake.index] .= Cm_est

    return resid
end

"""
"""
function update_system_residual!(
    solver_options::Union{SpeedMappingOptions,FixedPointOptions},
    resid,
    vz_est,
    vz_rotor,
    vtheta_est,
    vtheta_rotor,
    Cm_est,
    Cm_wake,
    solve_parameter_cache_dims,
)
    resid[solve_parameter_cache_dims.state_dims.vz_rotor.index] .= reshape(vz_est, :)
    resid[solve_parameter_cache_dims.state_dims.vtheta_rotor.index] .= reshape(
        vtheta_est, :,
    )
    resid[solve_parameter_cache_dims.state_dims.Cm_wake.index] .= Cm_est

    return resid
end

"""
"""
function estimate_states!(
    solve_containers,
    vz_rotor,
    vtheta_rotor,
    Cm_wake,
    operating_point,
    ivr,
    ivw,
    linsys,
    blade_elements,
    wakeK,
    idmaps;
    verbose=false,
)

    # - Rename for Convenience - #
    vz_est = solve_containers.vz_est
    vtheta_est = solve_containers.vtheta_est
    Cm_est = solve_containers.Cm_est

    # - Get Absolute Rotor Velocities - #
    # currently has 9 allocations
    reframe_rotor_velocities!(
        solve_containers.Cz_rotor,
        solve_containers.Ctheta_rotor,
        solve_containers.Cmag_rotor,
        vz_rotor, # state var
        vtheta_rotor, # state var
        operating_point.Vinf[1],
        operating_point.Omega,
        blade_elements.rotor_panel_centers,
    )

    # - Calculate Blade Element Values - #
    # currently has 141 allocations using DFDC airfoils
    calculate_blade_element_coefficients!(
        solve_containers.cl,
        solve_containers.cd,
        solve_containers.beta1,
        solve_containers.alpha,
        solve_containers.reynolds,
        solve_containers.mach,
        blade_elements,
        solve_containers.Cz_rotor,
        solve_containers.Ctheta_rotor,
        solve_containers.Cmag_rotor,
        operating_point;
        verbose=verbose,
    )

    # - Calculate Blade Element Circulation - #
    # currently has 9 allocations
    calculate_rotor_circulation_strengths!(
        solve_containers.Gamr,
        solve_containers.Cmag_rotor,
        blade_elements.chords,
        solve_containers.cl,
    )

    # - Calculate Rotor Panel Strengths - #
    # currently has 122 allocations
    calculate_rotor_source_strengths!(
        solve_containers.sigr,
        solve_containers.Cmag_rotor,
        blade_elements.chords,
        blade_elements.B,
        solve_containers.cd,
    )

    # - Get Average Wake Velocities - #
    # currently has 5 allocations
    average_wake_velocities!(
        solve_containers.Cm_avg,
        Cm_wake, # state var
        idmaps.wake_nodemap,
        idmaps.wake_endnodeidxs,
    )

    # - Calculate Wake Panel Strengths - #
    # in-place solve for gamw,
    # currently has 27 allocations
    calculate_wake_vortex_strengths!(
        solve_containers.gamw,
        solve_containers.Gamma_tilde,
        solve_containers.H_tilde,
        solve_containers.deltaGamma2,
        solve_containers.deltaH,
        solve_containers.Gamr,
        solve_containers.Cm_avg,
        blade_elements.B,
        operating_point.Omega,
        wakeK,
        idmaps.wake_node_sheet_be_map,
        idmaps.wake_node_ids_along_casing_wake_interface,
        idmaps.wake_node_ids_along_centerbody_wake_interface;
    )

    # - Solve Linear System for Body Strengths - #
    # currently has 18 allocations
    # the gamw and sigr match pretty well (there's likely an indexing difference in gamw, I think the old one was wrong)
    calculate_body_vortex_strengths!(
        solve_containers.gamb,
        linsys.A_bb_LU,
        linsys.b_bf,
        solve_containers.gamw,
        linsys.A_bw,
        linsys.A_pw,
        solve_containers.sigr,
        linsys.A_br,
        linsys.A_pr,
        linsys.A_bb,
        solve_containers.rhs,
    )

    # - Calcuate vz_est and vtheta_est- #
    # TODO: test this function
    # currently has 23 allocations
    calculate_induced_velocities_on_rotors!(
        vz_est,
        vtheta_est,
        solve_containers.Gamr,
        solve_containers.gamw,
        solve_containers.sigr,
        @view(solve_containers.gamb[1:(idmaps.body_totnodes)]),
        @view(ivr.v_rw[:, :, 1]),
        @view(ivr.v_rr[:, :, 1]),
        @view(ivr.v_rb[:, :, 1]),
        blade_elements.B,
        blade_elements.rotor_panel_centers,
    )

    # - Calculate Velocities on Wake Panels - #
    # TODO: test this function
    # currently has 23 allocations
    calculate_wake_velocities!(
        Cm_est,
        solve_containers.vz_wake,
        solve_containers.vr_wake,
        solve_containers.gamw,
        solve_containers.sigr,
        @view(solve_containers.gamb[1:(idmaps.body_totnodes)]),
        ivw,
        operating_point.Vinf[],
    )

    # return estimated states
    return vz_est, vtheta_est, Cm_est
end
