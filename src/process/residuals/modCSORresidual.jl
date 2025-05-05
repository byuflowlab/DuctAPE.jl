"""
    mod_CSOR_residual!(r, current_states, inputs, constants)

Modified DFDC-like CSOR residual that does not include any relaxation within the residual calculation.

# Arguments
- `r::Vector{Float}` : solve residual
- `current_states::Vector{Float}` : solve states
- `inputs::Vector{Float}` : solve inputs
- `constants::Vector{Float}` : solve constants
"""
function mod_CSOR_residual!(r, current_states, inputs, constants)

    # - Extract constants - #
    (;
        # solve options for dispatch
        solver_options,
        airfoils,                   # airfoils
        A_bb_LU,                    # linear system left hand side LU decomposition
        idmaps,                     # book keeping items
        solve_parameter_cache_dims, # dimensions for shaping the view of the parameter cache
        solve_container_cache,      # cache for solve_containers used in solve
        solve_container_cache_dims, # dimensions for shaping the view of the solve cache
    ) = constants

    # - Separate out the state variables - #
    Gamr, sigr, gamw = extract_state_variables(
        solver_options, current_states, solve_parameter_cache_dims.state_dims
    )

    # separate out sensitivity_parameters here as well
    solve_parameter_tuple = withdraw_solve_parameter_cache(
        solver_options, inputs, solve_parameter_cache_dims
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
        solve_container_cache, promote_type(eltype(current_states), eltype(inputs))(1.0)
    )
    # zero out contents of solve_containers to avoid any potential contamination issues
    solve_container_cache_vector .= 0
    solve_containers = withdraw_solve_container_cache(
        solver_options, solve_container_cache_vector, solve_container_cache_dims
    )

    estimated_states = estimate_CSOR_states!(
        solve_containers,
        Gamr,
        sigr,
        gamw,
        operating_point,
        ivr,
        ivw,
        (;linsys..., A_bb_LU),
        (;blade_elements..., airfoils...),
        wakeK,
        idmaps;
        verbose=false,
    )

    @. r = estimated_states - current_states

    return r
end

"""
    estimate_CSOR_states!(
        solve_containers,
        Gamr,
        sigr,
        gamw,
        operating_point,
        ivr,
        ivw,
        linsys,
        blade_elements,
        wakeK,
        idmaps;
        verbose=false,
    )

Estimate states for modified CSOR solver.

# Arguments:
- `solve_containers::NamedTuple` : cache for intermediate solve values
- `Gamr::type` : Blade element circulation strengths
- `sigr::type` : Rotor source panel strengths
- `gamw::type` : Wake vortex panel strengths
- `operating_point::NamedTuple` : Named tuple containing operating_point information
- `ivr::NamedTuple` : unit induced velocities on rotor(s)
- `ivw::NamedTuple` : unit induced velocities on wake
- `linsys::NamedTuple` : vectors and matricies comprising the panel method linear system
- `blade_elements::NamedTuple` : blade element geometry and airfoil polar information
- `wakeK::Vector{Float}` : geometric constants used in caculating wake strengths
- `idmaps::NamedTuple` : index maps used throughout solve
"""
function estimate_CSOR_states!(
    solve_containers,
    Gamr,
    sigr,
    gamw,
    operating_point,
    ivr,
    ivw,
    linsys,
    blade_elements,
    wakeK,
    idmaps;
    verbose=false,
)

    # - Solve Linear System - #
    # in place solve for gamb
    calculate_body_vortex_strengths!(
        solve_containers.gamb,
        linsys.A_bb_LU,
        linsys.b_bf,
        gamw,
        linsys.A_bw,
        linsys.A_pw,
        sigr,
        linsys.A_br,
        linsys.A_pr,
        linsys.A_bb,
        solve_containers.rhs;
        post=false,
    )

    # - Update rotor blade element velocities with body influence - #
    calculate_induced_velocities_on_rotors!(
        solve_containers.vz_rotor,
        solve_containers.vtheta_rotor,
        Gamr,
        gamw,
        sigr,
        @view(solve_containers.gamb[1:(idmaps.body_totnodes)]),
        @view(ivr.v_rw[:, :, 1]),
        @view(ivr.v_rr[:, :, 1]),
        @view(ivr.v_rb[:, :, 1]),
        blade_elements.B,
        blade_elements.rotor_panel_centers,
    )

    # - Get Absolute Rotor Velocities - #
    reframe_rotor_velocities!(
        solve_containers.Cz_rotor,
        solve_containers.Ctheta_rotor,
        solve_containers.Cmag_rotor,
        solve_containers.vz_rotor,
        solve_containers.vtheta_rotor,
        operating_point.Vinf[],
        operating_point.Omega,
        blade_elements.rotor_panel_centers,
    )

    # - Calculate Blade Element Values - #
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

    ##### ----- Estimate Gamr ----- #####

    # - Calculate Blade Element Circulation - #
    calculate_rotor_circulation_strengths!(
        solve_containers.Gamr_est,
        solve_containers.Cmag_rotor,
        blade_elements.chords,
        solve_containers.cl,
        blade_elements
    )

    ##### ----- Estimate gamw ----- #####

    # - Calculate Velocities on Wake Panels - #
    calculate_wake_velocities!(
        solve_containers.Cm_wake,
        solve_containers.vz_wake,
        solve_containers.vr_wake,
        gamw,
        sigr,
        @view(solve_containers.gamb[1:(idmaps.body_totnodes)]),
        ivw,
        operating_point.Vinf[],
    )

    # - Get Average Wake Velocities - #
    average_wake_velocities!(
        solve_containers.Cm_avg,
        solve_containers.Cm_wake,
        idmaps.wake_nodemap,
        idmaps.wake_endnodeidxs,
    )

    # - Estimate wake strengths - #
    # in-place solve for gamw_est
    calculate_wake_vortex_strengths!(
        solve_containers.gamw_est,
        solve_containers.Gamma_tilde,
        solve_containers.H_tilde,
        solve_containers.deltaGamma2,
        solve_containers.deltaH,
        Gamr,
        solve_containers.Cm_avg,
        blade_elements.B,
        operating_point.Omega,
        wakeK,
        idmaps.wake_node_sheet_be_map,
        idmaps.wake_node_ids_along_casing_wake_interface,
        idmaps.wake_node_ids_along_centerbody_wake_interface;
    )

    ##### ----- Update sigr ----- #####
    # - Update rotor blade element velocities with body influence - #
    calculate_induced_velocities_on_rotors!(
        solve_containers.vz_rotor,
        solve_containers.vtheta_rotor,
        Gamr,
        gamw,
        sigr,
        @view(solve_containers.gamb[1:(idmaps.body_totnodes)]),
        @view(ivr.v_rw[:, :, 1]),
        @view(ivr.v_rr[:, :, 1]),
        @view(ivr.v_rb[:, :, 1]),
        blade_elements.B,
        blade_elements.rotor_panel_centers,
    )

    # - Get Absolute Rotor Velocities - #
    reframe_rotor_velocities!(
        solve_containers.Cz_rotor,
        solve_containers.Ctheta_rotor,
        solve_containers.Cmag_rotor,
        solve_containers.vz_rotor,
        solve_containers.vtheta_rotor,
        operating_point.Vinf[],
        operating_point.Omega,
        blade_elements.rotor_panel_centers,
    )

    # - Calculate Rotor Panel Strengths - #
    calculate_rotor_source_strengths!(
        solve_containers.sigr_est,
        solve_containers.Cmag_rotor,
        blade_elements.chords,
        blade_elements.B,
        solve_containers.cd,
    )

    return [
        reshape(solve_containers.Gamr_est, :)
        reshape(solve_containers.sigr_est, :)
        reshape(solve_containers.gamw_est, :)
    ]
end
