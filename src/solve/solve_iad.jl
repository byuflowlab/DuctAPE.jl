#---------------------------------#
#              SOLVE              #
#---------------------------------#

"""
"""
function solve_iad(inputs, const_cache)

    ##### ----- Do everything in here ----- #####
    # Note: put all the precomputation in here, but not in the nlsolve residual.
    # Note: it will have to be done twice, but oh well.
    # - Extract constants - #
    (;
        #general
        verbose,
        #nlsolve options
        nlsolve_method,
        nlsolve_autodiff,
        nlsolve_linesearch_method,
        nlsolve_linesearch_kwargs,
        nlsolve_ftol,
        nlsolve_iteration_limit,
        nlsolve_converged,
        # Parameters
        operating_point,
        ivr,
        ivw,
        linsys,
        blade_elements,
        wakeK,
        idmaps,
        state_dims,
        #caches
        # solve_parameter_cache,
        # solve_parameter_cache_dims,
        solve_container_cache,
        solve_container_cache_dims,
    ) = const_cache

    # - Wrap residual - #
    function rwrap!(r, state_variables)
        #TODO: Test this function
        return residual!(
            r,
            state_variables,
            (;
                operating_point,# includes freestream and Omega
                ivr,            # induced velocities on rotor panels
                ivw,            # induced velocities on wake panels
                linsys,         # includes AIC's for linear system
                blade_elements, # includes blade element geometry
                wakeK,          # Geometry-based constant for wake strength calculation
                idmaps,         # book keeping items
                state_dims,     # dimensions for state variable extraction
                solve_container_cache,      # cache for solve_containers used in solve
                solve_container_cache_dims, # dimensions for shaping the view of the solve cache
            ),
        )
    end

    # - Wrap Jacobian - #
    # configure jacobian
    jconfig = ForwardDiff.JacobianConfig(
        rwrap!,
        inputs, # object of size of residual vector
        inputs, # object of size of input vector
        ForwardDiff.Chunk{12}(), #TODO: set the chunk size as an option and make sure the solver cache chunk and this chunk size match. PreallocationTools chooses poorly if not given a chunk size
    )
    # get allocated array for jacobian
    Jwork = DiffResults.JacobianResult(
        inputs, # object of size of residual vector
        inputs, # object of size of input vector
    )
    # store everything in a cache that can be accessed inside the solver
    jcache = (; rwrap!, Jwork, config=jconfig, r=zeros(size(inputs)))

    # wrap the jacobian
    function jwrap!(J, state_variables)
        ForwardDiff.jacobian!(
            jcache.Jwork, jcache.rwrap!, jcache.r, state_variables, jcache.config
        );
        J = DiffResults.jacobian(jcache.Jwork);
        return J
    end

    df = NLsolve.OnceDifferentiable(rwrap!, jwrap!, inputs, similar(inputs))
    result = NLsolve.nlsolve(
        df,
        inputs; # initial states guess
        method=nlsolve_method,
        # autodiff=nlsolve_autodiff,
        linesearch=nlsolve_linesearch_method(nlsolve_linesearch_kwargs...),
        linsolve=(x, A, b) -> x .= ImplicitAD.implicit_linear(A, b),
        ftol=nlsolve_ftol,
        iterations=nlsolve_iteration_limit,
        show_trace=verbose,
    )

    # update convergence flag
    nlsolve_converged[1] = NLsolve.converged(result)

    return result.zero
end

#---------------------------------#
#             RESIDUALS           #
#---------------------------------#

"""
This is the residual that gets passed into nlsolve, and is just for converging Gamr and gamw.
"""
function residual!(r, state_variables, parameters)
    # - Extract Parameters - #
    (;
        operating_point,             # includes freestream and Omega
        ivr,            # induced velocities on rotor panels
        ivw,            # induced velocities on wake panels
        linsys,         # includes AIC's for linear system
        blade_elements, # includes blade element geometry
        wakeK,
        idmaps,         # book keeping items
        state_dims,     # dimensions for state variable extraction
        solve_container_cache,      # cache for solve_containers used in solve
        solve_container_cache_dims, # dimensions for shaping the view of the solve cache
    ) = parameters

    # - Separate out the state variables - #
    # TODO: test this function
    vz_rotor, vtheta_rotor, Cm_wake = extract_state_vars(state_variables, state_dims)

    # - Extract and Reset Cache - #
    # get cache vector of correct types
    solve_container_cache_vec = @views PreallocationTools.get_tmp(
        solve_container_cache, state_variables
    )
    # TODO: test this function
    solve_containers = withdraw_solve_container_cache(
        solve_container_cache_vec, solve_container_cache_dims
    )
    # zero out contents of solve_containers to avoid any potential contamination issues
    # TODO: test this function
    # Note: there are 32 allocations for this function, can that be reduced to zero?
    reset_containers!(solve_containers) #note: also zeros out state estimates

    # - Estimate New States - #
    # TODO: test this function
    estimate_states!(
        solve_containers.vz_est,
        solve_containers.vtheta_est,
        solve_containers.Cm_est,
        vz_rotor,
        vtheta_rotor,
        Cm_wake,
        solve_containers,
        operating_point,
        ivr,
        ivw,
        linsys,
        blade_elements,
        wakeK,
        idmaps,
    )

    # - Get final Residual Values - #
    solve_containers.vz_est .-= vz_rotor
    solve_containers.vtheta_est .-= vtheta_rotor
    solve_containers.Cm_est .-= Cm_wake

    r[state_dims.vz_rotor.index] .= reshape(solve_containers.vz_est, :)
    r[state_dims.vtheta_rotor.index] .= reshape(solve_containers.vtheta_est, :)
    r[state_dims.Cm_wake.index] .= solve_containers.Cm_est

    ##NOTE: can't use a smaller residual with NLsolve
    #r[1] = solve_containers.vz_est[findmax(abs.(solve_containers.vz_est))[2]]
    #r[2] = solve_containers.vtheta_est[findmax(abs.(solve_containers.vtheta_est))[2]]
    #r[3] = solve_containers.Cm_est[findmax(abs.(solve_containers.Cm_est))[2]]

    return r
end

"""
"""
function estimate_states!(
    vz_est,
    vtheta_est,
    Cm_est,
    vz_rotor,
    vtheta_rotor,
    Cm_wake,
    solve_containers,
    operating_point,
    ivr,
    ivw,
    linsys,
    blade_elements,
    wakeK,
    idmaps;
    verbose=false,
)

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
        operating_point.Vinf[1],
    )

    # return estimated states
    return vz_est, vtheta_est, Cm_est
end
