include("../../test/test_utils.jl")
######################################################################
#                                                                    #
#                        FUNCTION SET HEADER                         #
#                                                                    #
######################################################################

"""
TODO: move to analysis.jl
"""
function analyze(
    propulsor;
    options=set_options(),
    # precomp_container_caching=nothing,
    # solve_parameter_caching=nothing,
    solve_container_caching=nothing,
    # post_caching=nothing
)

    # if isnothing(precomp_container_caching)
    #     precomp_container_caching = generate_precomp_container_cache(
    #         propulsor.paneling_constants
    #     )
    # end
    # if isnothing(solve_parameter_caching)
    #     solve_parameter_caching = generate_solve_parameter_cache(
    #         propulsor.paneling_constants
    #     )
    # end
    # if isnothing(post_caching)
    #     post_caching = generate_post_cache(
    #         propulsor.paneling_constants
    #     )
    # end
    if isnothing(solve_container_caching)
        solve_container_caching = allocate_solve_container_cache(
            propulsor.paneling_constants
        )
    end

    # TODO: write generate_caches function (need to go through and figure out what all goes in the cache)
    # - Pull out the Caches - #
    # (; precomp_container_cache, precomp_container_cache_dims) = precomp_container_caching
    # (; post_cache, post_cache_dims) = post_caching

    # # - Reshape precomp_container_cache - #
    # # TODO: get a different type to dispatch on here
    # precomp_container_cache_vec = @views pat.get_tmp(precomp_container_cache, TF(1.0))
    # # TODO: test this function
    # precomp_containers = withdraw_precomp_container_cache(
    #     precomp_container_cache, precomp_container_cache_dims
    # )

    # # - Reshape solve_parameter_cache - #
    # solve_parameter_cache_vec = @views pat.get_tmp(solve_parameter_cache, inputs)
    # # TODO: test this function
    # solve_parameter_containers = withdraw_solve_parameter_cache(
    #     solve_parameter_cache, solve_parameter_cache_dims
    # )

    # - Do precomputations - #
    # out-of-place version currently has 22,292,181 allocations.
    ivr, ivw, ivb, linsys, blade_elements, wakeK, idmaps, panels, problem_dimensions = precompute_parameters_iad(
        propulsor;
        wake_nlsolve_ftol=options.wake_nlsolve_ftol,
        wake_max_iter=options.wake_max_iter,
        max_wake_relax_iter=options.max_wake_relax_iter,
        wake_relax_tol=options.wake_relax_tol,
        autoshiftduct=options.autoshiftduct,
        itcpshift=options.itcpshift,
        axistol=options.axistol,
        tegaptol=options.tegaptol,
        finterp=options.finterp,
        silence_warnings=options.silence_warnings,
        verbose=options.verbose,
    )

    if iszero(linsys.lu_decomp_flag[1])
        if !silence_warnings
            @warn "Exiting.  LU decomposition of the LHS matrix for the linear system failed.  Please check your body geometry and ensure that there will be no panels lying directly atop eachother or other similar problematic geometry."
            #TODO: write this function that returns the same as outs below, but all zeros
            return zero_outputs(), false
        end
    end

    # - Initialize Aero - #
    # TODO; add some sort of unit test for this function
    # note: if using a cache for intermediate calcs here, doesn't need to be fancy.
    # out-of-place version currently has 2990 allocations
    vz_rotor, vtheta_rotor, Cm_wake = initialize_velocities(
        propulsor.operating_point,
        blade_elements,
        linsys,
        ivr,
        ivw,
        idmaps.body_totnodes,
        idmaps.wake_panel_sheet_be_map,
    )

    # - Vectorize Inputs - #
    inputs, state_dims = vectorize_inputs(vz_rotor, vtheta_rotor, Cm_wake)

    # - combine cache and constants - #
    const_cache = (;
        # - Constants - #
        options.verbose,
        options.silence_warnings,
        #nlsolve options
        options.nlsolve_method,
        options.nlsolve_autodiff,
        options.nlsolve_linesearch_method,
        options.nlsolve_linesearch_kwargs,
        options.nlsolve_ftol,
        options.nlsolve_iteration_limit,
        options.nlsolve_converged,
        # Parameters
        propulsor.operating_point,
        ivr,
        ivw,
        linsys,
        blade_elements,
        wakeK,
        idmaps,
        state_dims,
        # Cache(s)
        solve_container_caching...,
    )

    # - Initialize Aero and Converge Gamr and gamw - #
    velocity_states = ImplicitAD.implicit(solve_iad, residual!, inputs, const_cache)

    # - Post-Process - #
    # NOTE: post processing cache doesn't need to be fancy, since it's just here in the analysis and will always be the same type if derivative checks are turned off.
    # update tests for this function
    outs = post_process_iad!(
        solve_container_caching,
        velocity_states,
        state_dims,
        propulsor.operating_point,
        propulsor.reference_parameters,
        ivr,
        ivw,
        ivb,
        linsys,
        blade_elements,
        wakeK,
        panels.body_vortex_panels,
        panels.rotor_source_panels.influence_length,
        idmaps,
        problem_dimensions;
        write_outputs=options.write_outputs,
        outfile=options.outfile,
        checkoutfileexists=options.checkoutfileexists,
        tuple_name=options.tuple_name,
    )

    return outs, true
    # return velocity_states
end

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

    result = NLsolve.nlsolve(
        rwrap!,
        inputs; # initial states guess
        method=nlsolve_method,
        autodiff=nlsolve_autodiff,
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
        @view(solve_containers.Gamr[:, :]),
        solve_containers.Cmag_rotor,
        blade_elements.chords,
        solve_containers.cl,
    )

    # - Calculate Rotor Panel Strengths - #
    # currently has 122 allocations
    calculate_rotor_source_strengths!(
        @view(solve_containers.sigr[:, :]),
        solve_containers.Cmag_rotor,
        blade_elements.chords,
        blade_elements.B,
        solve_containers.cd,
    )

    # - Get Average Wake Velocities - #
    # currently has 5 allocations
    average_wake_velocities!(
        @view(solve_containers.Cm_avg[:]),
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
        idmaps.wake_panel_sheet_be_map,
        idmaps.wake_node_ids_along_casing_wake_interface,
        idmaps.wake_node_ids_along_centerbody_wake_interface;
    )

    # - Solve Linear System for Body Strengths - #
    # currently has 18 allocations
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
        operating_point.Vinf,
    )

    # return estimated states
    return vz_est, vtheta_est, Cm_est
end
