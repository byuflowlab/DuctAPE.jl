"""
"""
function analyze(
    propulsor::Propulsor,
    options=set_options();
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
        wake_solve_options=options.wake_options,
        # wake_nlsolve_ftol=options.wake_nlsolve_ftol,
        # wake_max_iter=options.wake_max_iter,
        # max_wake_relax_iter=options.max_wake_relax_iter,
        # wake_relax_tol=options.wake_relax_tol,
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

    return analyze(
        propulsor,
        ivr,
        ivw,
        ivb,
        linsys,
        blade_elements,
        wakeK,
        idmaps,
        panels,
        problem_dimensions,
        options;
        # precomp_container_caching=nothing,
        # solve_parameter_caching=nothing,
        solve_container_caching=solve_container_caching,
        # post_caching=nothing
    )
end

"""
"""
function analyze(
    propulsor::Propulsor,
    ivr,
    ivw,
    ivb,
    linsys,
    blade_elements,
    wakeK,
    idmaps,
    panels,
    problem_dimensions,
    options=set_options();
    # precomp_container_caching=nothing,
    # solve_parameter_caching=nothing,
    solve_container_caching=nothing,
    # post_caching=nothing
)

    # - Finish Pre-Processing - #
    # Set up Solve Container Cache
    if isnothing(solve_container_caching)
        solve_container_caching = allocate_solve_container_cache(
            options.solve_options, propulsor.paneling_constants
        )
    end

    # - Process - #
    velocity_states = process(
        options.solve_options,
        propulsor::Propulsor,
        ivr,
        ivw,
        ivb,
        linsys,
        blade_elements,
        wakeK,
        idmaps,
        panels,
        problem_dimensions,
        options,
        solve_container_caching,
    )

    # - Post-Process - #
    # NOTE: post processing cache doesn't need to be fancy, since it's just here in the analysis and will always be the same type if derivative checks are turned off.
    # update tests for this function

    return outs, (; ivr, ivw, ivb, linsys, blade_elements, wakeK, idmaps, panels), true
end

"""
"""
function process(
    process_type::NewtonSolve,
    propulsor::Propulsor,
    ivr,
    ivw,
    ivb,
    linsys,
    blade_elements,
    wakeK,
    idmaps,
    panels,
    problem_dimensions,
    options,
    solve_container_caching,
)

    # - Rename for Convenience - #
    (; solve_options) = options

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
        solve_options.nlsolve_method,
        solve_options.nlsolve_autodiff,
        solve_options.nlsolve_linesearch_method,
        solve_options.nlsolve_linesearch_kwargs,
        solve_options.nlsolve_ftol,
        solve_options.nlsolve_iteration_limit,
        solve_options.nlsolve_converged,
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
    return ImplicitAD.implicit(solve_iad, residual!, inputs, const_cache)
end

"""
"""
function process(
    process_type::QuickSolve,
    propulsor::Propulsor,
    ivr,
    ivw,
    ivb,
    linsys,
    blade_elements,
    wakeK,
    idmaps,
    panels,
    problem_dimensions,
    options,
    solve_container_caching,
)

    # - Rename for Convenience - #
    (; verbose, solve_options) = options

    # - Initialize States - #

    # vectorize inputs

    # - TODO: Set up const_cache for this function - #

    #TODO: update this function
    # need to withdraw cache inside this function to have the *_est values, as well as sigr and gamb available.
    return solve(
        inputs,
        const_cache;
        verbose=verbose,
        solve_options.nosource,
        solve_options.maxiter,
        solve_options.nrf,
        solve_options.bt1,
        solve_options.bt2,
        solve_options.pf1,
        solve_options.pf2,
        solve_options.btw,
        solve_options.pfw,
        solve_options.f_circ,
        solve_options.f_dgamw,
    )
end

"""
"""
function post_process(
    process_type::NewtonSolve,
    solve_container_caching,
    states,
    state_dims,
    operating_point,
    reference_parameters,
    ivr,
    ivw,
    ivb,
    linsys,
    blade_elements,
    wakeK,
    body_vortex_panels,
    influence_length,
    idmaps,
    problem_dimensions;
    write_outputs=false,
    outfile="outputs.jl",
    checkoutfileexists=false,
    output_tuple_name="outs",
)
    return post_process_iad!(
        solve_container_caching,
        states,
        state_dims,
        operating_point,
        reference_parameters,
        ivr,
        ivw,
        ivb,
        linsys,
        blade_elements,
        wakeK,
        body_vortex_panels,
        influence_length,
        idmaps,
        problem_dimensions;
        write_outputs=options.write_outputs,
        outfile=options.outfile,
        checkoutfileexists=options.checkoutfileexists,
        output_tuple_name=options.output_tuple_name,
    )
end

"""
"""
function post_process(
    process_type::QuickSolve,
    solve_container_caching,
    states,
    state_dims,
    operating_point,
    reference_parameters,
    ivr,
    ivw,
    ivb,
    linsys,
    blade_elements,
    wakeK,
    body_vortex_panels,
    influence_length,
    idmaps,
    problem_dimensions;
    write_outputs=false,
    outfile="outputs.jl",
    checkoutfileexists=false,
    output_tuple_name="outs",
)
    #TODO: update this function
    return post_process()
end
