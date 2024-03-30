"""
"""
function analyze(
    propulsor::Propulsor,
    options=set_options();
    # precomp_container_caching=nothing,
    solve_parameter_caching=nothing,
    solve_container_caching=nothing,
    # post_caching=nothing
)

    # - Get type to dispatch caches - #
    TF = promote_type(
        eltype(propulsor.duct_coordinates),
        eltype(propulsor.centerbody_coordinates),
        eltype(propulsor.operating_point.Vinf),
        eltype(propulsor.operating_point.Omega),
        eltype(propulsor.operating_point.rhoinf),
        eltype(propulsor.operating_point.muinf),
        eltype(propulsor.operating_point.asound),
        eltype(propulsor.rotorstator_parameters.B),
        eltype(propulsor.rotorstator_parameters.Rhub),
        eltype(propulsor.rotorstator_parameters.Rtip),
        eltype(propulsor.rotorstator_parameters.rotorzloc),
        eltype(propulsor.rotorstator_parameters.chords),
        eltype(propulsor.rotorstator_parameters.twists),
    )

    # if isnothing(precomp_container_caching)
    #     precomp_container_caching = allocate_precomp_container_cache(
    #         propulsor.paneling_constants
    #     )
    # end

    # if isnothing(post_caching)
    #     post_caching = allocate_post_cache(
    #         propulsor.paneling_constants
    #     )
    # end

    # - Pull out the Caches - #
    # (; precomp_container_cache, precomp_container_cache_dims) = precomp_container_caching
    # (; post_cache, post_cache_dims) = post_caching

    # # - Reshape precomp_container_cache - #
    # precomp_container_cache_vec = @views PreallocationTools.get_tmp(precomp_container_cache, TF(1.0))
    # precomp_containers = withdraw_precomp_container_cache(
    #     precomp_container_cache, precomp_container_cache_dims
    # )

    # - Set up Solver Sensitivity Paramter Cache - #

    # Allocate Cache
    if isnothing(solve_parameter_caching)
        solve_parameter_caching = allocate_solve_parameter_cache(
            options.solve_options, propulsor.paneling_constants
        )
    end

    # separate out caching items
    (; solve_parameter_cache, solve_parameter_cache_dims) = solve_parameter_caching

    # get correct cache type
    solve_parameter_cache_vector = @views PreallocationTools.get_tmp(
        solve_parameter_cache, TF(1.0)
    )
    # reset cache
    solve_parameter_cache_vector .= 0

    # reshape cache
    solve_parameter_tuple = withdraw_solve_parameter_cache(
        options.solve_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )

    # copy over operating point
    for f in fieldnames(typeof(propulsor.operating_point))
        solve_parameter_tuple.operating_point[f] .= getfield(propulsor.operating_point, f)
    end

    # - Do precomputations - #
    if options.verbose
        println("Pre-computing Parameters")
    end
    # out-of-place version currently has 22,292,181 allocations.
    # TODO: do this in place for the solve input cache items. eventually will want to have a post-processing and output cache too.
    ivb, A_bb_LU, lu_decomp_flag, airfoils, idmaps, panels, problem_dimensions = precompute_parameters_iad!(
        solve_parameter_tuple.ivr,
        solve_parameter_tuple.ivw,
        solve_parameter_tuple.blade_elements,
        solve_parameter_tuple.linsys,
        solve_parameter_tuple.wakeK,
        propulsor;
        wake_solve_options=options.wake_options,
        autoshiftduct=options.autoshiftduct,
        itcpshift=options.itcpshift,
        axistol=options.axistol,
        tegaptol=options.tegaptol,
        finterp=options.finterp,
        silence_warnings=options.silence_warnings,
        verbose=options.verbose,
    )

    # - Check that the precomputation went well - #
    #=
      NOTE: If the linear system or wake did not converge, there is likely a serious problem that would lead to an error in the solve, so we will exit here with a fail flag for an optimizer or user
    =#
    if iszero(lu_decomp_flag) || !options.wake_options.converged[1]
        if !options.silence_warnings
            if iszero(lu_decomp_flag)
                @warn "Exiting.  LU decomposition of the LHS matrix for the linear system failed.  Please check your body geometry and ensure that there will be no panels lying directly atop eachother or other similar problematic geometry."
            elseif !options.wake_options.converged[1]
                @warn "Exiting. Wake elliptic grid solve did not converge. Consider a looser convergence tolerance if the geometry looks good."
            end
        end
        #TODO: write a function that returns the same as outs below, but all zeros
        return [],#zero_outputs(),
        (; solve_parameter_tuple..., ivb, airfoils, idmaps, panels, problem_dimensions),
        false
    end

    # - Continue with Analysis - #
    return analyze(
        propulsor,
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        airfoils,
        ivb,
        A_bb_LU,
        panels,
        idmaps,
        problem_dimensions,
        options;
        solve_container_caching=solve_container_caching,
        # post_caching=nothing
    )
end

"""
"""
function analyze(
    propulsor::Propulsor,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
    airfoils,
    ivb,
    A_bb_LU,
    panels,
    idmaps,
    problem_dimensions,
    options=set_options();
    # precomp_container_caching=nothing,
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
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        airfoils,
        A_bb_LU,
        solve_container_caching,
        idmaps,
        options,
    )

    # - Post-Process - #
    # TODO: probably just move this bit inside the post-process to clean up inputs
    # actually probably want to do a second dispatch, one with all the inputs, and one that extracts all the caches and calls the one with all the inputs.
    solve_parameter_tuple = withdraw_solve_parameter_cache(
        options.solve_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )
    (; ivr, ivw, linsys, blade_elements, wakeK) = solve_parameter_tuple

    # NOTE: post processing cache doesn't need to be fancy, since it's just here in the analysis and will always be the same type if derivative checks are turned off.
    # update tests for this function

    # TODO: need to figure out how to set this up to work with the various residual options
#########################################################
##########################     ##########################
#####################     LOOK!    ######################
###########                                   ###########
#####     -----    TODO: YOU ARE HERE     -----     #####
###########                                   ###########
#####################     LOOK!    ######################
##########################     ##########################
#########################################################
    outs = post_process(
        options.solve_options,
        solve_container_caching,
        velocity_states,
        vectorize_velocity_states(
            solve_parameter_tuple.vz_rotor,
            solve_parameter_tuple.vtheta_rotor,
            solve_parameter_tuple.Cm_wake,
        )[2], # TODO: think of a better way to get the state_dims
        propulsor.operating_point,
        propulsor.reference_parameters,
        ivr,
        ivw,
        ivb,
        (; linsys..., A_bb_LU),
        (; blade_elements..., airfoils...),
        wakeK,
        panels.body_vortex_panels,
        panels.rotor_source_panels.influence_length,
        idmaps,
        problem_dimensions;
        write_outputs=options.write_outputs,
        outfile=options.outfile,
        checkoutfileexists=options.checkoutfileexists,
        output_tuple_name=options.output_tuple_name,
        verbose=options.verbose,
    )

    return outs,
    (; ivr, ivw, ivb, linsys, blade_elements, wakeK, idmaps, panels),
    options.solve_options.converged[1]
end

"""
"""
function process(
    solve_options::SolverOptions,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
    airfoils,
    A_bb_LU,
    solve_container_caching,
    idmaps,
    options,
)

    # - Initialize Aero - #
    if options.verbose
        println("\nInitializing Velocities")
    end
    # view the initial conditions out of the inputs cache
    solve_parameter_tuple = withdraw_solve_parameter_cache(
        solve_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )
    (; vz_rotor, vtheta_rotor, Cm_wake, operating_point, linsys, ivr, ivw) =
        solve_parameter_tuple

    # initialize velocities
    # TODO; add some sort of unit test for this function
    initialize_velocities!(
        vz_rotor,
        vtheta_rotor,
        Cm_wake,
        solve_parameter_tuple.operating_point,
        (; solve_parameter_tuple.blade_elements..., airfoils...),
        (; solve_parameter_tuple.linsys..., A_bb_LU),
        solve_parameter_tuple.ivr,
        solve_parameter_tuple.ivw,
        idmaps.body_totnodes,
        idmaps.wake_panel_sheet_be_map,
    )

    # TODO: find a better way to do this so that state_dims is accessible to implicitAD
    _, state_dims = vectorize_velocity_states(
        reshape(vz_rotor, :), reshape(vtheta_rotor, :), reshape(Cm_wake, :)
    )

    # - combine cache and constants - #
    const_cache = (;
        # - Constants - #
        options.verbose,
        options.silence_warnings,
        #nlsolve options
        solve_options,
        # Constant Parameters
        airfoils,
        A_bb_LU,
        idmaps,
        state_dims,
        solve_parameter_cache_dims,
        # Cache(s)
        solve_container_caching...,
    )

    # - Solve with ImplicitAD - #
    if options.verbose
        println("\nSolving Nonlinear System using Newton Method")
    end
    return ImplicitAD.implicit(
        solve_iad, residual!, solve_parameter_cache_vector, const_cache
    )
end

"""
"""
function process(
    solve_options::CSORSolverOptions,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
    airfoils,
    A_bb_LU,
    solve_container_caching,
    idmaps,
    options,
)

    # - Rename for Convenience - #
    (; verbose, solve_options) = options
    # view the initial conditions out of the inputs cache
    solve_parameter_tuple = withdraw_solve_parameter_cache(
        solve_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )
    (; Gamr, sigr, gamw, operating_point, blade_elements, linsys, ivr, ivw, wakeK) =
        solve_parameter_tuple

    # - Initialize States - #
    initialize_strengths!(
        Gamr,
        sigr,
        gamw,
        operating_point,
        (; blade_elements..., airfoils...),
        (; linsys..., A_bb_LU),
        ivr,
        ivw,
        wakeK,
        idmaps.body_totnodes,
        idmaps.wake_nodemap,
        idmaps.wake_endnodeidxs,
        idmaps.wake_panel_sheet_be_map,
        idmaps.wake_node_sheet_be_map,
        idmaps.wake_node_ids_along_casing_wake_interface,
        idmaps.wake_node_ids_along_centerbody_wake_interface,
    )

    sigr .= 0

    state_variables, state_dims = vectorize_strength_states(Gamr, sigr, gamw)

    # - combine cache and constants - #
    const_cache = (;
        # - Constants - #
        options.silence_warnings,
        #CSOR solve options
        solve_options,
        # Constant Parameters
        airfoils,
        A_bb_LU,
        idmaps,
        state_dims,
        # Cache(s)
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        solve_container_caching...,
    )

    #TODO: update this function
    # need to withdraw cache inside this function to have the *_est values, as well as sigr and gamb available.
    return solve!(state_variables, const_cache; verbose=options.verbose)
end

"""
"""
function post_process(
    process_type::SolverOptions,
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
    verbose=false,
)
    if verbose
        println("\nPost-Processing")
    end
    return post_process_iad!(
        process_type,
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
        write_outputs=write_outputs,
        outfile=outfile,
        checkoutfileexists=checkoutfileexists,
        output_tuple_name=output_tuple_name,
    )
end

"""
"""
function post_process(
    process_type::CSORSolverOptions,
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
