######################################################################
#                                                                    #
#                        FUNCTION SET HEADER                         #
#                                                                    #
######################################################################

# """
# notional objective function for an optimization
# """
# function objective!(constraints, design_vars, optimization_parameters; cache=make_a_cache)
#     parameters = taylors_fancy_package_function(optimization_parameters, design_vars)
#     propulsor = f(parameters)
#     aero_outs = analyze(propulsor; cache=cache)
#     objective = fun(aero_outs)
#     func!(constraints, aero_outs)
#     return objective
# end

"""
TODO: move to analysis.jl
"""
function define_constants(;
    # - General Parameters - #
    verbose=false,
    # - NLSolve Parameters - #
    # nlsolve parameters
    nlsolve_method=:newton,
    nlsolve_autodiff=:forward,
    nlsolve_ftol=1e-8, #1e-8 is nlsolve default
    nlsolve_iteration_limit=1e3, #1000 is nlsolve default
    show_trace=false,
    # line search parameters
    dolinesearch=true,
    maximum_linesearch_step_size=1e6,
    # - Post-processing Options - #
    write_outputs=false,
    outfile="outputs.jl",
    checkoutfileexists=false,
    tuple_name="outs",
)
    return (;
        # - General Parameters - #
        verbose,
        # - NLSolve Parameters - #
        # nlsolve parameters
        nlsolve_method,
        nlsolve_autodiff,
        nlsolve_ftol,
        nlsolve_iteration_limit,
        show_trace,
        # line search parameters
        nlsolve_linesearch=if dolinesearch
            BackTracking(; maxstep=maximum_linesearch_step_size)
        else
            LineSearches.Static()
        end,
    )
end

"""
Version of `analyze_propeller` designed for sensitivity analysis.  `x` is an input vector
and `fx` is a function which returns the standard input arguments to `analyze_propeller`.
"""
function analyze(x, fx=x -> x, constants=define_constants())

    # - combine cache and constants - #
    const_cache = (;
        constants...,
        # parameters,
        # TODO: write generate_cache function (need to go through and figure out what all goes in the cache) (also need to figure out the preallocationtools package.)
        cache_dims=cache_and_dims.cache_dims,
        cache=cache_and_dims.cache,
        fx,
    )

    # - Initialize Aero and Converge Gamr and gamw - #
    Gamr_gamw_sigr = implicit(solve_iad, solve_iad_res!, x, const_cache)

    # - Post-Process - #
    # de-vectorize Gamr and gamw
    # TODO: write this function
    Gamr, gamw, sigr = extract_state_vars(Gamr_gamw_sigr, const_cache.cache_dims)

    # do the rest of the post-processing
    # TODO: write this function
    outs = post_process_iad(Gamr, gamw, sigr, x, cache_and_dims)

    return outs
end

"""
TODO: move to analysis.jl
"""
function analyze(
    propulsor,
    constants=define_constants();
    cache_and_dims=generate_cache(repanel(propulsor)),
)

    # - Check that propulsor has the required fields - #
    @assert all(
        haskey.(
            Ref(propulsor),
            [
                :duct_coordinates,
                :centerbody_coordinates,
                :rotorstator_parameters,
                :paneling_constants,
                :freestream,
                :reference_parameters,
            ],
        ),
    )

    # - Vectorize propulsor - #
    # inputs, parameters = vectorize!(propulsor)
    inputs, _ = vectorize!(propulsor)

    # - combine cache and constants - #
    const_cache = (;
        constants...,
        # parameters,
        # TODO: write generate_cache function (need to go through and figure out what all goes in the cache)
        dimensions=cache_and_dims.cache_dims,
        cache=cache_and_dims.cache,
        fx=x -> (propulsor),
    )

    # - Initialize Aero and Converge Gamr and gamw - #
    velocity_states = implicit(solve_iad, solve_iad_res!, inputs, const_cache)

    # - Post-Process - #
    # reshape velocities
    Vz_rotor, Vtheta_rotor, Cm_wake = extract_state_vars(velocity_states, state_dims)

    # do the rest of the post-processing
    # TODO: write this function
    outs = post_process_iad(Vz_rotor, Vtheta_rotor, Cm_wake, propulsor, cache_and_dims)

    return outs
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
        fx,
        nlsolve_method,
        nlsolve_autodiff,
        nlsolve_linesearch,
        nlsolve_ftol,
        nlsolve_iteration_limit,
        verbose,
        converged,
        # TODO: work backwards from estimate_states and figure out what needs to get passed in here as far as cache's, dimensions, and parameters go.
        cache,
        parameters,
        dimensions,
    ) = const_cache

    # - Get propulsor back out using fx(x) in constants - #
    # TODO: figure out how to make sure these are floats in this context, probably need to look at how implicitAD strips things to floats and apply something similar when putting together the outputs of the fx function
    (; propulsor) = fx(inputs)

    # - Do precomputations - #
    # TODO: re-think how to set up this function and outputs.
    bp, rp, op = precompute_parameters_iad!(propulsor)

    # - Initialize Aero - #
    # TODO: write this function
    # Note: you can solve the isolated body system and use that for the initial guess here because it's just an initial guess, it doesn't matter what it is relative to the differentiated variables.
    # TODO: re-think the inputs to this function, don't just use mega tuples.
    Vz_rotor, Vtheta_rotor, Cm_wake = initialize_velocities(parameters, cache)

    # - Wrap residual - #
    function rwrap!(r, state_variables)
        # Note: the inputs going into ImplicitAD don't need to get passed in to NLsolve, right?
        # return nls_res!(r, state_variables, inputs, constants)
        return nls_res!(r, state_variables, (; parameters, dimensions, cache))
    end

    result = NLsolve.nlsolve(
        rwrap!,
        [Vz_rotor; Vtheta_rotor; Cm_wake]; # initial states guess
        method=nlsolve_method,
        autodiff=nlsolve_autodiff,
        linesearch=nlsolve_linesearch,
        linsolve=(x, A, b) -> x .= ImplicitAD.implicit_linear(A, b), #used in newton method, unused otherwise so fine to keep it here.
        ftol=nlsolve_ftol,
        iterations=nlsolve_iteration_limit,
        show_trace=verbose,
    )

    # update convergence flag
    converged[1] = NLsolve.converged(result)

    return result.zero
end

#---------------------------------#
#             RESIDUALS           #
#---------------------------------#

"""
This is the residual that gets passed into ImplicitAD, and includes all the precomputation stuff.
basically just copy and paste the solve function but replace the nlsolve with a simple call the the nlsolve residual function.
 - Do all the initialization and setup stuff
 - call the nlsolve residual function
 - return the residual and various outputs if you want to return them

"""
function solve_iad_res!(r, state_variables, inputs, const_cache)

    # - Extract constants - #
    (;
        fx,
        nlsolve_method,
        nlsolve_autodiff,
        nlsolve_linesearch,
        nlsolve_ftol,
        nlsolve_iteration_limit,
        verbose,
        parameters,
        dimensions,
        cache,
    ) = const_cache

    # - Get propulsor back out using fx(x) in constants - #
    (; propulsor) = fx(inputs)

    # - Do precomputations - #
    populate_cache!(cache, propulsor)

    # - Initialize Aero - #
    Gamr, gamw, sigr = extract_state_vars(state_variables)

    # - Call nlsolve residual - #
    return nls_res!(r, [Gamr; gamw; sigr], (; parameters, dimensions, cache))
end

"""
This is the residual that gets passed into nlsolve, and is just for converging Gamr and gamw.
"""
function nls_res!(r, state_variables, parameters)

    # - Extract Parameters - #
    (;
        op,               # includes freestream and Omega
        ivr,              # induced velocities on rotor panels
        ivw,              # induced velocities on wake panels
        linsys,           # includes AIC's for linear system
        blade_elements,   # includes blade element geometry
        idmaps,           # book keeping items
        state_dims,       # dimensions for state variable extraction
        solve_cache,      # cache for containers used in solve
        solve_cache_dims, # dimensions for shaping the view of the solve cache
    ) = parameters

    # - Separate out the state variables - #
    # TODO: test this function
    Vz_rotor, Vtheta_rotor, Cm_wake = extract_state_vars(state_variables, state_dims)

    # - Extract and Reset Cache - #
    # TODO: test this process.
    # get cache vector of correct types
    cache_vec = @views pat.get_tmp(solve_cache, state_variables)
    # view cache vector and get named tuple of containers used for intermediate variables in solve
    # TODO: write this function
    containers = withdraw_solve_cache(cache_vec, solve_cache_dims)
    # zero out contents of containers to avoid any potential contamination issues
    # TODO: test this function
    reset_containers!(containers) #note: also zeros out state estimates

    # - Estimate New States - #
    # TODO: test this function
    estimate_states!(
        containers.Vz_est,
        containers.Vtheta_est,
        containers.Cm_est,
        Vz_rotor,
        Vtheta_rotor,
        Cm_wake,
        containers,
        op,
        ivr,
        ivw,
        linsys,
        blade_elements,
        idmaps,
    )

    # - Get final Residual Values - #
    r .= [
        reshape(containers.Vz_est .- Vz_rotor, length(Vz_rotor))
        reshape(containers.Vtheta_est .- Vtheta_rotor, length(Vtheta_rotor))
        containers.Cm_est .- Cm_wake
    ]

    return r
end

"""
"""
function estimate_states!(
    Vz_est,
    Vtheta_est,
    Cm_est,
    Vz_rotor,
    Vtheta_rotor,
    Cm_wake,
    containers,
    op,
    ivr,
    ivw,
    linsys,
    blade_elements,
    idmaps,
)

    # - Get Absolute Rotor Velocities - #
    # TODO: test this function
    reframe_rotor_velocities!(
        containers.Cz_rotor,
        containers.Ctheta_rotor,
        containers.Cmag_rotor,
        Vz_rotor, # state var
        Vtheta_rotor, # state var
        op.Vinf,
        Omega,
        rotor_panel_centers,
    )

    # - Calculate Blade Element Values - #
    # TODO: test this function
    calculate_blade_element_coefficients!(
        containers.cl,
        containers.cd,
        blade_elements,
        containers.Cz_rotor,
        containers.Ctheta_rotor,
        containers.Cmag_rotor,
        op;
        post=false,
        verbose=false,
    )

    # - Calculate Blade Element Circulation - #
    # in-place solve for Gamr,
    # TODO: test this function
    calculate_rotor_circulation_strengths!(
        containers.Gamr, containers.Wmag_rotor, blade_elements.chords, containers.cl
    )

    # - Calculate Rotor Panel Strengths - #
    # in-place solver for sigr,
    # TODO: test this function
    calculate_rotor_source_strengths!(
        containers.sigr,
        containers.Cmag_rotor,
        blade_elements.chords,
        blade_elements.B,
        containers.cd,
    )

    # - Get Average Wake Velocities - #
    # TODO: test this function
    average_wake_velocities!(
        containers.Cm_avg,
        Cm_wake,
        idmaps.wake_nodemap,
        idmaps.wake_endnodeidxs,
        idmaps.rotorwakeid,
    )

    # - Calculate Wake Panel Strengths - #
    # in-place solve for gamw,
    # TODO: test this function
    calculate_wake_vortex_strengths!(
        containers.gamw,
        containers.Gamr,
        containers.Cm_avg,
        blade_elements.B,
        op.Omega,
        wakeK,
        idmaps.rotorwakeid,
        idmaps.ductwakeinterfacenodeid,
        idmaps.hubwakeinterfacenodeid;
        post=false,
    )

    # - Solve Linear System for Body Strengths - #
    calculate_body_vortex_strengths!(
        containers.gamb,
        linsys.A_bb,
        linsys.b_bf,
        containers.gamw,
        linsys.A_bw,
        linsys.A_pw,
        containers.sigr,
        linsys.A_br,
        linsys.A_pr,
        linsys.A_bb_LU;
    )

    # - Calcuate Vz_est and Vtheta_est- #
    # TODO: test this function
    calculate_induced_velocities_on_rotors!(
        Vz_est,
        Vtheta_est,
        containers.Gamr,
        containers.gamw,
        containers.sigr,
        containers.gamb[1:(idmaps.totnodes_body)],
        ivr,
        blade_elements.B,
        blade_elements.rotor_panel_center,
    )

    # - Calculate Velocities on Wake Panels - #
    # TODO: test this function
    calculate_wake_velocities!(
        Cm_est,
        containers.vz_wake,
        containers.vr_wake,
        containers.gamw,
        containers.sigr,
        containers.gamb[1:(idmaps.totnodes_body)],
        ivw,
        op.Vinf,
    )

    # return estimated states
    return Vz_est, Vtheta_est, Cm_est
end
