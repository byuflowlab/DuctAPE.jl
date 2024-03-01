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
    silence_warnings=true,
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
    # - Geometry Re-interpolation and generation options - #
    finterp=fm.akima,
    max_wake_relax_iter=100,
    wake_relax_tol=1e-9,
    itcpshift=0.05,
    axistol=1e-15,
    tegaptol=1e1 * eps(),
)
    return (;
        # - General Parameters - #
        verbose,
        silence_warnings,
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
        # - Post-processing Options - #
        write_outputs,
        outfile,
        checkoutfileexists,
        tuple_name,
        # - Geometry Re-interpolation and generation options - #
        finterp,
        max_wake_relax_iter,
        wake_relax_tol,
        itcpshift,
        axistol,
        tegaptol,
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
    velocity_states = implicit(solve_iad, solve_iad_res!, x, const_cache)

    # - Post-Process - #
    # do the rest of the post-processing
    # TODO: write this function
    # NOTE: post processing cache doesn't need to be fancy, since it's just here in the analysis and will always be the same type if derivative checks are turned off.
    outs = post_process_iad(velocity_states, fx(x), cache_and_dims)

    return outs
end

"""
TODO: move to analysis.jl
"""
function analyze(
    propulsor,
    constants=define_constants();
    precomp_container_caching=generate_precomp_container_cache(repanel(propulsor)),
    solve_parameter_caching=generate_solve_parameters_cache(repanel(propulsor)),
    solve_container_caching=generate_solve_containter_cache(repanel(propulsor)),
    post_caching=generate_post_cache(repanel(propulsor)),
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
    inputs = vectorize!(propulsor)

    # - combine cache and constants - #
    const_cache = (;
        # - Constants - #
        #nlsolve options
        constants.nlsolve_method,
        constants.nlsolve_autodiff,
        constants.nlsolve_linesearch,
        constants.nlsolve_ftol,
        constants.nlsolve_iteration_limit,
        constants.verbose,
        constants.silence_warnings,
        constants.converged,
        #geometry options
        constants.finterp,
        constants.max_wake_relax_iter,
        constants.wake_relax_tol,
        constants.itcpshift,
        constants.axistol,
        constants.tegaptol,
        # TODO: write generate_caches function (need to go through and figure out what all goes in the cache)
        # Caches
        precomp_container_caching...,
        # precomp_container_cache_dims,
        solve_parameter_caching...,
        # solve_parameter_cache_dims,
        solve_container_caching...,
        # solve_container_cache_dims,
        fx=x -> (propulsor),
    )

    # - Initialize Aero and Converge Gamr and gamw - #
    velocity_states = implicit(solve_iad, solve_iad_res!, inputs, const_cache)

    # - Post-Process - #

    # do the rest of the post-processing
    # TODO: write this function
    # NOTE: post processing cache doesn't need to be fancy, since it's just here in the analysis and will always be the same type if derivative checks are turned off.
    outs = post_process_iad(velocity_states, propulsor, post_caching)

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
        #general
        verbose,
        silence_warnings,
        #inputs
        fx,
        #solve
        nlsolve_method,
        nlsolve_autodiff,
        nlsolve_linesearch,
        nlsolve_ftol,
        nlsolve_iteration_limit,
        converged,
        #geometry
        finterp,
        max_wake_relax_iter,
        wake_relax_tol,
        itcpshift,
        axistol,
        tegaptol,
        #caches
        precomp_container_cache,
        precomp_container_cache_dims,
        solve_parameter_cache,
        solve_parameter_cache_dims,
        solve_container_cache,
        solve_container_cache_dims,
    ) = const_cache

    # - Get propulsor back out using fx(x) in constants - #
    # TODO: figure out how to make sure these are floats in this context, probably need to look at how implicitAD (or eduardo in flowpanel) strips things to floats and apply something similar when putting together the outputs of the fx function (dispatch on type of inputs maybe?)
    # propulsor includes: duct_coordinates, centerbody_coordinates, rotorstator_parameters, paneling_constants, operating_point, and reference_parameters
    (; propulsor) = fx(inputs)

    # - Extract precomp_container_cache - #
    precomp_container_cache_vec = @views pat.get_tmp(precomp_container_cache, inputs)
    # TODO: test this function
    precomp_containers = withdraw_precomp_container_cache(precomp_container_cache, precomp_container_cache_dims)

    # - Extract solve_parameter_cache - #
    # TODO; can these just be floats?  If so, want a different setup than using PreallocationTools.  Want something closer to the fx function for inputs. but return the parameters cache, also want to zero it out.
    solve_parameter_cache_vec = @views pat.get_tmp(solve_parameter_cache, inputs)
    # TODO: test this function
    solve_parameter_containers = withdraw_solve_parameter_cache(
        solve_parameter_cache, solve_parameter_cache_dims
    )

    # - Do precomputations - #
    # TODO: test this function
    precompute_solve_parameters_iad!(
        solve_parameter_containers.ivr,
        solve_parameter_containers.ivw,
        solve_parameter_containers.ivb,
        solve_parameter_containers.linsys,
        solve_parameter_containers.blade_elements,
        solve_parameter_containers.idmaps,
        propulsor,
        precomp_containers;
        max_wake_relax_iter=max_wake_relax_iter,
        wake_relax_tol=wake_relax_tol,
        itcpshift=itcpshift,
        axistol=axistol,
        tegaptol=tegaptol,
        finterp=finterp,
        silence_warnings=silence_warnings,
    )

    # - Initialize Aero - #
    # TODO; test this function
    Vz_rotor, Vtheta_rotor, Cm_wake = initialize_velocities(
        op,
        blade_elements,
        linsys,
        ivr,
        ivw,
        solve_parameter_containers.idmaps.body_totnodes,
        solve_parameter_containers.idmaps.rotorwakenodeid,
    )

    # - Wrap residual - #
    function rwrap!(r, state_variables)
        #TODO: Test this function
        return nls_res!(
            r,
            state_variables,
            (;
                op,             # includes freestream and Omega
                ivr,            # induced velocities on rotor panels
                ivw,            # induced velocities on wake panels
                linsys,         # includes AIC's for linear system
                blade_elements, # includes blade element geometry
                idmaps,         # book keeping items
                solve_container_cache,      # cache for containers used in solve
                solve_container_cache_dims, # dimensions for shaping the view of the solve cache
                solve_container_cache_dims.state_dims, # dimensions for state variable extraction
                # solve_parameter_containers: precomputed parameters used in state estimation
                solve_parameter_containers.ivr,
                solve_parameter_containers.ivw,
                solve_parameter_containers.linsys,
                solve_parameter_containers.blade_elements,
                solve_parameter_containers.idmaps,
            ),
        )
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
function solve_iad_res!(r, state_variables, inputs, const_cache) end

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
        solve_container_cache,      # cache for containers used in solve
        solve_container_cache_dims, # dimensions for shaping the view of the solve cache
        ivr,              # induced unit velocities on rotors
        ivw,              # induced unit velocities on wakes
        linsys,           # Linear System components
        blade_elements,   # Blade element geometry and airfoils
        idmaps,           # index maps for accessing some tricky things
    ) = parameters

    # - Separate out the state variables - #
    # TODO: test this function
    Vz_rotor, Vtheta_rotor, Cm_wake = extract_state_vars(state_variables, state_dims)

    # - Extract and Reset Cache - #
    # get cache vector of correct types
    solve_container_cache_vec = @views pat.get_tmp(solve_container_cache, state_variables)
    # TODO: test this function
    solve_containers = withdraw_solve_container_cache(
        solve_container_cache_vec, solve_container_cache_dims
    )
    # zero out contents of solve_containers to avoid any potential contamination issues
    # TODO: test this function
    reset_solve_containers!(solve_containers) #note: also zeros out state estimates

    # - Estimate New States - #
    # TODO: test this function
    estimate_states!(
        solve_containers.Vz_est,
        solve_containers.Vtheta_est,
        solve_containers.Cm_est,
        Vz_rotor,
        Vtheta_rotor,
        Cm_wake,
        solve_containers,
        op,
        ivr,
        ivw,
        linsys,
        blade_elements,
        idmaps,
    )

    # - Get final Residual Values - #
    r .= [
        reshape(solve_containers.Vz_est .- Vz_rotor, length(Vz_rotor))
        reshape(solve_containers.Vtheta_est .- Vtheta_rotor, length(Vtheta_rotor))
        solve_containers.Cm_est .- Cm_wake
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
    idmaps;
    verbose=false,
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
        verbose=verbose,
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
        idmaps.rotorwakenodeid,
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
        idmaps.rotorwakenodeid,
        idmaps.ductwakeinterfacenodeid,
        idmaps.cbwakeinterfacenodeid;
        post=false,
    )

    # - Solve Linear System for Body Strengths - #
    calculate_body_vortex_strengths!(
        containers.gamb,
        linsys.A_bb_LU;
        linsys.b_bf,
        containers.gamw,
        linsys.A_bw,
        linsys.A_pw,
        containers.sigr,
        linsys.A_br,
        linsys.A_pr,
        linsys.A_bb,
    )

    # - Calcuate Vz_est and Vtheta_est- #
    # TODO: test this function
    calculate_induced_velocities_on_rotors!(
        Vz_est,
        Vtheta_est,
        containers.Gamr,
        containers.gamw,
        containers.sigr,
        containers.gamb[1:(idmaps.body_totnodes)],
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
        containers.gamb[1:(idmaps.body_totnodes)],
        ivw,
        op.Vinf,
    )

    # return estimated states
    return Vz_est, Vtheta_est, Cm_est
end
