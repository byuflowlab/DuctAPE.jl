#---------------------------------#
#              SOLVE              #
#---------------------------------#

"""
"""
function solve_iad(sensitivity_parameters, const_cache)

    # - Extract constants - #
    (;
        # General
        verbose,
        silence_warnings,
        # nlsolve options
        solve_options,
        # Constant Parameters
        airfoils,
        A_bb_LU,
        idmaps,
        solve_parameter_cache_dims,
        # Cache(s)
        solve_container_cache,
        solve_container_cache_dims,
    ) = const_cache

    (;
        nlsolve_method,
        nlsolve_linesearch_method,
        nlsolve_linesearch_kwargs,
        nlsolve_ftol,
        nlsolve_iteration_limit,
        converged,
    ) = solve_options

    # - Extract Initial Guess Vector for State Variables - #
    # TODO: rename this to vectorize velocities since the state_dims are named tuples of the velocity names
    initial_guess, state_dims = vectorize_velocity_states(
        reshape(
            @view(sensitivity_parameters[solve_parameter_cache_dims.vz_rotor.index]),
            solve_parameter_cache_dims.vz_rotor.shape,
        ),
        reshape(
            @view(sensitivity_parameters[solve_parameter_cache_dims.vtheta_rotor.index]),
            solve_parameter_cache_dims.vtheta_rotor.shape,
        ),
        reshape(
            @view(sensitivity_parameters[solve_parameter_cache_dims.Cm_wake.index]),
            solve_parameter_cache_dims.Cm_wake.shape,
        ),
    )

    # - Wrap residual - #
    if verbose
        println("  " * "Wrapping Residual and Jacobian")
    end
    function rwrap!(resid, state_variables)
        return residual!(
            resid,
            state_variables,
            sensitivity_parameters,
            (;
                solve_options, # for dispatch
                airfoils,                   # inner and outer airfoil objects along blades
                A_bb_LU,                    # LU decomposition of linear system LHS
                idmaps,                     # book keeping items
                state_dims,                 # dimensions for state variable extraction
                solve_parameter_cache_dims, # dimensions for shaping sensitivity parameters
                solve_container_cache,      # cache for solve_containers used in solve
                solve_container_cache_dims, # dimensions for shaping the view of the solve cache
            ),
        )
    end

    # - Wrap Jacobian - #
    # configure jacobian
    jconfig = ForwardDiff.JacobianConfig(
        rwrap!,
        initial_guess, # object of size of residual vector
        initial_guess, # object of size of state_variable vector
        ForwardDiff.Chunk{12}(), #TODO: set the chunk size as an option and make sure the solver cache chunk and this chunk size match. PreallocationTools chooses poorly if not given a chunk size
    )
    # get allocated array for jacobian
    JR = DiffResults.JacobianResult(
        initial_guess, # object of size of residual vector
        initial_guess, # object of size of state_variable vector
    )
    # store everything in a cache that can be accessed inside the solver
    jcache = (; rwrap!, JR, config=jconfig, resid=zeros(size(initial_guess)))

    # wrap the jacobian
    function jwrap!(J, state_variables)
        ForwardDiff.jacobian!(
            jcache.JR, jcache.rwrap!, jcache.resid, state_variables, jcache.config
        )
        J .= DiffResults.jacobian(jcache.JR)
        return J
    end

    # build the OnceDifferentiable object
    df = NLsolve.OnceDifferentiable(
        rwrap!, jwrap!, initial_guess, similar(initial_guess) .= 0
    )

    # - SOLVE - #
    if verbose
        println("  " * "Newton Solve Trace:")
    end
    result = NLsolve.nlsolve(
        df,
        initial_guess; # initial states guess
        method=nlsolve_method,
        # autodiff=nlsolve_autodiff,
        linesearch=nlsolve_linesearch_method(nlsolve_linesearch_kwargs...),
        linsolve=(x, A, b) -> x .= ImplicitAD.implicit_linear(A, b),
        ftol=nlsolve_ftol,
        iterations=nlsolve_iteration_limit,
        show_trace=verbose,
    )

    # update convergence flag
    # note: need to do this complicated check to ensure that we're only checking |f(x)|<ftol rather than claiming convergences when the step size is zero but it's really just stuck.
    _, converged[1] = NLsolve.assess_convergence(NLsolve.value(df), nlsolve_ftol)

    return result.zero
end

#---------------------------------#
#             RESIDUALS           #
#---------------------------------#

"""
This is the residual that gets passed into nlsolve, and is just for converging Gamr and gamw.
"""
function residual!(resid, state_variables, sensitivity_parameters, constants)

    # - Extract constants - #
    (;
        # solve options for dispatch
        solve_options,
        # comp,
        airfoils,                   # airfoils
        A_bb_LU,                    # linear system left hand side LU decomposition
        idmaps,                     # book keeping items
        state_dims,                 # dimensions for state variable extraction
        solve_parameter_cache_dims, # dimensions for shaping the view of the parameter cache
        solve_container_cache,      # cache for solve_containers used in solve
        solve_container_cache_dims, # dimensions for shaping the view of the solve cache
    ) = constants

    # - Separate out the state variables - #
    # TODO: test this function
    vz_rotor, vtheta_rotor, Cm_wake = extract_state_variables(
        solve_options, state_variables, state_dims
    )

    # separate out sensitivity_parameters here as well
    solve_parameter_tuple = withdraw_solve_parameter_cache(
        solve_options, sensitivity_parameters, solve_parameter_cache_dims
    )

    (;
        operating_point, # freestream, Omega, etc.
        ivr,             # induced velocities on rotor panels
        ivw,             # induced velocities on wake panels
        linsys,          # includes AIC's for linear system
        blade_elements,  # includes blade element geometry
        wakeK,           # geometric "constants" for wake node strength calculation
    ) = solve_parameter_tuple

    # println(
    #     "comp: ",
    #     all(compare_arrays(Cm_wake, comp; verbose=true, parent="Cm_wake", tab="    ")),
    # )

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
        solve_options, solve_container_cache_vector, solve_container_cache_dims
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
    solve_containers.vz_est .-= vz_rotor
    solve_containers.vtheta_est .-= vtheta_rotor
    solve_containers.Cm_est .-= Cm_wake

    resid[state_dims.vz_rotor.index] .= reshape(solve_containers.vz_est, :)
    resid[state_dims.vtheta_rotor.index] .= reshape(solve_containers.vtheta_est, :)
    resid[state_dims.Cm_wake.index] .= solve_containers.Cm_est

    ##NOTE: can't use a smaller residual with NLsolve
    #resid[1] = solve_containers.vz_est[findmax(abs.(solve_containers.vz_est))[2]]
    #resid[2] = solve_containers.vtheta_est[findmax(abs.(solve_containers.vtheta_est))[2]]
    #resid[3] = solve_containers.Cm_est[findmax(abs.(solve_containers.Cm_est))[2]]

    # if eltype(Cm_wake) != Float64
    #     comp .= (p->p.value).(Cm_wake)
    # else
    #     comp .= Cm_wake
    # end

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
    # println("Cz type: ", eltype(solve_containers.Cz_rotor) != Float64 ? "Dual" : "Float")
    # println("vz type: ", eltype(vz_rotor) != Float64 ? "Dual" : "Float")
    # println("Vinf type: ", eltype(operating_point.Vinf) != Float64 ? "Dual" : "Float")
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
