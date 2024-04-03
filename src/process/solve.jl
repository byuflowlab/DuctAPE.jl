#---------------------------------#
#              SOLVE              #
#---------------------------------#

# TODO: add dispatches for SpeedMapping and FixedPointAcceleration

"""
"""
function solve(sensitivity_parameters, const_cache)
    return solve(const_cache.solver_options, sensitivity_parameters, const_cache)
end

"""
"""
function solve(solver_options::NonlinearSolveOptions, sensitivity_parameters, const_cache)

    # - Extract constants - #
    (;
        # General
        verbose,
        silence_warnings,
        # nlsolve options
        solver_options,
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
        nlsolve_algorithm,
        # # General Controls
        nlsolve_alias_initial_guess,
        # # Iteration Controls
        nlsolve_abstol,
        nlsolve_maxiters,
        converged,
    ) = solver_options

    # - Extract Initial Guess Vector for State Variables - #
    initial_guess = extract_initial_guess(
        solver_options, sensitivity_parameters, solve_parameter_cache_dims.state_dims
    )

    # - Wrap residual - #
    if verbose
        println("  " * "Wrapping Residual and Jacobian")
    end

    function rwrap!(resid, state_variables, p)
        return system_residual!(
            resid, state_variables, p.sensitivity_parameters, p.constants
        )
    end

    # build problem object
    prob = NonlinearSolve.NonlinearProblem(
        rwrap!,
        initial_guess,
        (;
            sensitivity_parameters,
            constants=(;
                solver_options, # for dispatch
                airfoils,                   # inner and outer airfoil objects along blades
                A_bb_LU,                    # LU decomposition of linear system LHS
                idmaps,                     # book keeping items
                solve_parameter_cache_dims, # dimensions for shaping sensitivity parameters
                solve_container_cache,      # cache for solve_containers used in solve
                solve_container_cache_dims, # dimensions for shaping the view of the solve cache
            ),
        ),
    )

    # - SOLVE - #
    if verbose
        println("  " * "Nonlinear Solve Trace:")
    end
    sol = NonlinearSolve.solve(
        prob, # problem
        nlsolve_algorithm();
        abstol=nlsolve_abstol,
        maxiters=nlsolve_maxiters,
    )

    # update convergence flag
    converged[1] = SciMLBase.successful_retcode(sol)

    return sol.u
end

function solve(solver_options::NLsolveOptions, sensitivity_parameters, const_cache)

    # - Extract constants - #
    (;
        # General
        verbose,
        silence_warnings,
        # nlsolve options
        solver_options,
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
    ) = solver_options

    # - Extract Initial Guess Vector for State Variables - #
    initial_guess = extract_initial_guess(
        solver_options, sensitivity_parameters, solve_parameter_cache_dims.state_dims
    )

    # - Wrap residual - #
    if verbose
        println("  " * "Wrapping Residual and Jacobian")
    end
    function rwrap!(resid, state_variables)
        return system_residual!(
            resid,
            state_variables,
            sensitivity_parameters,
            (;
                solver_options, # for dispatch
                airfoils,                   # inner and outer airfoil objects along blades
                A_bb_LU,                    # LU decomposition of linear system LHS
                idmaps,                     # book keeping items
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

function solve(solver_options::CSORSolverOptions, sensitivity_parameters, const_cache)

    # - Extract constants - #
    (;
        # General
        verbose,
        silence_warnings,
        # nlsolve options
        solver_options,
        # Constant Parameters
        airfoils,
        A_bb_LU,
        idmaps,
        solve_parameter_cache_dims,
        # Cache(s)
        solve_container_cache,
        solve_container_cache_dims,
    ) = const_cache

    # - Extract Initial Guess Vector for State Variables - #
    state_variables = extract_initial_guess(
        solver_options, sensitivity_parameters, solve_parameter_cache_dims.state_dims
    )

    # - Separate out the state variables - #
    Gamr, sigr, gamw = extract_state_variables(
        solver_options, state_variables, solve_parameter_cache_dims.state_dims
    )

    # - Extract and Reset Cache - #
    # get cache vector of correct types
    solve_container_cache_vector = @views PreallocationTools.get_tmp(
        solve_container_cache, eltype(state_variables)(1.0)
    )
    # reset cache
    solve_container_cache_vector .= 0
    solve_containers = withdraw_solve_container_cache(
        solver_options, solve_container_cache_vector, solve_container_cache_dims
    )

    # - initialize convergence criteria - #
    @. solve_containers.deltaG_prev = solve_containers.Gamr_est - Gamr
    @. solve_containers.deltaG = 0.0
    @. solve_containers.deltag_prev = solve_containers.gamw_est - gamw
    @. solve_containers.deltag = 0.0

    TF = eltype(state_variables)

    resid = MVector{2,TF}(999 * ones(TF, 2))
    conv = solver_options.converged
    iter = 0

    # - SOLVE - #
    if verbose
        println("  " * "CSOR Solve Trace:")
    end

    # loop until converged or max iterations are reached
    while !conv[] && iter <= solver_options.maxiter
        # update iteration number
        iter += 1
        if verbose
            println("Iteration $(iter):")
        end

        # Call residual
        CSOR_residual!(
            resid,
            state_variables,
            sensitivity_parameters,
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
            ),
        )

        # Check Convergence
        check_CSOR_convergence!(
            conv,
            resid;
            solver_options.f_circ,
            solver_options.f_dgamw,
            convergence_type=typeof(solver_options.convergence_type),
            verbose=verbose,
        )
    end

    return state_variables
end

