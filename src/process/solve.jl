"""
    solve(sensitivity_parameters, const_cache; initial_guess=nothing)

A compact dispatch of `solve` that automatically dispatches based on the solver_options contained in const_cache.
"""
function solve(sensitivity_parameters, const_cache; initial_guess=nothing)
    return solve(
        const_cache.solver_options,
        sensitivity_parameters,
        const_cache;
        initial_guess=initial_guess,
    )
end

#---------------------------------#
#       FIXED POINT SOLVERS       #
#---------------------------------#

"""
    solve(
        solver_options::SolverOptionsType,
        sensitivity_parameters,
        const_cache;
        initial_guess=nothing,
    )

Converge the residual, solving for the state variables that do so.

# Arguments
- `solver_options::SolverOptionsType` : SolverOptionsType used for dispatch
- `sensitivity_parameters::Vector{Float}` : Sensitivity parameters for solve (parameters passed in through ImplicitAD)
- `const_cache::NamedTuple` : A named tuple containing constants and caching helpers.

# Keyword Arguments
- `initial_guess=nothing::Vector{Float}` : An optional manually provided initial guess (contained in the sensitivity parameters anyway).

# Returns
- `converged_states::Vector{Float}` : the states for which the residual has converged.
"""
function solve(
    solver_options::CSORSolverOptions,
    sensitivity_parameters,
    const_cache;
    initial_guess=nothing,
)

    # - Extract constants - #
    (;
        # General
        verbose,
        silence_warnings,
        multipoint_index,
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
    if isnothing(initial_guess)
        state_variables = extract_initial_guess(
            solver_options, sensitivity_parameters, solve_parameter_cache_dims.state_dims
        )
    else
        state_variables = initial_guess
    end

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
    conv = @view(solver_options.converged[multipoint_index[]])
    iter = @view(solver_options.iterations[multipoint_index[]]) .= 0

    # - SOLVE - #
    if verbose
        println("  " * "CSOR Solve Trace:")
    end

    # loop until converged or max iterations are reached
    while !conv[] && iter[] <= solver_options.iteration_limit
        # update iteration number
        iter[] += 1
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

function solve(
    solver_options::SpeedMappingOptions,
    sensitivity_parameters,
    const_cache;
    initial_guess=nothing,
)

    # - Extract constants - #
    (;
        # General
        verbose,
        silence_warnings,
        multipoint_index,
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
        orders,
        sig_min,
        stabilize,
        check_obj,
        atol,
        iteration_limit,
        time_limit,
        lower,
        upper,
        buffer,
        Lp,
        converged,
    ) = solver_options

    # - Extract Initial Guess Vector for State Variables - #
    if isnothing(initial_guess)
        initial_guess = extract_initial_guess(
            solver_options, sensitivity_parameters, solve_parameter_cache_dims.state_dims
        )
    end

    # - Wrap residual - #
    if verbose
        println("  " * "Wrapping Residual")
    end

    function rwrap!(state_estimates, state_variables)
        return system_residual!(
            state_estimates,
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

    # - SOLVE - #
    if verbose
        println("  " * "SpeedMapping Solve:")
    end
    sol = SpeedMapping.speedmapping(
        initial_guess;
        (m!)=rwrap!,
        orders=orders,
        σ_min=sig_min,
        stabilize=stabilize,
        check_obj=check_obj,
        tol=atol,
        maps_limit=iteration_limit,
        time_limit=time_limit,
        lower=lower,
        upper=upper,
        buffer=buffer,
        Lp=Lp,
    )

    # update convergence flag
    converged[multipoint_index[]] = sol.converged

    return sol.minimizer
end

function solve(
    solver_options::FixedPointOptions,
    sensitivity_parameters,
    const_cache;
    initial_guess=nothing,
)

    # - Extract constants - #
    (;
        # General
        verbose,
        silence_warnings,
        multipoint_index,
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
    if isnothing(initial_guess)
        initial_guess = extract_initial_guess(
            solver_options, sensitivity_parameters, solve_parameter_cache_dims.state_dims
        )
    end

    # - Wrap residual - #
    if verbose
        println("  " * "Wrapping Residual")
    end

    function rwrap!(state_estimates, state_variables)
        return system_residual!(
            state_estimates,
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

    # - SOLVE - #
    if verbose
        println("  " * "FixedPoint Solve:")
    end
    sol = FixedPoint.afps!(
        rwrap!,
        initial_guess;
        iters=solver_options.iteration_limit,
        vel=solver_options.vel,
        ep=solver_options.ep,
        tol=solver_options.atol,
    )

    # update convergence flag
    solver_options.converged[multipoint_index[]] = sol.error <= solver_options.atol

    return sol.x
end

#---------------------------------#
#      QUASI-NEWTON SOLVERS       #
#---------------------------------#

function solve(
    solver_options::MinpackOptions,
    sensitivity_parameters,
    const_cache;
    initial_guess=nothing,
)

    # - Extract constants - #
    (;
        # General
        verbose,
        silence_warnings,
        multipoint_index,
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

    (; algorithm, atol, iteration_limit, converged, iterations) = solver_options

    # - Extract Initial Guess Vector for State Variables - #
    if isnothing(initial_guess)
        initial_guess = extract_initial_guess(
            solver_options, sensitivity_parameters, solve_parameter_cache_dims.state_dims
        )
    end

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

    # - SOLVE - #
    if verbose
        println("  " * "MINPACK Solve Trace:")
    end
    result = MINPACK.fsolve(
        rwrap!,
        jwrap!,
        copy(initial_guess);
        tol=atol,
        tracing=true,
        show_trace=verbose,
        method=algorithm,
        iterations=iteration_limit,
    )

    # update convergence flag
    converged[multipoint_index[]] = result.converged
    iterations[multipoint_index[]] = result.trace.trace[end].iteration

    return result.x
end

function solve(
    solver_options::SIAMFANLEOptions,
    sensitivity_parameters,
    const_cache;
    initial_guess=nothing,
)

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
        resid_cache_vec,
        krylov_cache_vec,
        jvp_cache_vec,
        solve_parameter_cache_dims,
        # Cache(s)
        solve_container_cache,
        solve_container_cache_dims,
    ) = const_cache

    # - Extract Initial Guess Vector for State Variables - #
    if isnothing(initial_guess)
        initial_guess = extract_initial_guess(
            solver_options, sensitivity_parameters, solve_parameter_cache_dims.state_dims
        )
    end

    # - Wrap residual - #
    if verbose
        println("  " * "Wrapping Residual")
    end

    function rwrap!(resid, state_variables, p)
        return system_residual!(
            resid, state_variables, p.sensitivity_parameters, p.constants
        )
    end

    function rwrap(state_variables, p)
        return system_residual(state_variables, p.sensitivity_parameters, p.constants)
    end

    """
        jvp(v,r,x,p)

    Caculate the Jacobian-Vector Product such that JVP = F'(x)*v, or in other words: the directional derivative of F in the direction of v.

    The JVP is calculated as the directional derivative by:

    `g(t) = F(x + t * v, p)`

    `JVP = ForwardDiff.derivative!(p.jvp_cache_vec, g, 0.0)`

    Note that F is defined in the same scope as this sub-function.
    Also note that the input structure of this function is that required by SIAMFANLEquations.jl

    # Arguments:
    - `v::Vector` : the vector v by which F'(x) is multiplied
    - `r::Vector` : unused, but requried by SIAMFANLEquations.jl, it is the storage vector used for an in-place verson of F
    - `x::Vector` : the vector x at which F is evaluated
    - `p::NamedTuple` : the parameters needed by F, as well as a vector storage used in an in-place derivative call from ForwardDiff.derivative!
    """
    function jvp(v, r, x, p)
        g(t) = rwrap(x + t * v, p)
        return ForwardDiff.derivative!(p.jvp_cache_vec, g, 0.0)
    end

    # - SOLVE - #
    if verbose
        println("  " * "Nonlinear Solve Trace:")
    end
    result = solver_options.algorithm(
        rwrap!,
        initial_guess,
        resid_cache_vec,
        krylov_cache_vec,
        jvp;
        atol=solver_options.atol,
        rtol=0.0,
        iteration_limit=solver_options.iteration_limit,
        linear_iteration_limit=solver_options.linear_iteration_limit,
        pdata=(;
            sensitivity_parameters,
            jvp_cache_vec, # TODO: add this to solve_parameter_cache
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
        solver_options.additional_kwargs...,
    )

    # update convergence flag
    # errorcode is 0 if everything is good and could be a bunch of other stuff if not.
    solver_options.converged[multipoint_index[]] = Bool(iszero(result.errcode))

    return result.solution
end

#---------------------------------#
#         NEWTON+ SOLVERS         #
#---------------------------------#
#=
  NOTE: Both SimpleNonlinearSolve and NLsolve also contain fixed-point solvers that do well
=#

function solve(
    solver_options::NonlinearSolveOptions,
    sensitivity_parameters,
    const_cache;
    initial_guess=nothing,
)

    # - Extract constants - #
    (;
        # General
        verbose,
        silence_warnings,
        multipoint_index,
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
        algorithm,
        additional_kwargs,
        # Iteration Controls
        atol,
        iteration_limit,
        converged,
        iterations,
    ) = solver_options

    # - Extract Initial Guess Vector for State Variables - #
    if isnothing(initial_guess)
        initial_guess = extract_initial_guess(
            solver_options, sensitivity_parameters, solve_parameter_cache_dims.state_dims
        )
    end

    # - Wrap residual - #
    if verbose
        println("  " * "Wrapping Residual")
    end

    function rwrap!(resid, state_variables, p)
        return system_residual!(
            resid, state_variables, p.sensitivity_parameters, p.constants
        )
    end

    # build problem object
    prob = SimpleNonlinearSolve.NonlinearProblem(
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
    sol = SimpleNonlinearSolve.solve(
        prob, # problem
        algorithm(additional_kwargs...);
        abstol=atol,
        iteration_limit=iteration_limit,
    )

    # update convergence flag
    converged[multipoint_index[]] = SciMLBase.successful_retcode(sol)

    return sol.u
end

function solve(
    solver_options::NLsolveOptions,
    sensitivity_parameters,
    const_cache;
    initial_guess=nothing,
)

    # - Extract constants - #
    (;
        # General
        verbose,
        silence_warnings,
        multipoint_index,
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
        algorithm,
        linesearch_method,
        linesearch_kwargs,
        atol,
        iteration_limit,
        converged,
        iterations,
    ) = solver_options

    # - Extract Initial Guess Vector for State Variables - #
    if isnothing(initial_guess)
        initial_guess = extract_initial_guess(
            solver_options, sensitivity_parameters, solve_parameter_cache_dims.state_dims
        )
    end

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
        if algorithm == :newton
            println("  " * "Newton Solve Trace:")
        elseif algorithm == :anderson
            println("  " * "Anderson Solve Trace:")
        else
            println("  " * "$(algorithm) Solve Trace:")
        end
    end
    result = NLsolve.nlsolve(
        df,
        initial_guess; # initial states guess
        method=algorithm,
        # autodiff=autodiff,
        linesearch=linesearch_method(linesearch_kwargs...),
        linsolve=(x, A, b) -> x .= ImplicitAD.implicit_linear(A, b),
        ftol=atol,
        iterations=iteration_limit,
        show_trace=verbose,
    )

    # update convergence flag
    # note: need to do this complicated check to ensure that we're only checking |f(x)|<atol rather than claiming convergences when the step size is zero but it's really just stuck.
    _, converged[multipoint_index[]] = NLsolve.assess_convergence(NLsolve.value(df), atol)
    iterations[multipoint_index[]] = result.iterations

    return result.zero
end

#---------------------------------#
#     POLY-ALGORITHM SOLVERS      #
#---------------------------------#

function solve(
    solver_options::CompositeSolverOptions,
    sensitivity_parameters,
    const_cache;
    initial_guess=nothing,
)

    # - Run First Solver - #
    solution = solve(
        solver_options.solvers[1],
        sensitivity_parameters,
        (; const_cache..., solver_options=solver_options.solvers[1]);
        initial_guess=initial_guess,
    )

    # If there is only one solver, return the solution
    if length(solver_options.solvers) == 1
        solver_options.converged[const_cache.multipoint_index[]] = solver_options.solvers[1].converged[const_cache.multipoint_index[]]
        return solution
    end

    # - If there are more than 2 solvers, then loop through the middle ones and return the solution from the last one- #
    if length(solver_options.solvers) > 2
        for (s, sopt) in enumerate(solver_options.solvers[2:(end - 1)])
            solution[:] .= solve(
                sopt,
                sensitivity_parameters,
                (; const_cache..., solver_options=sopt);
                initial_guess=solution,
            )
        end

        solution .= solve(
            solver_options.solvers[end],
            sensitivity_parameters,
            (; const_cache..., solver_options=solver_options.solvers[end]);
            initial_guess=solution,
        )

        solver_options.converged[const_cache.multipoint_index[]] = solver_options.solvers[end].converged[const_cache.multipoint_index[]]
        return solution
    else
        # - If there are only 2 solvers, return the solution from the second one - #
        solution .= solve(
            solver_options.solvers[end],
            sensitivity_parameters,
            (; const_cache..., solver_options=solver_options.solvers[end]);
            initial_guess=solution,
        )

        solver_options.converged[const_cache.multipoint_index[]] = solver_options.solvers[end].converged[const_cache.multipoint_index[]]
        return solution
    end
end

function solve(
    solver_options::ChainSolverOptions,
    sensitivity_parameters,
    const_cache;
    initial_guess=nothing,
)

    # - Run First Solver - #
    solution = solve(
        solver_options.solvers[1],
        sensitivity_parameters,
        (; const_cache..., solver_options=solver_options.solvers[1]);
        initial_guess=initial_guess,
    )

    # If there is only one solver, or if the first solver converged, return the solution
    if length(solver_options.solvers) == 1 ||
        solver_options.solvers[1].converged[const_cache.multipoint_index[]]
        solver_options.converged[const_cache.multipoint_index[]] = solver_options.solvers[1].converged[const_cache.multipoint_index[]]
        return solution
    end

    # - If there are more than 2 solvers, then loop through the middle ones and return the solution from the last one- #
    if length(solver_options.solvers) > 2
        for (s, sopt) in enumerate(solver_options.solvers[2:(end - 1)])
            solution = solve(
                sopt,
                sensitivity_parameters,
                (; const_cache..., solver_options=sopt);
                initial_guess=initial_guess,
            )

            if sopt.converged[1]
                solver_options.converged[const_cache.multipoint_index[]] = sopt.converged[const_cache.multipoint_index[]]
                return solution
            end
        end

        solution = solve(
            solver_options.solvers[end],
            sensitivity_parameters,
            (; const_cache..., solver_options=solver_options.solvers[end]);
            initial_guess=initial_guess,
        )

        solver_options.converged[const_cache.multipoint_index[]] = solver_options.solvers[end].converged[const_cache.multipoint_index[]]
        return solution
    else

        # - If there are only 2 solvers, return the solution from the second one - #
        solution = solve(
            solver_options.solvers[end],
            sensitivity_parameters,
            (; const_cache..., solver_options=solver_options.solvers[end]);
            initial_guess=initial_guess,
        )

        solver_options.converged[const_cache.multipoint_index[]] = solver_options.solvers[end].converged[const_cache.multipoint_index[]]
        return solution
    end
end
