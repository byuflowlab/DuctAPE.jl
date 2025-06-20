"""
    process(
        solver_options::SolverOptionsType,
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        A_bb_LU,
        solve_container_caching,
        idmaps,
        options,
    )

Process (the step between pre-process and post-process) the solution, in other words: call the solver(s).

# Arguments
- `solver_options::SolverOptionsType` : the solver options contained in the options object, used for dispatch.
- `solve_parameter_cache_vector::Vector{Float}` : The vector cache for parameters used in the solve.
- `solve_parameter_cache_dims::NamedTuple` : A named tuple containing the dimensions of the solve parameters.
- `airfoils::NamedTuple` : The airfoils to be interpolated that are associated with each blade element
- `A_bb_LU::LinearAlgebra.LU` : The LU decomposition of the panel method LHS matrix
- `solve_container_caching::NamedTuple` : A named tuple containing the cache and dimensions for the intermediate solve values.
- `idmaps::NamedTuple` : The set of index maps used in various solve sub-functions
- `options::Options` : User options

# Returns
- `converged_states::Vector{Float}` : The output of a call to `ImplicitAD.implicit`
"""
function process(
    solver_options::TS,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
    A_bb_LU,
    solve_container_caching,
    idmaps,
    options,
) where {TS<:ExternalSolverOptions}

    # - Initialize Aero - #
    if options.verbose
        println("\nInitializing Velocities")
    end
    # view the initial conditions out of the inputs cache
    solve_parameter_tuple = withdraw_solve_parameter_cache(
        solver_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )
    (; vz_rotor, vtheta_rotor, Cm_wake, operating_point, linsys, ivr, ivw) =
        solve_parameter_tuple

    # initialize velocities
    # TODO; add some sort of unit test for this function
    initialize_velocities!(
        solver_options,
        vz_rotor,
        vtheta_rotor,
        Cm_wake,
        solve_parameter_tuple.operating_point,
        solve_parameter_tuple.blade_elements,
        (; solve_parameter_tuple.linsys..., A_bb_LU),
        solve_parameter_tuple.ivr,
        solve_parameter_tuple.ivw,
        idmaps.body_totnodes,
        idmaps.wake_panel_sheet_be_map,
    )

    # - combine cache and constants - #
    const_cache = (;
        # - Constants - #
        options.verbose,
        options.silence_warnings,
        options.multipoint_index,
        #nlsolve options
        solver_options,
        # Constant Parameters
        A_bb_LU,
        idmaps,
        solve_parameter_cache_dims,
        # Cache(s)
        solve_container_caching...,
        solve_parameter_tuple..., # contains SIAMFANLEOptions containers needed for solver definition
    )

    # - Solve with ImplicitAD - #
    if options.verbose
        println("\nSolving Nonlinear System")
    end

    return ImplicitAD.implicit(
        solve, system_residual!, solve_parameter_cache_vector, const_cache
    )
end

"""
    process(
            solver_options::CSORSolverOptions,
            solve_parameter_cache_vector,
            solve_parameter_cache_dims,
            A_bb_LU,
            solve_container_caching,
            idmaps,
            options,
        )

Runs the full initialization and nonlinear solve sequence for the CSOR aerodynamic solver.

# Arguments
- `solver_options::CSORSolverOptions`: Configuration parameters for the CSOR solver.
- `solve_parameter_cache_vector`: Vector containing cached solve parameters.
- `solve_parameter_cache_dims`: Dimensions or metadata describing the cache vector structure.
- `A_bb_LU`: LU factorization of the body-body influence matrix for efficient linear solves.
- `solve_container_caching`: Additional cached containers needed for the solver.
- `idmaps`: Named tuple of index maps for panel and wake bookkeeping.
- `options`: User options including verbosity flags and multipoint indexing.

# Returns
- `converged_states::Vector{Float}` : the states for which the residual has converged
"""
function process(
    solver_options::CSORSolverOptions,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
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
        solver_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )
    (; Gamr, sigr, gamw, operating_point, blade_elements, linsys, ivr, ivw, wakeK) =
        solve_parameter_tuple

    # - Initialize States - #
    initialize_strengths!(
        solver_options,
        Gamr,
        sigr,
        gamw,
        operating_point,
        blade_elements,
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
        idmaps.wake_node_ids_along_center_body_wake_interface,
    )

    # - combine cache and constants - #
    const_cache = (;
        # - Constants - #
        options.silence_warnings,
        options.verbose,
        options.multipoint_index,
        #CSOR solve options
        solver_options,
        # Constant Parameters
        A_bb_LU,
        idmaps,
        # Cache(s)
        # solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        solve_container_caching...,
    )

    # - Solve with ImplicitAD - #
    if options.verbose
        println("\nSolving Nonlinear System using CSOR Method")
    end

    return solve(solver_options, solve_parameter_cache_vector, const_cache)
end

"""
    process(
        solver_options::ModCSORSolverOptions,
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        A_bb_LU,
        solve_container_caching,
        idmaps,
        options,
    )

Runs the full initialization and nonlinear solve sequence for the modified CSOR aerodynamic solver
using the Implicit Automatic Differentiation (ImplicitAD) framework.

# Description
This function performs these steps:
1. Extracts initial flow state variables and parameters from a cached input vector.
2. Initializes vortex and source strength states for rotor and wake panels.
3. Combines solver constants, options, and cached data into a single tuple for the solver.
4. Invokes the nonlinear solver using the ImplicitAD framework with a modified CSOR residual.

# Arguments
- `solver_options::ModCSORSolverOptions`: Configuration parameters for the modified CSOR solver.
- `solve_parameter_cache_vector`: Vector holding cached solver parameters and states.
- `solve_parameter_cache_dims`: Metadata describing the structure of the cached vector.
- `A_bb_LU`: LU factorization of the body-body influence matrix for efficient solves.
- `solve_container_caching`: Additional cached containers required by the solver.
- `idmaps`: Named tuple of index maps for panel and wake bookkeeping.
- `options`: User-defined options controlling verbosity, multipoint indexing, and warning behavior.

# Returns
- `converged_states::Vector{Float}` : The output of a call to `ImplicitAD.implicit`
"""
function process(
    solver_options::ModCSORSolverOptions,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
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
        solver_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )

    (; Gamr, sigr, gamw, operating_point, blade_elements, linsys, ivr, ivw, wakeK) =
        solve_parameter_tuple

    # - Initialize States - #
    initialize_strengths!(
        solver_options,
        Gamr,
        sigr,
        gamw,
        operating_point,
        blade_elements,
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
        idmaps.wake_node_ids_along_center_body_wake_interface,
    )

    # - combine cache and constants - #
    const_cache = (;
        # - Constants - #
        options.silence_warnings,
        options.verbose,
        options.multipoint_index,
        #solve options
        solver_options,
        # Constant Parameters
        A_bb_LU,
        idmaps,
        # Cache(s)
        solve_parameter_cache_dims,
        solve_container_caching...,
    )

    # - Solve with ImplicitAD - #
    if options.verbose
        println("\nSolving Nonlinear System using Modified CSOR Method")
    end

    return ImplicitAD.implicit(
        solve, mod_CSOR_residual!, solve_parameter_cache_vector, const_cache
    )
end
