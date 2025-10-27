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
    if !solver_options.warm_start[options.multipoint_index[]]
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

        # if skip solve, just return states from initialization
        if solver_options.skip_solve[options.multipoint_index[]]
            return extract_initial_guess(
                solver_options,
                solve_parameter_cache_vector,
                solve_parameter_cache_dims.state_dims,
            )
        end
    end

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
    if !solver_options.warm_start[options.multipoint_index[]]
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

        # if skip solve, just return states from initialization
        if solver_options.skip_solve[options.multipoint_index[]]
            return extract_initial_guess(
                solver_options,
                solve_parameter_cache_vector,
                solve_parameter_cache_dims.state_dims,
            )
        end
    end

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
    if !solver_options.warm_start[options.multipoint_index[]]
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

        # if skip solve, just return states from initialization
        if solver_options.skip_solve[options.multipoint_index[]]
            return extract_initial_guess(
                solver_options,
                solve_parameter_cache_vector,
                solve_parameter_cache_dims.state_dims,
            )
        end
    end

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

    # println("checking for inf's and nan's in inputs to solver in process.jl line 223:")
    if any(isinf.(solve_parameter_cache_vector))
        print(
            "Infs found in solve parameter cache vector at indices: ",
            findall(isinf, solve_parameter_cache_vector),
        )
    end
    if any(isnan.(solve_parameter_cache_vector))
        print(
            "NaNs found in solve parameter cache vector at indices: ",
            findall(isnan, solve_parameter_cache_vector),
        )
    end

    return ImplicitAD.implicit(
        solve, mod_CSOR_residual!, solve_parameter_cache_vector, const_cache
    )
end
