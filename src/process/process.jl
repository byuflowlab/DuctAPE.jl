"""
"""
function process(
    solver_options::TS,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
    airfoils,
    A_bb_LU,
    solve_container_caching,
    idmaps,
    options,
) where {TS<:Union{ExternalSolverOptions,PolyAlgorithmOptions}}

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

    # - combine cache and constants - #
    const_cache = (;
        # - Constants - #
        options.verbose,
        options.silence_warnings,
        options.multipoint_index,
        #nlsolve options
        solver_options,
        # Constant Parameters
        airfoils,
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
"""
function process(
    solver_options::CSORSolverOptions,
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
        solver_options, solve_parameter_cache_vector, solve_parameter_cache_dims
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

    # - combine cache and constants - #
    const_cache = (;
        # - Constants - #
        options.silence_warnings,
        options.verbose,
        options.multipoint_index,
        #CSOR solve options
        solver_options,
        # Constant Parameters
        airfoils,
        A_bb_LU,
        idmaps,
        # Cache(s)
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        solve_container_caching...,
    )

    # - Solve with ImplicitAD - #
    if options.verbose
        println("\nSolving Nonlinear System using CSOR Method")
    end
    return ImplicitAD.implicit(
        solve, CSOR_residual!, solve_parameter_cache_vector, const_cache
    )
end
