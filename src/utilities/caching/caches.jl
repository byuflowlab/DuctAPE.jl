"""
"""
function initialize_all_caches(solver_options, paneling_constants)

    # - Pre/Post Containers Cache - #
    prepost_container_caching = allocate_prepost_container_cache(paneling_constants)

    # - Solve (Sensitivity) Parameters Cache - #
    solve_parameter_caching = allocate_solve_parameter_cache(
        solver_options, paneling_constants
    )

    # - Solve Intermediate Containers Cache - #
    solve_container_caching = allocate_solve_container_cache(
        solver_options, paneling_constants
    )

    return prepost_container_caching, solve_parameter_caching, solve_container_caching
end