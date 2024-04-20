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

"""
"""
function cache_dims!(total_length, l, s)
    dims = (; index=(total_length[] + 1):(total_length[] + l), shape=s)
    total_length[] += l
    return dims
end
