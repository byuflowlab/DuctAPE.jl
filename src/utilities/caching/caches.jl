"""
    initialize_all_caches(solver_options, paneling_constants)

Convenience function to initialize all caches before calling analysis.

# Arguments
- `solver_options::SolverOptionsType` : solver options used for cache allocation dispatch
- `paneling_constants::PanelingConstants` : PanelingConstants object upon which all cache sizing depends

# Returns
- `prepost_container_caching::NamedTuple` : A named tuple containing the PreallocationTools DiffCache and a named tuple with relevant dimensions for accessing the cache.
- `solve_parameter_caching::NamedTuple` : A named tuple containing the PreallocationTools DiffCache and a named tuple with relevant dimensions for accessing the cache.
- `solve_container_caching::NamedTuple` : A named tuple containing the PreallocationTools DiffCache and a named tuple with relevant dimensions for accessing the cache.
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
