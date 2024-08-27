"""
    initialize_all_caches(solver_options, paneling_constants)

Convenience function to initialize all caches before calling analysis.

# Arguments
- `solver_options::SolverOptionsType` : solver options used for cache allocation dispatch
- `paneling_constants::PanelingConstants` : PanelingConstants object upon which all cache sizing depends

# Keyword Arguments
- `fd_chunk_size::Int=12` : chunk size to use for PreallocationTools caches.  Note that the automated chunk size for DuctAPE will always be the ForwardDiff threshold of 12 due to the size of the system, so it will be best to leave this at the default unless further development allows for chunk size selection for individual solvers.
- `levels::Int=1` : levels for nested duals.  Note that since ImplicitAD is being used for all solves, there should be no need for more than 1 level.


# Returns
- `prepost_container_caching::NamedTuple` : A named tuple containing the PreallocationTools DiffCache and a named tuple with relevant dimensions for accessing the cache.
- `solve_parameter_caching::NamedTuple` : A named tuple containing the PreallocationTools DiffCache and a named tuple with relevant dimensions for accessing the cache.
- `solve_container_caching::NamedTuple` : A named tuple containing the PreallocationTools DiffCache and a named tuple with relevant dimensions for accessing the cache.
"""
function initialize_all_caches(
    solver_options, paneling_constants; fd_chunk_size=12, levels=1
)

    # - Pre/Post Containers Cache - #
    prepost_container_caching = allocate_prepost_container_cache(
        paneling_constants; fd_chunk_size=fd_chunk_size, levels=levels
    )

    # - Solve (Sensitivity) Parameters Cache - #
    solve_parameter_caching = allocate_solve_parameter_cache(
        solver_options, paneling_constants; fd_chunk_size=fd_chunk_size, levels=levels
    )

    # - Solve Intermediate Containers Cache - #
    solve_container_caching = allocate_solve_container_cache(
        solver_options, paneling_constants; fd_chunk_size=fd_chunk_size, levels=levels
    )

    return prepost_container_caching, solve_parameter_caching, solve_container_caching
end
