"""
    allocate_grid_parameter_cache(pg, x, n)

Allocate a cache used inside the elliptic grid solve.

# Arguments
- `pg::AbstractArray{Float,3}` : the proposed grid array
- `x::AbstractVector{Float}` : the array of ξ values used in the solve
- `n::AbstractVector{Float}` : the array of η values used in the solve

# Returns
- `grid_parameter_cache::NamedTuple` : A named tuple containing the PreallocationTools DiffCache and dimensions for accessing it.
"""
function allocate_grid_parameter_cache(pg, x, n)
    total_length = [0]

    s = size(pg)
    l = lfs(s)
    proposed_grid = cache_dims!(total_length, l, s)

    s = size(x)
    l = lfs(s)
    xi = cache_dims!(total_length, l, s)

    s = size(n)
    l = lfs(s)
    eta = cache_dims!(total_length, l, s)

    return (;
        x_cache=PreallocationTools.DiffCache(zeros(total_length[])),
        x_dims=(; proposed_grid, xi, eta),
    )
end

"""
    withdraw_grid_parameter_cache(vec, dims)

Reshape the cache used inside the elliptic grid solve.

# Arguments
- `vec::Vector{Float}` : the cache vector
- `dims::NamedTuple` : the named tuple of dimensions used to reshape the cache vector

# Returns
- `proposed_grid::AbstractArray{Float,3}` : the proposed grid array
- `xi::AbstractVector{Float}` : the array of ξ values used in the solve
- `eta::AbstractVector{Float}` : the array of η values used in the solve
"""
function withdraw_grid_parameter_cache(vec, dims)
    proposed_grid = reshape(@view(vec[dims.proposed_grid.index]), dims.proposed_grid.shape)

    xi = reshape(@view(vec[dims.xi.index]), dims.xi.shape)

    eta = reshape(@view(vec[dims.eta.index]), dims.eta.shape)

    return proposed_grid, xi, eta
end
