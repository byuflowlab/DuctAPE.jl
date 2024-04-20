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

function withdraw_grid_parameter_cache(vec, dims)
    proposed_grid = reshape(@view(vec[dims.proposed_grid.index]), dims.proposed_grid.shape)

    xi = reshape(@view(vec[dims.xi.index]), dims.xi.shape)

    eta = reshape(@view(vec[dims.eta.index]), dims.eta.shape)

    return proposed_grid, xi, eta
end
