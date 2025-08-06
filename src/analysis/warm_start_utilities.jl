function find_large_and_nearest_small(residuals)
    # Indices of large and small residuals
    large_indices = findall(>(1), residuals)
    small_indices = findall(<(1), residuals)

    # Precompute for speed
    small_array = collect(small_indices)  # ensure it's indexable

    # Result: tuple of (large_index, closest_small_index)
    pairs = [(li, small_array[argmin(abs.(small_array .- li))]) for li in large_indices]

    return sort(pairs; by=x -> abs(x[1] - x[2]))
end

function find_false_and_nearest_true_sorted(converged)
    false_indices = findall(!, converged)
    true_indices = findall(identity, converged)
    true_array = collect(true_indices)  # make sure it's indexable

    pairs = [(fi, true_array[argmin(abs.(true_array .- fi))]) for fi in false_indices]

    # Sort by proximity
    return sort(pairs; by=x -> abs(x[1] - x[2]))
end

