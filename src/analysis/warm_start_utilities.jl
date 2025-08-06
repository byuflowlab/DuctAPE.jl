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
    flat = vec(converged)  # Convert to Vector{Bool}

    false_indices = findall(!, flat)
    true_indices = findall(identity, flat)

    pairs = [
        (
            fi,
            true_indices[argmin(abs.(true_indices .- fi))]
        )
        for fi in false_indices
    ]

    # Sort by proximity
    sort(pairs, by = x -> abs(x[1] - x[2]))
end
