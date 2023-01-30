"""
assumes inputs are tuples
"""
function calculate_induced_velocity(coefficient_matrix, singularity_strengths)

    # Initialize output
    induced_velocity = 0.0

    for i in 1:length(coefficient_matrix)
        induced_velocity += coefficient_matrix[i] * singularity_strengths[i]
    end

    return induced_velocity
end

