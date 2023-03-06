"""
assumes inputs are tuples
"""
function calculate_induced_velocity(vxd, vrd, singularity_strengths)

    # Initialize output
    Vx = 0.0
    Vr = 0.0

    for i in 1:length(vxd)
        Vx += vxd[i] * singularity_strengths[i]
        Vr += vrd[i] * singularity_strengths[i]
    end

    return Vx, Vr
end

