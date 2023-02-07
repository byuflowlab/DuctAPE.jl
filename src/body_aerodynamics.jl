#=

Functions regarding body aerodyanmics

=#

"""
"""
function calculate_body_vortex_strengths!(
    body_vortex_strengths,
    A_body_to_body,
    body_boundary_conditions,
    wake_gammas,
    A_wake_to_body,
)#, Sigmas, A_rotor_to_body)

    # - Initialize Right Hand Side - #
    TF = eltype(wake_gammas)
    b = zeros(TF, length(body_boundary_conditions))
    b .+= body_boundary_conditions

    # - Add Wake Contributions to Boundary Conditions - #
    for i in 1:length(A_wake_to_body)
        b[1:(end - 1)] .-= A_wake_to_body[i] * wake_gammas[i, :]
    end

    # # - Add Rotor Contributions to Boundary Conditions - #
    # for i in 1:length(Sigmas[1, :])
    #     b .-= A_rotor_to_body[i] * Sigmas[:, i]'
    # end

    # - Return Solution to Linear System - #

    body_vortex_strengths .= ImplicitAD.implicit_linear(A_body_to_body, b)[1:(end - 1)]

    return nothing
end
