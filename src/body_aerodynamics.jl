#=

Functions regarding body aerodyanmics

=#

"""
"""
function solve_linear_system!(
    body_vortex_strengths,
    A_body_to_body,
    Vinf,
    body_boundary_conditions,
    wake_gammas,
    A_wake_to_bodies,
)#, Sigmas, A_rotor_to_body)

    # - Initialize Right Hand Side - #
    b = -body_boundary_conditions .* Vinf

    # - Add Wake Contributions to Boundary Conditions - #
    for i in 1:length(A_wake_to_bodies)
        b .-= A_wake_to_bodies[i] * wake_gammas[i, :]
    end

    # # - Add Rotor Contributions to Boundary Conditions - #
    # for i in 1:length(Sigmas[1, :])
    #     b .-= A_rotor_to_body[i] * Sigmas[:, i]'
    # end

    # - Return Solution to Linear System - #
    body_vortex_strengths .= ImplicitAD.implicit_linear(A_body_to_body, b)

    return nothing
end
