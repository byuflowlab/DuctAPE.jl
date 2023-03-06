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
    vxd_wake_to_body,
    vrd_wake_to_body,
    # vxd_rotor_to_body,
    # vrd_rotor_to_body,
    # rotor_sources,
)
    # - Initialize Right Hand Side - #
    TF = eltype(wake_gammas)
    b = zeros(TF, length(body_boundary_conditions))
    b .+= body_boundary_conditions

    # - Add Wake Contributions to Boundary Conditions - #
    for i in 1:length(vxd_wake_to_body)
        # need to add sum_w[Vx*cos(beta) + Vr*sin(beta)]

        # Vx and Vr are the induced velocities, times the panel length (that product was precomputed as vdx and vrx) times the singularity strength.
        # This product here accomplishes that sum since it's a matrix product.
        Vx = vxd_wake_to_body[i] .* wake_gammas[i, :]
        Vr = vrd_wake_to_body[i] .* wake_gammas[i, :]

        #the added boundary condition is parallel to the panels
        bc = Vx .* cos.(body_panels.panel_angle) .+ Vr .* sin.(body_panels.panel_angle)

        #we subtract since it's on the RHS
        b[1:(end - 1)] .-= bc
    end

    # # - Add Rotor Contributions to Boundary Conditions - #
    # very similar to wakes, but using the source matrices
    for i in 1:length(vxd_rotor_to_body)
        # need to add sum_r[Vx*cos(beta) + Vr*sin(beta)]

        # Vx and Vr are the induced velocities, times the panel length (that product was precomputed as vdx and vrx) times the singularity strength.
        # This product here accomplishes that sum since it's a matrix product.
        Vx = vxd_rotor_to_body[i] .* source_sigmas[i, :]
        Vr = vrd_rotor_to_body[i] .* source_sigmas[i, :]

        #the added boundary condition is parallel to the panels
        bc = Vx .* cos.(body_panels.panel_angle) .+ Vr .* sin.(body_panels.panel_angle)

        #we subtract since it's on the RHS
        b[1:(end - 1)] .-= bc
    end

    # - Return Solution to Linear System - #

    body_vortex_strengths .= ImplicitAD.implicit_linear(A_body_to_body, b)[1:(end - 1)]

    return nothing
end
