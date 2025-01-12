function get_induced_velocities!(
    velocities,
    points,
    v_bs,
    v_ws,
    v_rs,
    body_vortex_panels,
    body_vortex_strengths,
    wake_vortex_panels,
    wake_vortex_strengths,
    rotor_source_panels,
    rotor_source_strengths,
    Vinf,
    integration_options,
)
    v_bs .= 0.0
    v_ws .= 0.0
    v_rs .= 0.0

    # body-induced velocity
    DuctAPE.induced_velocities_from_vortex_panels_on_points!(
        v_bs,
        points,
        body_vortex_panels.node,
        body_vortex_panels.nodemap,
        body_vortex_panels.influence_length,
        integration_options,
    )

    velocities[:, 1] .+= v_bs[:, :, 1] * body_vortex_strengths
    velocities[:, 2] .+= v_bs[:, :, 2] * body_vortex_strengths

    # wake-induced velocity
    DuctAPE.induced_velocities_from_vortex_panels_on_points!(
        v_ws,
        @view(points[:, :]),
        wake_vortex_panels.node,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.influence_length,
        integration_options,
    )

    velocities[:, 1] .+= v_ws[:, :, 1] * wake_vortex_strengths
    velocities[:, 2] .+= v_ws[:, :, 2] * wake_vortex_strengths

    # rotor-induced velocity
    DuctAPE.induced_velocities_from_source_panels_on_points!(
        v_rs,
        @view(points[:, :]),
        rotor_source_panels.node,
        rotor_source_panels.nodemap,
        rotor_source_panels.influence_length,
        integration_options,
    )

    velocities[:, 1] .+= v_rs[:, :, 1] * rotor_source_strengths
    velocities[:, 2] .+= v_rs[:, :, 2] * rotor_source_strengths

    return nothing
end

function calculate_velocity_field(
    field_points,
    body_vortex_panels,
    body_vortex_strengths,
    wake_vortex_panels,
    wake_vortex_strengths,
    rotor_source_panels,
    rotor_source_strengths,
    Vinf;
    integration_options=IntegrationOptions(),
)

    TF = promote_type(
        eltype(body_vortex_panels.node),
        eltype(wake_vortex_panels.node),
        eltype(rotor_source_panels.node),
    )

    # Initialize
    velocities = zeros(TF, length(field_points[1, :]), 2)

    # So we don't have to re-allocate all the memory:
    v_bs = zeros(TF, length(field_points[1,:]), Int(body_vortex_panels.totnode[]), 2)
    v_ws = zeros(TF, length(field_points[1,:]), Int(wake_vortex_panels.totnode[]), 2)
    v_rs = zeros(TF, length(field_points[1,:]), Int(rotor_source_panels.totnode[]), 2)

    velocities[:, 1] .= Vinf

    # - Calculate Velocity at Current Point - #
    get_induced_velocities!(
        velocities,
        field_points,
        v_bs,
        v_ws,
        v_rs,
        body_vortex_panels,
        body_vortex_strengths,
        wake_vortex_panels,
        wake_vortex_strengths,
        rotor_source_panels,
        rotor_source_strengths,
        Vinf,
        integration_options,
    )

    return velocities[:, 1], velocities[:, 2]
end
