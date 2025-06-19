function wake_streamlines()

    # set body streamlines
    # determine center_body wake streamline
    # determine inner rotor streamlines
    # determine duct wake streamline

    return nothing
end

function get_streamline_velocities!(
    velocities,
    points,
    z,
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
        @view(points[:, z, :]),
        body_vortex_panels.node,
        body_vortex_panels.nodemap,
        body_vortex_panels.influence_length,
        integration_options,
    )

    velocities[z, :, 1] .+= v_bs[:, :, 1] * body_vortex_strengths
    velocities[z, :, 2] .+= v_bs[:, :, 2] * body_vortex_strengths

    # wake-induced velocity
    DuctAPE.induced_velocities_from_vortex_panels_on_points!(
        v_ws,
        @view(points[:, z, :]),
        wake_vortex_panels.node,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.influence_length,
        integration_options,
    )

    velocities[z, :, 1] .+= v_ws[:, :, 1] * wake_vortex_strengths
    velocities[z, :, 2] .+= v_ws[:, :, 2] * wake_vortex_strengths

    # rotor-induced velocity
    DuctAPE.induced_velocities_from_source_panels_on_points!(
        v_rs,
        @view(points[:, z, :]),
        rotor_source_panels.node,
        rotor_source_panels.nodemap,
        rotor_source_panels.influence_length,
        integration_options,
    )

    velocities[z, :, 1] .+= v_rs[:, :, 1] * rotor_source_strengths
    velocities[z, :, 2] .+= v_rs[:, :, 2] * rotor_source_strengths

    return nothing
end

function calculate_streamlines(
    body_vortex_panels,
    body_vortex_strengths,
    wake_vortex_panels,
    wake_vortex_strengths,
    rotor_source_panels,
    rotor_source_strengths,
    Vinf;
    starting_radial_points=range(0.001, 1.0; length=20),
    axial_range=[0.5, 1.0],
    nominal_step_size=1e-2,
    step_limit=Int(1e2),
    integration_options=IntegrationOptions(),
    dot_tol=0.999,
    stag_tol=0.4,
)

    # TODO;
    # figure out how to ignore things after a stagnation point is hit, maybe set it to NaN?

    TF = promote_type(
        eltype(body_vortex_panels.node),
        eltype(wake_vortex_panels.node),
        eltype(rotor_source_panels.node),
    )

    # Initialize
    points = zeros(TF, 2, step_limit, length(starting_radial_points))
    points[1, 1, :] .= axial_range[1]
    points[2, 1, :] = starting_radial_points
    velocities = zeros(TF, step_limit, length(starting_radial_points), 2)

    # So we don't have to allocate all the memory:
    v_bs = zeros(TF, length(starting_radial_points), Int(body_vortex_panels.totnode[]), 2)
    v_ws = zeros(TF, length(starting_radial_points), Int(wake_vortex_panels.totnode[]), 2)
    v_rs = zeros(TF, length(starting_radial_points), Int(rotor_source_panels.totnode[]), 2)

    # March along each streamline
    for z in 1:step_limit
        velocities[z, :, 1] .= Vinf

        # - Calculate Velocity at Current Point - #
        get_streamline_velocities!(
            velocities,
            points,
            z,
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

        if z > 1
            dp = [
                dot(
                    velocities[z, r, :] / norm(velocities[z, r, :]),
                    velocities[z - 1, r, :] / norm(velocities[z - 1, r, :]),
                ) for r in eachindex(starting_radial_points)
            ]
            while all(dp .< dot_tol)
                points[:, z, :] .= 0.5 .* (points[:, z, :] .+ points[:, z - 1, :])
                velocities[z, :, 1] .= Vinf

                get_streamline_velocities!(
                    velocities,
                    points,
                    z,
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

                dp .= [
                    dot(
                        velocities[z, r, :] / norm(velocities[z, r, :]),
                        velocities[z - 1, r, :] / norm(velocities[z - 1, r, :]),
                    ) for r in eachindex(starting_radial_points)
                ]
            end
        end

        if z != step_limit
            # - Take Step and Set Next Point - #
            vmag = sqrt.(velocities[z, :, 2] .^ 2 .+ velocities[z, :, 1] .^ 2)
            points[1, z + 1, :] .=
                points[1, z, :] .+ velocities[z, :, 1] .* nominal_step_size ./ vmag
            points[2, z + 1, :] .=
                points[2, z, :] .+ velocities[z, :, 2] .* nominal_step_size ./ vmag
            if z > 1
                for r in eachindex(starting_radial_points)
                    if vmag[r] < stag_tol
                        points[:, z + 1, r] .= NaN
                    end
                end
            end
        end
    end

    return points, velocities
end
