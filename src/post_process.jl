#=

Various Post-processing functions

=#

"""
NEED TO TEST. NO GUARENTEES THIS OR CONSTITUENT FUNCTIONS ARE CORRECT
"""
function probe_velocity_field(
    field_points,
    Vinf;
    body_strengths=nothing,
    body_panels=nothing,
    wake_strengths=nothing,
    wake_panels=nothing,
    source_strengths=nothing,
    source_panels=nothing,
)

    # Initialize velocities
    Vfield_x = similar(field_points) .= Vinf
    body_induced_x = similar(field_points) .= 0.0
    wake_induced_x = similar(field_points) .= 0.0
    source_induced_x = similar(field_points) .= 0.0
    Vfield_r = similar(field_points) .= 0.0
    body_induced_r = similar(field_points) .= 0.0
    wake_induced_r = similar(field_points) .= 0.0
    source_induced_r = similar(field_points) .= 0.0


    # - add body induced velocity to total
    if !isnothing(body_strengths) && !isnothing(body_panels)
        # set up influence coefficients based on field points and body panels
        mesh_fpb = generate_field_mesh(body_panels, field_points)

        A_fpb = assemble_induced_velocity_matrices_infield(
            mesh_fpb, body_panels, field_points
        )

        # axial components
        vx_s = A_fpb[1]

        # radial components
        vr_s = A_fpb[2]

        # Mutliply things out to get induced velocities
        body_induced_x = vx_s * body_strengths
        body_induced_r = vx_s * body_strengths

        # Add to total velocity field
        Vfield_x .+= body_induced_x
        Vfield_r .+= body_induced_r
    end

    # - add wake induced velocity to total
    if !isnothing(wake_strengths) && !isnothing(wake_panels)
        # set up influence coefficients based on field points and body panels
        mesh_fpw = [
            generate_field_mesh([wake_panels[j]], field_points) for i in 1:1,
            j in 1:length(wake_panels)
        ]

        A_fpw = [
            assemble_induced_velocity_matrices_infield(
                mesh_fpw[i, j], [wake_panels[j]], field_points
            ) for i in 1:1, j in 1:length(wake_panels)
        ]

        # axial components
        vx_s = [A_fpw[i, j][1] for i in 1:1, j in 1:length(wake_panels)]

        # radial components
        vr_s = [A_fpw[i, j][2] for i in 1:1, j in 1:length(wake_panels)]

        # Mutliply things out to get induced velocities
        wake_induced_x = vx_s * wake_strengths
        wake_induced_r = vx_s * wake_strengths

        # Add to total velocity field
        Vfield_x .+= wake_induced_x
        Vfield_r .+= wake_induced_r
    end

    # - add source induced velocity to total
    if !isnothing(source_strengths) && !isnothing(source_panels)
        # set up influence coefficients based on field points and body panels
        mesh_fps = [
            generate_field_mesh([source_panels[j]], field_points) for i in 1:1,
            j in 1:length(source_panels)
        ]

        A_fps = [
            assemble_induced_velocity_matrices_infield(
                mesh_fps[i, j], [source_panels[j]], field_points
            ) for i in 1:1, j in 1:length(source_panels)
        ]

        # axial components
        vx_s = [A_fps[i, j][1] for i in 1:1, j in 1:length(source_panels)]

        # radial components
        vr_s = [A_fps[i, j][2] for i in 1:1, j in 1:length(source_panels)]

        # Mutliply things out to get induced velocities
        source_induced_x = vx_s * source_strengths
        source_induced_r = vx_s * source_strengths

        # Add to total velocity field
        Vfield_x .+= source_induced_x
        Vfield_r .+= source_induced_r
    end

    return Vfield_x,
    Vfield_r,
    body_induced_x,
    body_induced_r,
    wake_induced_x,
    wake_induced_r,
    source_induced_x,
    source_induced_r
end
