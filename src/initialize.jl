
"""
"""
function initialize_parameters(
    duct_coordinates, hub_coordinates, rotor_parameters, stator_parameters, freestream
)

    #---------------------------------#
    #             Poblem              #
    #---------------------------------#
    method = ff.AxisymmetricProblem(Vortex(Constant()), Neumann(), [false, true])
    problem = ff.define_problem(
        method, (duct_coordinates, hub_coordinates), 0.0, -1.0, -1.0
    )

    #---------------------------------#
    #             Bodies              #
    #---------------------------------#
    body_geometry, body_panels = generate_body_geometry(duct_coordinates, hub_coordinates)
    mesh_body_to_body = ff.generate_mesh(method, body_panels)
    body_system = ff.generate_inviscid_system(method, body_panels, mesh_body_to_body)

    #---------------------------------#
    #           First Rotor           #
    #---------------------------------#
    rotor_blade_elements, rotor_panels = generate_blade_elements(
        rotor_parameters.xpos,
        rotor_parameters.radial_positions,
        rotor_parameters.chords,
        rotor_parameters.twists,
        rotor_parameters.airfoils,
        rotor_parameters.num_blade_elements,
        rotor_parameters.num_blades,
        rotor_parameters.omega,
        body_geometry,
    )

    #---------------------------------#
    #            Wake Points          #
    #---------------------------------#
    x_grid_points, r_grid_points, nx, nr, rotoridxs, wake_panels = generate_wake_grid(
        body_geometry,
        [rotor_parameters.xpos; stator_parameters.xpos],
        rotor_blade_elements.radial_positions;
        wake_length=1.0,
        debug=false,
    )

    #---------------------------------#
    #          Other Rotor(s)         #
    #---------------------------------#
    stator_blade_elements, stator_panels = generate_blade_elements(
        stator_parameters.xpos,
        stator_parameters.radial_positions,
        stator_parameters.chords,
        stator_parameters.twists,
        stator_parameters.airfoils,
        stator_parameters.num_blade_elements,
        stator_parameters.num_blades,
        stator_parameters.omega,
        body_geometry;
        updated_radial_positions=r_grid_points[rotoridxs[2], :],
    )

    return (
        converged=[false], # Initialize Convergence Flag
        body_geometry=body_geometry,
        body_panels=body_panels,
        rotor_blade_elements=rotor_blade_elements,
        rotor_panels=rotor_panels,
        stator_blade_elements=stator_blade_elements,
        stator_panels=stator_panels,
        rotoridxs=rotoridxs,
        wake_panels=wake_panels,
        #---------------------------------#
        #             Meshes              #
        #---------------------------------#
        # - Body -> Body - #
        mesh_body_to_body=mesh_body_to_body,
        # - Body -> Rotor - #
        mesh_bodies_to_rotor=generate_one_way_mesh(body_panels, rotor_panels),
        # - Body -> Stator - #
        mesh_bodies_to_stator=generate_one_way_mesh(body_panels, stator_panels),
        # - Wake -> Body - #
        mesh_wake_to_body=generate_one_way_mesh(wake_panels, body_panels),
        # - Wake -> Rotor - #
        mesh_wake_to_rotor=generate_one_way_mesh(wake_panels, rotor_panels),
        # - Wake -> Stator - #
        mesh_wake_to_stator=generate_one_way_mesh(wake_panels, stator_panels),
        # - Rotor -> Body - #
        mesh_rotor_to_body=generate_one_way_mesh(
            rotor_panels, body_panels; singularity="source"
        ),
        # - Rotor -> Rotor - #
        mesh_rotor_to_rotor=generate_one_way_mesh(
            rotor_panels, rotor_panels; singularity="source"
        ),
        # - Rotor -> Stator - #
        mesh_rotor_to_stator=generate_one_way_mesh(
            rotor_panels, stator_panels; singularity="source"
        ),
        # - Stator -> Body - #
        mesh_stator_to_body=generate_one_way_mesh(
            stator_panels, body_panels; singularity="source"
        ),
        # - Stator -> Rotor - #
        mesh_stator_to_rotor=generate_one_way_mesh(
            stator_panels, rotor_panels; singularity="source"
        ),
        # - Stator -> Stator - #
        mesh_stator_to_stator=generate_one_way_mesh(
            stator_panels, stator_panels; singularity="source"
        ),
        #---------------------------------#
        #       Coefficient Matrices      #
        #---------------------------------#
        # - Body -> Body - #
        A_body_to_body=body_system.A,
        bc_freestream_to_body=body_system.b,
        # - Body -> Rotor - #
        A_bodies_to_rotor=assemble_one_way_coefficient_matrix(
            mesh_bodies_to_rotor, body_panels, rotor_panels
        ),
        # - Body -> Stator - #
        A_bodies_to_stator=assemble_one_way_coefficient_matrix(
            mesh_bodies_to_stator, body_panels, stator_panels
        ),
        # - Wake -> Body - #
        A_wake_to_bodies=assemble_one_way_coefficient_matrix(
            mesh_wake_to_body, wake_panels, body_panels
        ),
        # - Wake -> Rotor - #
        A_wake_to_rotor=assemble_one_way_coefficient_matrix(
            mesh_wake_to_rotor, wake_panels, stator_panels
        ),
        # - Wake -> Stator - #
        A_wake_to_stator=assemble_one_way_coefficient_matrix(
            mesh_wake_to_stator, wake_panels, stator_panels
        ),
        # - Rotor -> Body - #
        A_rotor_to_body=assemble_one_way_coefficient_matrix(
            mesh_rotor_to_body, rotor_panels, body_panels; singularity="source"
        ),
        # - Rotor -> Rotor - #
        A_rotor_to_rotor=assemble_one_way_coefficient_matrix(
            mesh_rotor_to_rotor, rotor_panels, rotor_panels; singularity="source"
        ),
        # - Rotor -> Stator - #
        A_rotor_to_stator=assemble_one_way_coefficient_matrix(
            mesh_rotor_to_stator, rotor_panels, stator_panels; singularity="source"
        ),
        # - Stator -> Body - #
        A_stator_to_body=assemble_one_way_coefficient_matrix(
            mesh_stator_to_body, stator_panels, body_panels; singularity="source"
        ),
        # - Stator -> Rotor - #
        A_stator_to_rotor=assemble_one_way_coefficient_matrix(
            mesh_stator_to_rotor, stator_panels, rotor_panels; singularity="source"
        ),
        # - Stator -> Stator - #
        A_stator_to_stator=assemble_one_way_coefficient_matrix(
            mesh_stator_to_stator, stator_panels, stator_panels; singularity="source"
        ),
    )
end

"""
get first guess for gamma and sigma values
"""
function gamma_sigma_guess()
    return nothing
end
