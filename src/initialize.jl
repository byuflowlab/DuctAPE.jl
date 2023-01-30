
"""
"""
function initialize_parameters(
    duct_coordinates, hub_coordinates, rotor_parameters, freestream
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
    num_rotors = length(rotor_parameters)
    blade_elements, rotor_panels = initialize_blade_elements(
        rotor_parameters[1], body_geometry, num_rotors
    )

    # - Generate Rotor Panels - #
    rotor_panels = initialize_rotor_panels(blade_elements[1], num_rotors; method=method)

    #---------------------------------#
    #            Wake Points          #
    #---------------------------------#
    x_grid_points, r_grid_points, nx, nr, rotoridxs, wake_panels = generate_wake_grid(
        body_geometry,
        (p -> p.xpos).(rotor_parameters),
        rotor_blade_elements.radial_positions;
        wake_length=1.0,
        debug=false,
    )

    wake_grid = [[x_grid_points[i, j] r_grid_points[i, j]] for i in 1:nx, j in 1:nr]

    #---------------------------------#
    #          Other Rotor(s)         #
    #---------------------------------#
    if num_rotors > 1
        for i in 2:num_rotors
            generate_blade_elements!(
                blade_elements[i],
                rotor_panels[i],
                rotor_parameters[i].xpos,
                rotor_parameters[i].radial_positions,
                rotor_parameters[i].chords,
                rotor_parameters[i].twists,
                rotor_parameters[i].airfoils,
                rotor_parameters[i].num_blade_elements,
                rotor_parameters[i].num_blades,
                rotor_parameters[i].omega,
                body_geometry;
                method=method,
                updated_radial_positions=r_grid_points[:, rotoridxs[i]],
            )

            # - Generate Rotor Panels - #
            ff.generate_panels!(
                method,
                rotor_panels[i],
                [rotor_parameters[i].rotor_x_position .*
                 ones(blade_elements[1].num_radial_stations) r_grid_points[:, rotoridxs[i]]],
            )
        end
    end

    #---------------------------------#
    #             Meshes              #
    #---------------------------------#
    # - Body -> Body - #
    mesh_body_to_body = mesh_body_to_body

    # - Body -> Rotor - #
    mesh_bodies_to_rotor = [
        generate_one_way_mesh(body_panels, rotor_panels[i]) for i in 1:num_rotors
    ]

    # - Wake -> Body - #
    mesh_wake_to_body = [generate_one_way_mesh(wake_panels[i], body_panels) for i in 1:nr]

    # - Wake -> Rotor - #
    mesh_wake_to_rotor = [
        generate_one_way_mesh(wake_panels[i], rotor_panels[j]) for i in 1:nr,
        j in 1:num_rotors
    ]

    # - Rotor -> Body - #
    mesh_rotor_to_body = [
        generate_one_way_mesh(rotor_panels[i], body_panels; singularity="source") for
        i in 1:num_rotors
    ]

    # - Rotor -> Rotor - #
    mesh_rotor_to_rotor = [
        generate_one_way_mesh(rotor_panels[i], rotor_panels[j]; singularity="source") for
        i in 1:num_rotors, j in 1:num_rotors
    ]

    #---------------------------------#
    #       Coefficient Matrices      #
    #---------------------------------#
    # - Body -> Body - #
    A_body_to_body = body_system.A
    bc_freestream_to_body = body_system.b .* freestream.Vinf
    body_vortex_strengths = ImplicitAD.implicit_linear(
        A_body_to_body, bc_freestream_to_body
    )

    # - Body -> Rotor - #
    A_bodies_to_rotor = [
        assemble_one_way_coefficient_matrix(
            mesh_bodies_to_rotor[i], body_panels, rotor_panels[i]
        ) for i in 1:num_rotors
    ]

    # - Wake -> Body - #
    A_wake_to_bodies = [
        assemble_one_way_coefficient_matrix(
            mesh_wake_to_body[i], wake_panels[i], body_panels
        ) for i in 1:nr
    ]

    # - Wake -> Rotor - #
    A_wake_to_rotor = [
        assemble_one_way_coefficient_matrix(
            mesh_wake_to_rotor[i, j], wake_panels[i], rotor_panels[j]
        ) for i in 1:nr, j in 1:num_rotors
    ]

    # - Rotor -> Body - #
    A_rotor_to_body = [
        assemble_one_way_coefficient_matrix(
            mesh_rotor_to_body[i], rotor_panels[i], body_panels; singularity="source"
        ) for i in 1:num_rotors
    ]

    # - Rotor -> Rotor - #
    A_rotor_to_rotor = [
        assemble_one_way_coefficient_matrix(
            mesh_rotor_to_rotor[i, j],
            rotor_panels[i],
            rotor_panels[j];
            singularity="source",
        ) for i in 1:num_rotors, j in 1:num_rotors
    ]

    return (
        # - General - #
        converged=[false], # Initialize Convergence Flag
        body_geometry=body_geometry,
        wake_grid=wake_grid,
        blade_elements=blade_elements,
        rotoridxs=rotoridxs,
        freestream=freestream,
        # - Panels - #
        body_panels=body_panels,
        rotor_panels=rotor_panels,
        wake_panels=wake_panels,
        # - Meshes - #
        mesh_body_to_body=mesh_body_to_body,
        mesh_bodies_to_rotor=mesh_bodies_to_rotor,
        mesh_wake_to_body=mesh_wake_to_body,
        mesh_wake_to_rotor=mesh_wake_to_rotor,
        mesh_wake_to_body=mesh_wake_to_body,
        mesh_wake_to_rotor=mesh_wake_to_rotor,
        mesh_rotor_to_body=mesh_rotor_to_body,
        mesh_rotor_to_rotor=mesh_rotor_to_rotor,
        # - Coefficients - #
        A_body_to_body=A_body_to_body,
        A_A_wake_to_bodyA_wake_to_bodybodies_to_rotor,
        A_wake_to_body=A_wake_to_body,
        A_wake_to_rotor=A_wake_to_rotor,
        A_rotor_to_body=A_rotor_to_body,
        A_rotor_to_rotor=A_rotor_to_rotor,
    )
end

"""
get first guess for gamma and sigma values
"""
function gamma_sigma_guess()
    return nothing
end
