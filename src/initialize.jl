
"""
"""
function initialize_parameters(
    duct_coordinates, hub_coordinates, rotor_parameters, freestream
)

    #---------------------------------#
    #             Poblem              #
    #---------------------------------#
    #These are the structs used by FLOWFoil to use for dispatch purposes
    method_body = ff.AxisymmetricProblem(Vortex(Constant()), Neumann(), [false, true])
    method_vortex = ff.AxisymmetricProblem(Vortex(Constant()), Neumann(), [true])
    method_source = ff.AxisymmetricProblem(Source(Constant()), Neumann(), [true])

    problem_body = ff.define_problem(
        method_body, (duct_coordinates, hub_coordinates), 0.0, -1.0, -1.0
    )

    #---------------------------------#
    #             Bodies              #
    #---------------------------------#
    #Uses FLOWFoil to generate the panel objects
    #BodyGeometry is predominantly splines of the duct and hub surfaces
    body_geometry, body_panels = generate_body_geometry(duct_coordinates, hub_coordinates)
    #Use inherent FLOWFoil functions to generate the body_to_body mesh and linear system used for the no-rotor solution
    mesh_body_to_body = ff.generate_mesh(method_body, body_panels)
    body_system = ff.generate_inviscid_system(method_body, body_panels, mesh_body_to_body)

    #---------------------------------#
    #           First Rotor           #
    #---------------------------------#
    # - Rename for Convenience - #
    num_rotors = length(rotor_parameters)

    # - Generate blade element structs with blade element distributions as well as other useful rotor parameters - #
    # Note that this initializes all the rotors to the same values as the first rotor. Aft rotors are updated after the wake generation
    blade_elements = initialize_blade_elements(
        rotor_parameters[1], body_geometry, num_rotors
    )

    nbe = blade_elements[1].num_radial_stations[1]

    # - Generate Rotor Panels - #
    # Again, this contains information for all of the rotors, using the first rotor data. Updates to aft rotors take place after the wake generation.
    rotor_panels = initialize_rotor_panels(
        blade_elements[1], num_rotors; method=method_source
    )
    dummy_rotor_panels = initialize_dummy_rotor_panels(
        blade_elements[1], num_rotors; method=method_source
    )

    #---------------------------------#
    #            Wake Points          #
    #---------------------------------#
    #TODO: need to update grid to emminate from center of rotor panels
    #TODO: not sure what this will entail, but could be as simple as changing the starting r-positions.
    #TODO: also output number of wake objects (length of wake panels vector)
    wake_grid, wake_panels, num_wakes, rotoridxs = generate_wake_grid(
        body_geometry,
        (p -> p.rotor_x_position).(rotor_parameters),
        rotor_panels[1].panel_center[:, 2];
        method=method_vortex,
        wake_length=1.0,
        debug=false,
    )

    #---------------------------------#
    #          Other Rotor(s)         #
    #---------------------------------#
    #These need to be updated here, rather than generated before the wake generation because the rotor panels need to align with the wake grid.  The first rotor defines the spacing of the foremost panels, but the grid generation designates the radial positions of any aft rotors.
    if num_rotors > 1
        for i in 2:num_rotors

            # - Update Aft Rotor Geometries - #
            generate_blade_elements!(
                blade_elements[i],
                rotor_parameters[i].xpos,
                rotor_parameters[i].radial_positions,
                rotor_parameters[i].chords,
                rotor_parameters[i].twists,
                rotor_parameters[i].airfoils,
                rotor_parameters[i].num_blade_elements,
                rotor_parameters[i].num_blades,
                rotor_parameters[i].omega,
                body_geometry;
                updated_radial_positions=r_grid_points[:, rotoridxs[i]],
            )

            # - Update Aft Rotor Panels - #
            ff.generate_panels!(
                method_source,
                rotor_panels[i],
                [rotor_parameters[i].rotor_x_position .*
                 ones(blade_elements[1].num_radial_stations) r_grid_points[:, rotoridxs[i]]],
            )

            update_dummy_rotor_panels!(
                blade_elements[i],
                dummy_rotor_panels[i],
                r_grid_points[:, rotoridxs[i]];
                method=method_source,
            )
        end
    end

    #---------------------------------#
    #             Meshes              #
    #---------------------------------#
    #For cases other than the body-on-body case, a slightly different mesh implementation needed to be created so as to allow one-way interactions since all the other coefficient matrices to be created end up on the right hand side of the linear system.

    # - Body -> Rotor - #
    mesh_body_to_rotor = [
        generate_one_way_mesh(body_panels, dummy_rotor_panels[i]) for i in 1:num_rotors
    ]

    # - Body -> Wake - #
    mesh_body_to_wake = [
        generate_one_way_mesh(body_panels, wake_panels[i]) for i in 1:num_wakes
    ]

    # - Rotor -> Body - #
    mesh_rotor_to_body = [
        generate_one_way_mesh(rotor_panels[i], body_panels; singularity="source") for
        i in 1:num_rotors
    ]

    # - Rotor -> Rotor - #
    mesh_rotor_to_rotor = [
        generate_one_way_mesh(rotor_panels[i], dummy_rotor_panels[j]; singularity="source")
        for i in 1:num_rotors, j in 1:num_rotors
    ]

    # - rotor -> Wake - #
    mesh_rotor_to_wake = [
        generate_one_way_mesh(rotor_panels[i], wake_panels[j]) for i in 1:num_rotors,
        j in 1:num_wakes
    ]

    # - Wake -> Body - #
    mesh_wake_to_body = [
        generate_one_way_mesh(wake_panels[i], body_panels) for i in 1:num_wakes
    ]

    # - Wake -> Rotor - #
    mesh_wake_to_rotor = [
        generate_one_way_mesh(wake_panels[i], dummy_rotor_panels[j]) for i in 1:num_wakes,
        j in 1:num_rotors
    ]

    # - Wake -> Wake - #
    mesh_wake_to_wake = [
        generate_one_way_mesh(wake_panels[i], wake_panels[j]) for i in 1:num_wakes,
        j in 1:num_wakes
    ]

    #---------------------------------#
    #       Coefficient Matrices      #
    #---------------------------------#
    # - Body -> Body - #
    # The body-to-body case is simply the default FLOWFoil functionality.
    # Note that we apply the freestream velocity here, as it will be needed to be a part of the system throughout rather than just being applied in post-processing as is done in FLOWFoil.

    # Get Left Hand Side Matrix
    A_body_to_body = body_system.A

    # Get Right Hand Side vector for the body-freestream interactions
    bc_freestream_to_body = body_system.b .* freestream.Vinf

    # We put the remaining matrices into vectors for easy accessiblity when we don't know beforehand how many rotors will be used.  This allows us to use a single rotor, or a rotor + stator, or as many rotors as we want without having to change implementations or make any declarations ahead of time.
    # - Body -> Rotor - #
    A_body_to_rotor = [
        assemble_one_way_coefficient_matrix(
            mesh_body_to_rotor[i], body_panels, dummy_rotor_panels[i]
        ) for i in 1:num_rotors
    ]

    A_body_to_wake = [
        assemble_one_way_coefficient_matrix(
            mesh_body_to_wake[i], body_panels, wake_panels[i]
        ) for i in 1:num_wakes
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
            dummy_rotor_panels[j];
            singularity="source",
        ) for i in 1:num_rotors, j in 1:num_rotors
    ]

    # - Rotor -> Wake - #
    A_rotor_to_wake = [
        assemble_one_way_coefficient_matrix(
            mesh_rotor_to_wake[i, j], rotor_panels[i], wake_panels[j]; singularity="source"
        ) for i in 1:num_rotors, j in 1:num_wakes
    ]

    # - Wake -> Body - #
    A_wake_to_body = [
        assemble_one_way_coefficient_matrix(
            mesh_wake_to_body[i], wake_panels[i], body_panels
        ) for i in 1:num_wakes
    ]

    # - Wake -> Rotor - #
    A_wake_to_rotor = [
        assemble_one_way_coefficient_matrix(
            mesh_wake_to_rotor[i, j], wake_panels[i], dummy_rotor_panels[j]
        ) for i in 1:num_wakes, j in 1:num_rotors
    ]

    # - Wake -> Wake - #
    A_wake_to_wake = [
        assemble_one_way_coefficient_matrix(
            mesh_wake_to_wake[i, j], wake_panels[i], dummy_wake_panels[j]
        ) for i in 1:num_wakes, j in 1:num_wakes
    ]

    #TODO: return everything in a named tuple so that order doesn't matter during development.  May want to make this a struct later, but probably not.  Need to consult with others who have more experience to see if/what changes need to be made here.
    return (
        # - General - #
        converged=[false], # Initialize Convergence Flag updated in non-linear solve
        freestream=freestream,
        # - Body Only - #
        body_geometry=body_geometry,
        bc_freestream_to_body=bc_freestream_to_body,
        num_body_panels=length(bc_freestream_to_body) - 1,
        wake_grid=wake_grid,
        # - Rotors - #
        num_rotors=num_rotors,
        blade_elements=blade_elements, # Portions of these will be updated in the coupling with GXBeam potentially
        rotoridxs=rotoridxs,
        # - Panels - #
        body_panels=body_panels,
        rotor_panels=rotor_panels,
        wake_panels=wake_panels,
        # - Meshes - #
        # TODO: not sure if these are needed to output or not.  Probably clean them out later if not.
        # mesh_body_to_body=mesh_body_to_body,
        # mesh_body_to_rotor=mesh_body_to_rotor,
        # mesh_wake_to_body=mesh_wake_to_body,
        # mesh_wake_to_rotor=mesh_wake_to_rotor,
        # mesh_rotor_to_body=mesh_rotor_to_body,
        # mesh_rotor_to_rotor=mesh_rotor_to_rotor,
        # - Coefficients - #
        A_body_to_body=A_body_to_body,
        A_body_to_rotor=A_body_to_rotor,
        A_body_to_wake=A_body_to_wake,
        A_rotor_to_body=A_rotor_to_body,
        A_rotor_to_rotor=A_rotor_to_rotor,
        A_rotor_to_wake=A_rotor_to_wake,
        A_wake_to_body=A_wake_to_body,
        A_wake_to_rotor=A_wake_to_rotor,
        A_wake_to_wake=A_wake_to_wake,
    )
end
