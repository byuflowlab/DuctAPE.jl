
"""
    initialize_parameters(duct_coordinates, hub_coordinates, rotor_parameters, freestream)

Initializes parameters which are used throughout the analysis.
"""
function initialize_parameters(
    duct_coordinates, hub_coordinates, rotor_parameters, freestream
)
    # number of rotors
    num_rotors = length(rotor_parameters)

    #----------------------------------#
    #             Problem              #
    #----------------------------------#

    # These are the structs used by FLOWFoil for dispatch
    body_method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true])
    vortex_method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true])
    source_method = ff.AxisymmetricProblem(Source(Constant()), Neumann(), [true])

    #---------------------------------#
    #             Bodies              #
    #---------------------------------#

    # Uses FLOWFoil to generate the panel objects
    body_geometry, body_panels = generate_body_geometry(duct_coordinates, hub_coordinates)

    # Use FLOWFoil to generate the body to body mesh
    body_mesh = ff.generate_mesh(body_method, body_panels)

    # Use FLOWFoil to generate the linear system for the no-rotor solution
    body_system = ff.generate_inviscid_system(body_method, body_panels, body_mesh)

    #---------------------------------#
    #           First Rotor           #
    #---------------------------------#

    # Generate blade element structs and initialize all rotors using the first rotor
    blade_elements = initialize_blade_elements(rotor_parameters[1], body_geometry, num_rotors)

    # Generate rotor panels and initialize all rotors using the first rotor
    rotor_panels = initialize_rotor_panels(blade_elements[1], num_rotors; method=source_method)
    dummy_rotor_panels = initialize_dummy_rotor_panels(blade_elements[1], num_rotors; method=source_method)

    #---------------------------------#
    #            Wake Points          #
    #---------------------------------#

    # Generate a streamline-aligned wake grid
    wake_panels, num_wakes, rotoridxs, r_grid_points = generate_wake_grid(
        body_geometry,
        (p -> p.rotor_x_position).(rotor_parameters),
        blade_elements[1].radial_positions;
        method=vortex_method,
        wake_length=1.0,
        debug=false,
    )

    #---------------------------------#
    #          Other Rotor(s)         #
    #---------------------------------#

    # Update aft rotors. The radial positions of these rotors are defined by the radial
    # positions of the wake panels.

    if num_rotors > 1

        # loop through each aft rotor
        for i in 2:num_rotors

            # update aft rotor geometries
            generate_blade_elements!(
                blade_elements[i],
                rotor_parameters[i].rotor_x_position,
                rotor_parameters[i].radial_positions,
                rotor_parameters[i].chords,
                rotor_parameters[i].twists,
                rotor_parameters[i].airfoils,
                nothing, #not actually used when updated radial positions are given
                rotor_parameters[i].num_blades,
                rotor_parameters[i].omega,
                body_geometry;
                updated_radial_positions=r_grid_points[rotoridxs[i], :],
            )

            x_grid_points = fill(rotor_parameters[i].rotor_x_position, blade_elements[1].num_radial_stations[1])

            # update aft rotor panels
            ff.generate_panels!(
                source_method,
                rotor_panels[i],
                [x_grid_points, r_grid_points[rotoridxs[i], :]],
            )

            update_dummy_rotor_panels!(
                blade_elements[i],
                dummy_rotor_panels[i],
                r_grid_points[rotoridxs[i], :];
                method=source_method,
            )
        end
    end

    #---------------------------------#
    #             Meshes              #
    #---------------------------------#

    # For cases other than the body-on-body case, a slightly different mesh implementation
    # needs to be created to allow one-way interactions since all the other coefficient
    # matrices to be created end up on the right hand side of the linear system.

    # - Body -> Rotor - #
    # These meshes are from the body to the rotor blade elements and are used to generate
    # the coefficient matrices that are part of finding the body-induced velocity at the
    # blade elements.
    mesh_body_to_rotor = [
        generate_one_way_mesh(body_panels, dummy_rotor_panels[i]) for i in 1:num_rotors
    ]

    # - Body -> Wake - #
    # These meshes are from the body to the wake elements and are used to generate the
    # coefficient matrices that are part of finding the average velocity in the wakes.
    mesh_body_to_wake = [
        generate_one_way_mesh(body_panels, wake_panels[i]) for i in 1:num_wakes
    ]

    # # - Rotor -> Body - #
    # # For coefficients used in linear solve for body vortex strength residual
    # mesh_rb = [
    #     generate_one_way_mesh(rotor_panels[i], body_panels; singularity="source") for
    #     i in 1:num_rotors
    # ]

    # # - Rotor -> Rotor - #
    # # Meshes from rotor panel centers to blade element locations to account for
    # # rotor-induced velocities, meaning the source panel induced velocities
    # mesh_rr = [
    #     generate_one_way_mesh(rotor_panels[i], dummy_rotor_panels[j]; singularity="source")
    #     for i in 1:num_rotors, j in 1:num_rotors
    # ]

    # # - Rotor -> Wake - #
    # # rotor source panel to wake panels for use in calculating induced velocity on wake
    # mesh_rw = [
    #     generate_one_way_mesh(rotor_panels[i], wake_panels[j]; singularity="source") for
    #     i in 1:num_rotors, j in 1:num_wakes
    # ]

    # - Wake -> Body - #
    # used in generating coefficients used in linear solve for body strength residuals
    mesh_wb = [
        generate_one_way_mesh(wake_panels[i], body_panels) for i in 1:num_wakes
    ]

    # - Wake -> Rotor - #
    # used in generaing coefficients for induced velocity on blade sections
    mesh_wr = [
        generate_one_way_mesh(wake_panels[i], dummy_rotor_panels[j]) for i in 1:num_wakes,
        j in 1:num_rotors
    ]

    # - Wake -> Wake - #
    # self induction of wakes for finding wake velocities
    mesh_ww = [
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
    A_bb = body_system.A

    # Get Right Hand Side vector for the body-freestream interactions
    bc_freestream_to_body = body_system.b .* freestream.Vinf

    # We put the remaining matrices into vectors for easy accessiblity when we don't know beforehand how many rotors will be used.  This allows us to use a single rotor, or a rotor + stator, or as many rotors as we want without having to change implementations or make any declarations ahead of time.
    # - Body -> Rotor - #
    A_body_to_rotor = [
        assemble_one_way_coefficient_matrix(
            mesh_body_to_rotor[i], body_panels, dummy_rotor_panels[i]
        ) for i in 1:num_rotors
    ]
    Ax_br = [A_body_to_rotor[i][1] for i in 1:num_rotors]
    Ar_body_to_rotor = [A_body_to_rotor[i][2] for i in 1:num_rotors]

    # - Body -> Wake - #
    # for finding wake velocities
    A_body_to_wake = [
        assemble_one_way_coefficient_matrix(
            mesh_body_to_wake[i], body_panels, wake_panels[i]
        ) for i in 1:num_wakes
    ]
    Ax_body_to_wake = [A_body_to_wake[i][1] for i in 1:num_wakes]
    Ar_body_to_wake = [A_body_to_wake[i][2] for i in 1:num_wakes]

    # - Rotor -> Body - #
    # for linear solve
    A_rb = [
        assemble_one_way_coefficient_matrix(
            mesh_rb[i], rotor_panels[i], body_panels; singularity="source"
        ) for i in 1:num_rotors
    ]
    Ax_rb = [A_rb[i][1] for i in 1:num_rotors]
    Ar_rb = [A_rb[i][2] for i in 1:num_rotors]

    # - Rotor -> Rotor - #
    # for finding blade element velocities
    # NOTE: since the rotor panels are all aligned, there shouldn't be any self induction in the meridional direction.  consider setting the self-induced coefficient matrices to zeros
    A_rr = [
        assemble_one_way_coefficient_matrix(
            mesh_rr[i, j],
            rotor_panels[i],
            dummy_rotor_panels[j];
            singularity="source",
        ) for i in 1:num_rotors, j in 1:num_rotors
    ]
    Ax_rr = [
        A_rr[i, j][1] for i in 1:num_rotors, j in 1:num_rotors
    ]
    Ar_rr = [
        A_rr[i, j][2] for i in 1:num_rotors, j in 1:num_rotors
    ]

    # - Rotor -> Wake - #
    # for finding wake velocities
    A_rw = [
        assemble_one_way_coefficient_matrix(
            mesh_rw[i, j], rotor_panels[i], wake_panels[j]; singularity="source"
        ) for i in 1:num_rotors, j in 1:num_wakes
    ]
    Ax_rw = [A_rw[i, j][1] for i in 1:num_rotors, j in 1:num_wakes]
    Ar_rw = [A_rw[i, j][2] for i in 1:num_rotors, j in 1:num_wakes]

    # - Wake -> Body - #
    # for linear solve
    A_wb = [
        assemble_one_way_coefficient_matrix(
            mesh_wb[i], wake_panels[i], body_panels
        ) for i in 1:num_wakes
    ]
    Ax_wb = [A_wb[i][1] for i in 1:num_wakes]
    Ar_wb = [A_wb[i][2] for i in 1:num_wakes]

    # - Wake -> Rotor - #
    # for finding blade element velocities
    A_wr = [
        assemble_one_way_coefficient_matrix(
            mesh_wr[i, j], wake_panels[i], dummy_rotor_panels[j]
        ) for i in 1:num_wakes, j in 1:num_rotors
    ]
    Ax_wr = [A_wr[i, j][1] for i in 1:num_wakes, j in 1:num_rotors]
    Ar_wr = [A_wr[i, j][2] for i in 1:num_wakes, j in 1:num_rotors]

    # - Wake -> Wake - #
    # for finding wake velocities
    A_ww = [
        assemble_one_way_coefficient_matrix(
            mesh_ww[i, j], wake_panels[i], wake_panels[j]
        ) for i in 1:num_wakes, j in 1:num_wakes
    ]
    Ax_ww = [A_ww[i, j][1] for i in 1:num_wakes, j in 1:num_wakes]
    Ar_ww = [A_ww[i, j][2] for i in 1:num_wakes, j in 1:num_wakes]

    return (
        # - General - #
        converged=[false], # Initialize Convergence Flag updated in non-linear solve
        freestream=freestream,
        # - Body Only - #
        body_geometry=body_geometry,
        bc_freestream_to_body=bc_freestream_to_body,
        num_body_panels=length(bc_freestream_to_body) - 1,
        # - Rotors - #
        num_rotors=num_rotors,
        blade_elements=blade_elements, # Portions of these will be updated in the coupling with GXBeam potentially
        rotoridxs=rotoridxs,
        # - Wakes - #
        # wake_grid=wake_grid,
        num_wakes=num_wakes,
        num_wake_x_panels=wake_panels[1].npanels,
        # - Panels - #
        body_panels=body_panels,
        rotor_panels=rotor_panels,
        wake_panels=wake_panels,
        # - Coefficients - #
        A_bb=A_bb, # bb
        Ax_br=Ax_br, Ar_br=Ar_br, # body to rotor (x-direction, r-direction)
        Ax_bw=Ax_bw, Ar_bw=Ar_bw, # body to wake (x-direction, r-direction)
        # Ax_rb=Ax_rb, Ar_rb=Ar_rb, # rotor to body (x-direction, r-direction)
        # Ax_rr=Ax_rr, Ar_rr=Ar_rr, # rotor to rotor (x-direction, r-direction)
        # Ax_rw=Ax_rw, Ar_rw=Ar_rw, # rotor to wake (x-direction, r-direction)
        Ax_wb=Ax_wb, Ar_wb=Ar_wb, # wake to body (x-direction, r-direction)
        Ax_wr=Ax_wr, Ar_wr=Ar_wr, # wake to rotor (x-direction, r-direction)
        Ax_ww=Ax_ww, Ar_ww=Ar_ww, # wake to wake (x-direction, r-direction)
    )
end
