
"""
    initialize_parameters(duct_coordinates, hub_coordinates, rotor_parameters, freestream)

Initializes parameters which are used throughout the analysis.
"""
function initialize_parameters(
    duct_coordinates, hub_coordinates, rotor_parameters, freestream
)
    # number of rotors
    nrotor = length(rotor_parameters)

    # --- Methods --- #

    # method for the body to body problem
    body_method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true])

    # method for vortex panels
    vortex_method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true])
    
    # method for source panels
    source_method = ff.AxisymmetricProblem(Source(Constant()), Neumann(), [true])

    # --- Bodies --- #

    # generate duct and hub geometry and paneling
    body_geometry, body_panels = generate_body_geometry(duct_coordinates, hub_coordinates)

    # --- First Rotor --- #

    # initialize blade elements for first rotor
    first_blade_elements = generate_blade_elements(
        rotor_parameters.B,
        rotor_parameters.omega,
        rotor_parameters.xrotor,
        rotor_parameters.rblade,
        rotor_parameters.chords,
        rotor_parameters.twists,
        rotor_parameters.airfoils,
        body_geometry,
        rotor_parameters.nr)

    # initialize panels for first rotor
    first_rotor_panels = generate_rotor_panels(first_blade_elements, source_method)
    first_dummy_rotor_panels = generate_dummy_rotor_panels(first_blade_elements, source_method)

    # --- Wake Grid --- #

    # Generate a streamline-aligned wake grid
    wake_panels, num_wakes, rotoridxs, r_grid_points = generate_wake_grid(
        body_geometry,
        (p -> p.rotor_x_position).(rotor_parameters),
        blade_elements[1].radial_positions;
        method=vortex_method,
        wake_length=1.0,
        debug=false,
    )

    # --- Other Rotors --- #

    # allocate space for all rotors
    blade_elements = fill(first_blade_elements, nrotor)
    rotor_panels = fill(first_rotor_panels, nrotor)
    dummy_rotor_panels = fill(first_dummy_rotor_panels, nrotor)

    # calculate properties for remaining rotors
    for i in 2:nrotor

        # update aft rotor geometries
        blade_elements[i] = generate_blade_elements(
            rotor_parameters[i].B,
            rotor_parameters[i].omega,
            rotor_parameters[i].xrotor,
            rotor_parameters[i].rblade,
            rotor_parameters[i].chords,
            rotor_parameters[i].twists,
            rotor_parameters[i].airfoils,
            body_geometry,
            size(r_grid_points, 2), # number of grid points matches wake
            view(r_grid_points, rotor_indices[i], :) # use wake coordinates for rotor
        )

        # update aft rotor panels
        rotor_panels[i] = generate_rotor_panels(blade_elements[i], source_method)
        dummy_rotor_panels[i] = generate_dummy_rotor_panels(blade_elements[i], source_method)

    end

    # --- Meshes --- #

    # body to body
    mesh_bb = ff.generate_mesh(body_method, body_panels)

    # body to rotor
    mesh_br = generate_one_way_mesh.(Ref(body_panels), dummy_rotor_panels)

    # body to wake
    mesh_bw = generate_one_way_mesh.(Ref(body_panels), wake_panels)

    # rotor to body
    mesh_rb = generate_one_way_mesh.(rotor_panels, Ref(body_panels); singularity="source")

    # rotor to rotor
    mesh_rr = generate_one_way_mesh.(rotor_panels, dummy_rotor_panels'; singularity="source")

    # rotor to wake
    mesh_rw = generate_one_way_mesh.(rotor_panels, wake_panels'; singularity="source")

    # wake to body
    mesh_wb = generate_one_way_mesh.(wake_panels, Ref(body_panels))

    # wake to rotor
    mesh_wr = generate_one_way_mesh.(wake_panels, dummy_rotor_panels')

    # wake to wake
    mesh_ww = generate_one_way_mesh.(wake_panels, wake_panels')

    # --- Coefficient Matrices --- #
    
    body_system = ff.generate_inviscid_system(body_method, body_panels, body_mesh)

    # body to body
    A_bb = body_system.A

    # freestream to body
    b_fb = body_system.b .* freestream.Vinf

    # body to rotor
    A_br = assemble_one_way_coefficient_matrix.(mesh_br, Ref(body_panels), dummy_rotor_panels)
    Ax_br = getindex.(A_br, 1)
    Ar_br = getindex.(A_br, 2)

    # body to wake
    A_bw = assemble_one_way_coefficient_matrix.(mesh_bw, Ref(body_panels), wake_panels)
    Ax_bw = getindex.(A_bw, 1)
    Ar_bw = getindex.(A_bw, 2)

    # rotor to body
    A_rb = assemble_one_way_coefficient_matrix.(mesh_rb, rotor_panels, Ref(body_panels); singularity="source")
    Ax_rb = getindex.(A_rb, 1)
    Ar_rb = getindex.(A_rb, 2)

    # rotor to rotor
    A_rr = assemble_one_way_coefficient_matrix.(mesh_rr, rotor_panels, dummy_rotor_panels'; singularity="source")
    Ax_rr = getindex.(A_rr, 1)
    Ar_rr = getindex.(A_rr, 2)

    # rotor to wake
    A_rw = assemble_one_way_coefficient_matrix.(mesh_rw, rotor_panels, wake_panels'; singularity="source")
    Ax_rw = getindex.(A_rw, 1)
    Ar_rw = getindex.(A_rw, 2)

    # wake to body
    A_wb = assemble_one_way_coefficient_matrix.(mesh_wb, wake_panels, Ref(body_panels))
    Ax_wb = getindex.(A_wb, 1)
    Ar_wb = getindex.(A_wb, 2)

    # wake to rotor
    A_wr = assemble_one_way_coefficient_matrix.(mesh_wr, wake_panels, dummy_rotor_panels')
    Ax_wr = getindex.(A_wr, 1)
    Ar_wr = getindex.(A_wb, 2)

    # wake to wake
    A_ww = assemble_one_way_coefficient_matrix.(mesh_ww, wake_panels, wake_panels')
    Ax_ww = getindex.(A_ww, 1)
    Ar_ww = getindex.(A_ww, 2)


    return (
        # body
        body_geometry = body_geometry, # body geometry
        # rotors
        blade_elements = blade_elements, # blade elements
        rotor_indices = rotor_indices, # rotor locations
        # panels
        body_panels = body_panels, # body paneling
        rotor_panels = rotor_panels, # rotor paneling
        wake_panels = wake_panels, # wake paneling
        # aerodynamic influence coefficients
        A_bb = A_bb, # body to body
        b_fb = b_fb, # freestream contribution to body boundary conditions
        Ax_br = Ax_br, Ar_br = Ar_br, # body to rotor (x-direction, r-direction)
        Ax_bw = Ax_bw, Ar_bw = Ar_bw, # body to wake (x-direction, r-direction)
        Ax_rb = Ax_rb, Ar_rb = Ar_rb, # rotor to body (x-direction, r-direction)
        Ax_rr = Ax_rr, Ar_rr = Ar_rr, # rotor to rotor (x-direction, r-direction)
        Ax_rw = Ax_rw, Ar_rw = Ar_rw, # rotor to wake (x-direction, r-direction)
        Ax_wb = Ax_wb, Ar_wb = Ar_wb, # wake to body (x-direction, r-direction)
        Ax_wr = Ax_wr, Ar_wr = Ar_wr, # wake to rotor (x-direction, r-direction)
        Ax_ww = Ax_ww, Ar_ww = Ar_ww, # wake to wake (x-direction, r-direction)
        # operating conditions
        freestream = freestream, # freestream parameters
    )
end
