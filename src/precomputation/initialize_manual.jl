
function manual_precomputed_inputs(
    duct_coordinates,      # panel node locations
    hub_coordinates,       # panel node locations
    wake_coordinates,      # tuple containing xgrid[x,r], rgrid[x,r]
    rotor_parameters,      # tuple with rotor paramters
    freestream,            # tuple containing Vinf, rho, mu, asound
    reference_parameters;  # tuple containing Vref, Rref
    debug=false,
)

    ## -- rename for convenience -- ##
    # Save number of wake panels in x-direction
    num_wake_x_panels = length(wake_coordinates.xgrid[1, :]) - 1
    nwakes = length(wake_coordinates.xgrid[:, 1])

    #rotor panel edges
    rpe = rotor_parameters[1].rotor_panel_edges

    #--------------------------------#
    #      Generate Body Panels      #
    #--------------------------------#

    body_panels = generate_body_panels(duct_coordinates, hub_coordinates)

    #--------------------------------#
    #      Generate Wake Panels      #
    #--------------------------------#

    # generate wake sheet paneling
    wake_vortex_panels = generate_wake_panels(
        wake_coordinates.xgrid', wake_coordinates.rgrid'
    )

    # - generate "rotor" panels to receive wake induced velocities - #
    wake_affect_panels = [
        generate_wake_affect_panels(
            (p -> p.panel_center[i, 1]).(wake_vortex_panels),
            (p -> p.panel_center[i, 2]).(wake_vortex_panels),
        ) for i in 1:num_wake_x_panels
    ]

    #------------------------------------------#
    # Generate Rotor Panels and Blade Elements #
    #------------------------------------------#

    # rotor source panel objects
    rotor_source_panels = [
        generate_rotor_panels(rotor_parameters[1].rotorzloc, rpe) for i in 1:1
    ]

    # rotor blade element objects
    blade_elements = [(;
        B=rotor_parameters[1].B,
        Omega=rotor_parameters[1].Omega,
        rotorzloc=rotor_parameters[1].rotorzloc,
        rbe=rotor_source_panels[1].panel_center[:, 2],
        chords=rotor_parameters[1].chords,
        twists=rotor_parameters[1].twists,
        stagger=rotor_parameters[1].stagger,
        solidity=rotor_parameters[1].solidity,
        outer_airfoil=rotor_parameters[1].outer_airfoil,
        inner_airfoil=rotor_parameters[1].inner_airfoil,
        inner_fraction=rotor_parameters[1].inner_fraction,
        Rtip=rotor_parameters[1].Rtip,
        Rhub=rotor_parameters[1].Rhub,
    )]

    #--------------------------------#
    # Generate Relational Geometries #
    #--------------------------------#

    # - Relative to Body - #
    # body to body
    # body_meshff = ff.generate_mesh(body_method, body_panels)
    mesh_bb = generate_one_way_mesh(body_panels, body_panels)

    # body to rotor
    mesh_rb = [
        generate_one_way_mesh(body_panels, rotor_source_panels[i]) for i in 1:1, j in 1:1
    ]

    # rotor to body
    mesh_br = [
        generate_one_way_mesh(rotor_source_panels[j], body_panels) for i in 1:1, j in 1:1
    ]

    # - Relative to Rotor - #
    # rotor to rotor
    # note: broadcasting like this throws an error. so using comprehension instead
    # mesh_rr = generate_one_way_mesh.(rotor_source_panels, rotor_source_panels')
    mesh_rr = [
        generate_one_way_mesh(rotor_source_panels[j], rotor_source_panels[i]) for i in 1:1,
        j in 1:1
    ]

    # wake to body
    mesh_bw = [
        generate_one_way_mesh(wake_vortex_panels[j], body_panels) for i in 1:1,
        j in 1:length(wake_vortex_panels)
    ]

    # wake to rotor
    # mesh_rw = generate_one_way_mesh.(wake_vortex_panels, rotor_source_panels')
    mesh_rw = [
        generate_one_way_mesh(wake_vortex_panels[j], rotor_source_panels[i]) for i in 1:1,
        j in 1:length(wake_vortex_panels)
    ]

    # - Relative to Wake - #
    # body to wake

    mesh_wba = [
        generate_one_way_mesh(body_panels, wake_affect_panels[i]) for
        i in 1:length(wake_affect_panels), j in 1:1
    ]

    mesh_wb = [
        generate_one_way_mesh(body_panels, wake_vortex_panels[i]) for
        i in 1:length(wake_vortex_panels), j in 1:1
    ]

    # rotor to wake
    # mesh_rw = generate_one_way_mesh.(wake_vortex_panels, rotor_source_panels')
    mesh_wr = [
        generate_one_way_mesh(rotor_source_panels[j], wake_vortex_panels[i]) for
        i in 1:length(wake_vortex_panels), j in 1:1
    ]

    # wake to wake
    mesh_ww = [
        generate_one_way_mesh(wake_vortex_panels[j], wake_affect_panels[i]) for
        i in 1:length(wake_affect_panels), j in 1:length(wake_vortex_panels)
    ]

    #---------------------------------#
    # Calculate Coefficient Matrices  #
    #---------------------------------#

    ##### ----- Induced Velcocities on Bodies ----- #####

    # - body to body - #
    # A_bbff = ff.assemble_ring_vortex_matrix_raw(ff.Constant(), [false; true], body_panels, body_meshff)
    A_bb = assemble_induced_velocity_on_body_matrix(
        mesh_bb, body_panels, body_panels; singularity="vortex"
    )

    # apply back-diagonal correction to duct portions of coefficient matrix
    apply_back_diagonal_correction!(
        A_bb, body_panels[1], mesh_bb.affect_panel_indices[1], mesh_bb.mesh2panel_affect
    )

    # - freestream to body - #
    # b_bfff = ff.assemble_ring_boundary_conditions_raw(
    #     ff.Dirichlet(), [false; true], body_panels, body_meshff
    # )
    b_bf =
        freestream.Vinf .*
        assemble_body_freestream_boundary_conditions(body_panels, mesh_bb)

    # - rotor to body - #
    A_br = [
        assemble_induced_velocity_on_body_matrix(
            mesh_br[i, j], [rotor_source_panels[j]], body_panels; singularity="source"
        ) for i in 1:1, j in 1:1
    ]

    # - wake to body - #
    A_bw = [
        assemble_induced_velocity_on_body_matrix(
            mesh_bw[i, j],
            [wake_vortex_panels[j]],
            body_panels;
            singularity="vortex",
            debug=debug,
        ) for i in 1:1, j in 1:length(wake_vortex_panels)
    ]

    ##### ----- Induced Velcocities on Rotors ----- #####
    # - rotor to rotor - #
    A_rr = [
        assemble_induced_velocity_matrices(
            mesh_rr[i, j], rotor_source_panels, rotor_source_panels; singularity="source"
        ) for i in 1:1, j in 1:1
    ]

    # axial components
    vx_rr = [A_rr[i, j][1] for i in 1:1, j in 1:1]

    # radial components
    vr_rr = [A_rr[i, j][2] for i in 1:1, j in 1:1]

    # - body to rotor - #
    A_rb = [
        assemble_induced_velocity_matrices(mesh_rb[i, j], body_panels, rotor_source_panels)
        for i in 1:1, j in 1:1
    ]

    # axial components
    vx_rb = [A_rb[i, j][1] for i in 1:1, j in 1:1]

    # radial components
    vr_rb = [A_rb[i, j][2] for i in 1:1, j in 1:1]

    # - wake to rotor - #
    A_rw = [
        assemble_induced_velocity_matrices(
            mesh_rw[i, j], wake_vortex_panels[j], rotor_source_panels
        ) for i in 1:1, j in 1:length(wake_vortex_panels)
    ]

    # axial components
    vx_rw = [A_rw[i, j][1] for i in 1:1, j in 1:length(wake_vortex_panels)]

    # radial components
    vr_rw = [A_rw[i, j][2] for i in 1:1, j in 1:length(wake_vortex_panels)]

    ##### ----- Induced Velcocities on Wake ----- #####
    # - body to wake - #

    A_wba = [
        assemble_induced_velocity_matrices(
            mesh_wba[i, j], body_panels, wake_affect_panels[i]
        ) for i in 1:length(wake_affect_panels), j in 1:1
    ]

    A_wb = [
        assemble_induced_velocity_matrices(
            mesh_wb[i, j], body_panels, wake_vortex_panels[i]
        ) for i in 1:length(wake_vortex_panels), j in 1:1
    ]

    # axial components

    vx_wba = [A_wba[i, j][1] for i in 1:length(wake_affect_panels), j in 1:1]

    vx_wb = [A_wb[i, j][1] for i in 1:length(wake_vortex_panels), j in 1:1]

    # radial components

    vr_wba = [A_wba[i, j][2] for i in 1:length(wake_affect_panels), j in 1:1]

    vr_wb = [A_wb[i, j][2] for i in 1:length(wake_vortex_panels), j in 1:1]

    # - rotor to wake - #
    A_wr = [
        assemble_induced_velocity_matrices(
            mesh_wr[i, j], rotor_source_panels, wake_vortex_panels[i]
        ) for i in 1:length(wake_vortex_panels), j in 1:1
    ]

    # axial components
    vx_wr = [A_wr[i, j][1] for i in 1:length(wake_vortex_panels), j in 1:1]

    # radial components
    vr_wr = [A_wr[i, j][2] for i in 1:length(wake_vortex_panels), j in 1:1]

    # - wake to wake - #
    A_ww = [
        assemble_induced_velocity_matrices(
            mesh_ww[i, j], wake_vortex_panels[j], wake_affect_panels[i]
        ) for i in 1:length(wake_affect_panels), j in 1:length(wake_vortex_panels)
    ]

    # axial components
    vx_ww = [
        A_ww[i, j][1] for i in 1:length(wake_affect_panels),
        j in 1:length(wake_vortex_panels)
    ]

    # radial components
    vr_ww = [
        A_ww[i, j][2] for i in 1:length(wake_affect_panels),
        j in 1:length(wake_vortex_panels)
    ]
    ## -- Miscellaneous Values for Indexing -- ##

    # - Get rotor panel edges and centers - #
    rotor_panel_edges = [wake_coordinates.rgrid[1:length(rpe), 1] for i in 1:1]
    rotor_panel_centers = [rotor_source_panels[i].panel_center[:, 2] for i in 1:1]

    rotor_panel_edges = reduce(hcat, rotor_panel_edges)
    rotor_panel_centers = reduce(hcat, rotor_panel_centers)

    ## -- Final indexing/bookkeeping -- ##
    # get the total number of vortex panels on the bodies
    num_body_panels = length(b_bf)

    #wake index of duct TE location
    ductTE_index = findfirst(
        x -> x >= maximum(duct_coordinates[:, 1]), wake_coordinates.xgrid[end, :]
    )
    #wake index of hub TE location
    hubTE_index = findfirst(
        x -> x >= maximum(hub_coordinates[:, 1]), wake_coordinates.xgrid[1, :]
    )
    # index of duct panel where rotor lies
    rotor_indices_on_duct =
        findfirst(x -> x < rotor_parameters[1].rotorzloc, body_panels[1].panel_center[:, 1]) -
        1
    # indices on duct after rotor (where deltacp is applied)
    ductwakeidx = [1:rotor_indices_on_duct[end]]
    # index of hub panel where rotor lies
    rotor_indices_on_hub = findfirst(
        x -> x > rotor_parameters[1].rotorzloc, body_panels[2].panel_center[:, 1]
    )
    # indices on hub aft of rotor
    hubwakeidx = [rotor_indices_on_hub[1]:length(body_panels[2].panel_center[:, 1])]

    return (;
        converged=[false],
        #freestream
        freestream,
        #reference for post process
        reference_parameters,
        # body_geometry, # body geometry
        # - rotors
        blade_elements, # blade elements
        num_rotors=1,
        rotor_panel_edges,
        rotor_panel_centers,
        # panels
        rotor_indices_in_wake=[1],
        rotor_indices_on_duct,
        rotor_indices_on_hub,
        num_wake_x_panels=length(wake_vortex_panels[1].panel_center[:, 1]),
        num_body_panels,
        ductTE_index,
        hubTE_index,
        ductwakeidx=ductwakeidx[1],
        hubwakeidx=hubwakeidx[1],
        # - unit induced velocities (INCLUDING PANEL LENGTH)
        A_bb, # body to body
        b_bf, # freestream contribution to body boundary conditions
        A_br, # rotor to body (total)
        A_bw, # wake to body (total)
        A_rb, # body to rotor
        A_wb, # body to wake
        A_wba, # body to "wake"
        A_wr, # rotor to wake
        A_ww, # rotor to wake
        vx_rb, # body to rotor (x-direction)
        vr_rb, # body to rotor (r-direction)
        vx_rr, # rotor to rotor (x-direction)
        vr_rr, # rotor to rotor ( r-direction)
        vx_rw, # wake to rotor (x-direction)
        vr_rw, # wake to rotor ( r-direction)
        vx_wb, # body to wake (x-direction)
        vr_wb, # body to wake ( r-direction)
        vx_wba, # body to "wake" (x-direction)
        vr_wba, # body to "wake" ( r-direction)
        vx_wr, # rotor to wake (x-direction)
        vr_wr, # rotor to wake ( r-direction)
        vx_ww, # wake to wake (x-direction)
        vr_ww, # wake to wake ( r-direction)
        # kutta condition indices
        kutta_idxs=get_kutta_indices([false; true], mesh_bb),
        # operating conditions
        Vinf=freestream.Vinf, # freestream parameters
        # - Debugging/Plotting
        duct_coordinates,
        hub_coordinates,
        isduct=true,
        ishub=true,
        body_panels,
        rotor_source_panels,
        wakexgrid=wake_coordinates.xgrid,
        wakergrid=wake_coordinates.rgrid,
        wake_vortex_panels,
        wake_affect_panels,
        mesh_bb,
        mesh_rb,
        mesh_br,
        mesh_bw,
        mesh_rw,
    )
end
