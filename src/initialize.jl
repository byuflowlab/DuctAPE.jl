"""
    precomputed_inputs(duct_coordinates, hub_coordinates, rotor_parameters, freestream; kwargs...)

Initializes the geometry, panels, and aerodynamic influence coefficient matrices.

# Arguments:
- `duct_coordinates::Array{Float64,2}` : duct coordinates, starting from trailing edge, going clockwise
- `hub_coordinates::Array{Float64,2}` : hub coordinates, starting from leading edge
- `rotor_parameters`: named tuple of rotor parameters
- `freestream`: freestream parameters

# Keyword Arguments
- `paneling_constants.wake_length=1.0` : non-dimensional length (based on maximum duct chord) that the wake
    extends past the furthest trailing edge.
- `nwake_sheets=length(rotor_parameters[1].rblade)`: Number of radial stations to use when defining
    the wake
- `finterp=FLOWMath.akima`: Method used to interpolate the duct and hub geometries.
- `xwake`: May be used to define a custom set of x-coordinates for the wake.
- `nhub`: Number of panels between the hub leading edge and the first rotor.
- `nduct_inner`: Number of panels between the duct leading edge and the first rotor.
- `nduct_outer`: Number of panels on the duct outer surface.

NOTE: The `paneling_constants.nwake_sheets`, `xwake`, `nhub`, `nduct_inner`, and `nduct_outer` arguments should
always be provided when performing gradient-based optimizatoin in order to ensure that the
resulting paneling is consistent between optimization iterations.
"""
function precomputed_inputs(
    duct_coordinates,
    hub_coordinates,
    paneling_constants,
    rotor_parameters, #vector of named tuples
    freestream,
    reference_parameters;
    finterp=fm.akima,
    debug=false,
)

    # ## -- Rename for Convenience -- ##
    # # promoted type
    # TF = promote_type(
    #     eltype(rotor_parameters.chords),
    #     eltype(rotor_parameters.twists),
    #     eltype(rotor_parameters.Omega),
    #     eltype(rotor_parameters.tip_gap),
    #     eltype(duct_coordinates),
    #     eltype(freestream.Vinf),
    # )

    # number of rotors
    nrotor = length(rotor_parameters)

    #------------------------------------#
    # Discretize Wake and Repanel Bodies #
    #------------------------------------#
    rotoronly = false
    nohub = false
    noduct = false
    if hub_coordinates == nothing && duct_coordinates != nothing
        nohub = true
        hub_coordinates = [
            minimum(duct_coordinates[:, 1]) 0.0
            maximum(duct_coordinates[:, 1]) 0.0
        ]
    elseif hub_coordinates != nothing && duct_coordinates == nothing
        noduct = true
        hub_coordinates = [
            minimum(hub_coordinates[:, 1]) 0.0
            maximum(hub_coordinates[:, 1]) 0.0
        ]
    elseif hub_coordinates == nothing && duct_coordinates == nothing
        rotoronly = true
    end

    # - Discretize Wake x-coordinates - #
    # also returns indices of rotor locations in the wake
    xwake, rotor_indices_in_wake, ductTE_index, hubTE_index = discretize_wake(
        duct_coordinates,
        hub_coordinates,
        rotor_parameters.xrotor,
        paneling_constants.wake_length,
        paneling_constants.npanels,
    )

    # Save number of wake panels in x-direction
    num_wake_x_panels = length(xwake) - 1

    # - Repanel Bodies - #
    rp_duct_coordinates, rp_hub_coordinates = update_body_geometry(
        duct_coordinates,
        hub_coordinates,
        xwake,
        paneling_constants.nhub_inlet,
        paneling_constants.nduct_inlet;
        finterp=finterp,
    )

    #-----------------------------------#
    # Position Duct and Get Rotor Radii #
    #-----------------------------------#
    # check that tip gap isn't too small
    # set vector of tip gaps to zero for initialization (any overrides to zero thus taken c of automatically)
    tip_gaps = zeros(eltype(rotor_parameters.tip_gap), nrotor)
    if rotor_parameters[1].tip_gap != 0.0
        if rotor_parameters[1].tip_gap < 1e-4
            @warn "You have selected a tip gap for the foremost rotor that is smaller than 1e-4. Overriding to 0.0 to avoid near singularity issues."
        else
            tip_gaps[1] = rotor_parameters[1].tip_gap
        end
    end

    # can't have non-zero tip gaps for aft rotors
    for ir in 2:nrotor
        if rotor_parameters[ir].tip_gap != 0.0
            @warn "DuctTAPE does not currently have capabilities for adding tip gap to any but the foremost rotor. Overriding to 0.0."
        else
            tip_gaps[ir] = rotor_parameters[ir].tip_gap
        end
    end

    # if hub was nothing, set hub radius to dimensional inner rotor radius
    if nohub
        rp_hub_coordinates[:, 2] .= rotor_parameters[1].r[1] * rotor_parameters[1].Rtip
    end
    t_duct_coordinates, Rtips, Rhubs = place_duct(
        rp_duct_coordinates,
        rp_hub_coordinates,
        rotor_parameters[1].Rtip,
        tip_gaps,
        rotor_parameters.xrotor,
    )

    # generate body paneling
    if nohub
        body_panels = generate_body_panels(t_duct_coordinates, nothing)
    else
        body_panels = generate_body_panels(t_duct_coordinates, rp_hub_coordinates)
    end

    # book keep wake wall interfaces
    if nohub
        rotor_indices_on_hub = nothing
        hubwakeidx = nothing
    else
        rotor_indices_on_hub = [
            findfirst(
                x -> x > rotor_parameters.xrotor[i], body_panels[2].panel_center[:, 1]
            ) for i in 1:length(rotor_parameters.xrotor)
        ]
        hwidraw = sort(
            [length(body_panels[2].panel_center[:, 1]); rotor_indices_on_hub .- 1]
        )
        hubwakeidx = reduce(vcat, [[(hwidraw[i] + 1):hwidraw[i + 1]] for i in 1:nrotor])
    end

    if noduct
        rotor_indices_on_duct = nothing
        ductwakeidx = nothing
    else
        rotor_indices_on_duct = [
            findfirst(
                x -> x < rotor_parameters.xrotor[i], body_panels[1].panel_center[:, 1]
            ) - 1 for i in 1:length(rotor_parameters.xrotor)
        ]
        dwidraw = sort([0; rotor_indices_on_duct])
        ductwakeidx = reduce(vcat, [[(dwidraw[i] + 1):dwidraw[i + 1]] for i in 1:nrotor])
    end

    #----------------------------------#
    # Generate Discretized Wake Sheets #
    #----------------------------------#
    # get discretization of wakes at leading rotor position

    for i in 1:nrotor
        @assert Rtips[i] > Rhubs[i] "Rotor #$i Tip Radius is set to be less than its Hub Radius."
    end

    #rotor panel edges
    rpe = range(Rhubs[1], Rtips[1]; length=paneling_constants.nwake_sheets)

    # wake sheet starting radius including dummy sheets for tip gap.
    if tip_gaps[1] == 0.0
        rwake = rpe
    else
        rwake = [rpe; Rtips[1] + tip_gaps[1]]
    end

    # Initialize wake "grid"
    xgrid, rgrid = initialize_wake_grid(
        t_duct_coordinates, rp_hub_coordinates, xwake, rwake
    )

    # Relax "Grid"
    relax_grid!(xgrid, rgrid; max_iterations=100, tol=1e-9, verbose=false)

    # generate wake sheet paneling
    wake_vortex_panels = generate_wake_panels(
        xgrid[:, 1:length(rpe)], rgrid[:, 1:length(rpe)]
    )

    # generate body wake panels for convenience
    # TODO: move panel generation into DuctTAPE
    duct_wake_panels = ff.generate_panels(
        ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true]),
        [xgrid[ductTE_index:end, end] rgrid[ductTE_index:end, end]],
    )
    hub_wake_panels = ff.generate_panels(
        ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true]),
        [xgrid[hubTE_index:end, 1] rgrid[hubTE_index:end, 1]],
    )

    #TODO: probably remove
    # # - generate "rotor" panels to receive wake induced velocities - #
    # wake_affect_panels = [
    #     generate_wake_affect_panels(
    #         (p -> p.panel_center[i, 1]).(wake_vortex_panels),
    #         (p -> p.panel_center[i, 2]).(wake_vortex_panels),
    #     ) for i in 1:num_wake_x_panels
    # ]

    #------------------------------------------#
    # Generate Rotor Panels and Blade Elements #
    #------------------------------------------#

    # rotor source panel objects
    rotor_source_panels = [
        generate_rotor_panels(
            rotor_parameters[i].xrotor, rgrid[rotor_indices_in_wake[i], 1:length(rpe)]
        ) for i in 1:nrotor
    ]

    # rotor blade element objects
    blade_elements = [
        generate_blade_elements(
            rotor_parameters[i].B,
            rotor_parameters[i].Omega,
            rotor_parameters[i].xrotor,
            rotor_parameters[i].r,
            rotor_parameters[i].chords,
            rotor_parameters[i].twists,
            rotor_parameters[i].airfoils,
            Rtips[i],
            Rhubs[i],
            rotor_source_panels[i].panel_center[:, 2],
        ) for i in 1:nrotor
    ]

    #--------------------------------#
    # Generate Relational Geometries #
    #--------------------------------#

    # - Relative to Body - #
    # body to body
    # body_meshff = ff.generate_mesh(body_method, body_panels)
    mesh_bb = generate_one_way_mesh(body_panels, body_panels)

    # body to rotor
    mesh_rb = [
        generate_one_way_mesh(body_panels, rotor_source_panels[i]) for
        i in 1:length(rotor_source_panels), j in 1:1
    ]

    # rotor to body
    mesh_br = [
        generate_one_way_mesh(rotor_source_panels[j], body_panels) for i in 1:1,
        j in 1:length(rotor_source_panels)
    ]

    # - Relative to Rotor - #
    # rotor to rotor
    # note: broadcasting like this throws an error. so using comprehension instead
    # mesh_rr = generate_one_way_mesh.(rotor_source_panels, rotor_source_panels')
    mesh_rr = [
        generate_one_way_mesh(rotor_source_panels[j], rotor_source_panels[i]) for
        i in 1:length(rotor_source_panels), j in 1:length(rotor_source_panels)
    ]

    # wake to body
    mesh_bw = [
        generate_one_way_mesh(wake_vortex_panels[j], body_panels) for i in 1:1,
        j in 1:length(wake_vortex_panels)
    ]

    # wake to rotor
    # mesh_rw = generate_one_way_mesh.(wake_vortex_panels, rotor_source_panels')
    mesh_rw = [
        generate_one_way_mesh(wake_vortex_panels[j], rotor_source_panels[i]) for
        i in 1:length(rotor_source_panels), j in 1:length(wake_vortex_panels)
    ]

    # - Relative to Wake - #
    # body to wake

    #TODO: probably remove
    # mesh_wba = [
    #     generate_one_way_mesh(body_panels, wake_affect_panels[i]) for
    #     i in 1:length(wake_affect_panels), j in 1:1
    # ]

    mesh_wb = [
        generate_one_way_mesh(body_panels, wake_vortex_panels[i]) for
        i in 1:length(wake_vortex_panels), j in 1:1
    ]

    # rotor to wake
    # mesh_rw = generate_one_way_mesh.(wake_vortex_panels, rotor_source_panels')
    mesh_wr = [
        generate_one_way_mesh(rotor_source_panels[j], wake_vortex_panels[i]) for
        i in 1:length(wake_vortex_panels), j in 1:length(rotor_source_panels)
    ]

    # wake to wake
    mesh_ww = [
        generate_one_way_mesh(wake_vortex_panels[j], wake_vortex_panels[i]) for
        i in 1:length(wake_vortex_panels), j in 1:length(wake_vortex_panels)
    ]

    # - Relative to duct wake - #
    # body to duct wake
    mesh_dwb = [generate_one_way_mesh(body_panels, duct_wake_panels) for i in 1:1, j in 1:1]

    # rotor to duct wake
    mesh_dwr = [
        generate_one_way_mesh(rotor_source_panels[j], duct_wake_panels) for i in 1:1,
        j in 1:length(rotor_source_panels)
    ]

    # rotor wake to duct wake
    mesh_dww = [
        generate_one_way_mesh(wake_vortex_panels[j], duct_wake_panels) for i in 1:1,
        j in 1:length(wake_vortex_panels)
    ]

    # - Relative to hub wake - #
    # body to hub wake
    mesh_hwb = [generate_one_way_mesh(body_panels, hub_wake_panels) for i in 1:1, j in 1:1]

    # rotor to hub wake
    mesh_hwr = [
        generate_one_way_mesh(rotor_source_panels[j], hub_wake_panels) for i in 1:1,
        j in 1:length(rotor_source_panels)
    ]

    # rotor wake to hub wake
    mesh_hww = [
        generate_one_way_mesh(wake_vortex_panels[j], hub_wake_panels) for i in 1:1,
        j in 1:length(wake_vortex_panels)
    ]

    #---------------------------------#
    # Calculate Coefficient Matrices  #
    #---------------------------------#

    ##### ----- Induced Velcocities on Bodies ----- #####

    # - body to body - #
    # A_bbff = ff.assemble_ring_vortex_matrix_raw(ff.Constant(), [false; true], body_panels, body_meshff)
    A_bb_raw = assemble_induced_velocity_on_body_matrix(
        mesh_bb, body_panels, body_panels; singularity="vortex"
    )

    # apply back-diagonal correction to duct portions of coefficient matrix
    apply_back_diagonal_correction!(
        A_bb_raw, body_panels[1], mesh_bb.affect_panel_indices[1], mesh_bb.mesh2panel_affect
    )

    # add kutta condition to coefficient matrix
    kidx = get_kutta_indices(false, mesh_bb)
    A_bb_kr = [
        A_bb_raw[1:kidx[end], :]
        zeros(eltype(A_bb_raw), length(A_bb_raw[1, :]))'
        A_bb_raw[(kidx[end] + 1):end, :]
    ]
    A_bb_kc = [A_bb_kr[:, 1:kidx[end]] ones(eltype(A_bb_kr), length(A_bb_kr[:, 1])) A_bb_kr[
        :, (kidx[end] + 1):end
    ]]
    A_bb_kc[kidx[end] + 1, kidx[end]] = eltype(A_bb_kc)(1.0)
    A_bb_kc[kidx[end] + 1, 1] = eltype(A_bb_kc)(1.0)
    A_bb_kc[kidx[end] + 1, kidx[end] + 1] = eltype(A_bb_kc)(0.0)

    # Factorize Matrix using LinearAlgebra
    A_bb = factorize(A_bb_kc)

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
        ) for i in 1:1, j in 1:length(rotor_source_panels)
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
            mesh_rr[i, j],
            rotor_source_panels[j],
            rotor_source_panels[i];
            singularity="source",
        ) for i in 1:length(rotor_source_panels), j in 1:length(rotor_source_panels)
    ]

    # axial components
    vx_rr = [
        A_rr[i, j][1] for i in 1:length(rotor_source_panels),
        j in 1:length(rotor_source_panels)
    ]

    # radial components
    vr_rr = [
        A_rr[i, j][2] for i in 1:length(rotor_source_panels),
        j in 1:length(rotor_source_panels)
    ]

    # - body to rotor - #
    A_rb = [
        assemble_induced_velocity_matrices(
            mesh_rb[i, j], body_panels, rotor_source_panels[i]
        ) for i in 1:length(rotor_source_panels), j in 1:1
    ]

    # axial components
    vx_rb = [A_rb[i, j][1] for i in 1:length(rotor_source_panels), j in 1:1]

    # radial components
    vr_rb = [A_rb[i, j][2] for i in 1:length(rotor_source_panels), j in 1:1]

    # - wake to rotor - #
    A_rw = [
        assemble_induced_velocity_matrices(
            mesh_rw[i, j], wake_vortex_panels[j], rotor_source_panels[i]
        ) for i in 1:length(rotor_source_panels), j in 1:length(wake_vortex_panels)
    ]

    # axial components
    vx_rw = [
        A_rw[i, j][1] for i in 1:length(rotor_source_panels),
        j in 1:length(wake_vortex_panels)
    ]

    # radial components
    vr_rw = [
        A_rw[i, j][2] for i in 1:length(rotor_source_panels),
        j in 1:length(wake_vortex_panels)
    ]

    ##### ----- Induced Velcocities on Wake ----- #####
    # - body to wake - #

    #TODO: probably remove
    # A_wba = [
    #     assemble_induced_velocity_matrices(
    #         mesh_wba[i, j], body_panels, wake_affect_panels[i]
    #     ) for i in 1:length(wake_affect_panels), j in 1:1
    # ]

    A_wb = [
        assemble_induced_velocity_matrices(
            mesh_wb[i, j], body_panels, wake_vortex_panels[i]
        ) for i in 1:length(wake_vortex_panels), j in 1:1
    ]

    # axial components
    #TODO: probably remove
    # vx_wba = [A_wba[i, j][1] for i in 1:length(wake_affect_panels), j in 1:1]

    vx_wb = [A_wb[i, j][1] for i in 1:length(wake_vortex_panels), j in 1:1]

    # radial components

    #TODO: probably remove
    # vr_wba = [A_wba[i, j][2] for i in 1:length(wake_affect_panels), j in 1:1]

    vr_wb = [A_wb[i, j][2] for i in 1:length(wake_vortex_panels), j in 1:1]

    # - rotor to wake - #
    A_wr = [
        assemble_induced_velocity_matrices(
            mesh_wr[i, j], rotor_source_panels[j], wake_vortex_panels[i]
        ) for i in 1:length(wake_vortex_panels), j in 1:length(rotor_source_panels)
    ]

    # axial components
    vx_wr = [
        A_wr[i, j][1] for i in 1:length(wake_vortex_panels),
        j in 1:length(rotor_source_panels)
    ]

    # radial components
    vr_wr = [
        A_wr[i, j][2] for i in 1:length(wake_vortex_panels),
        j in 1:length(rotor_source_panels)
    ]

    # - wake to wake - #
    A_ww = [
        assemble_induced_velocity_matrices(
            mesh_ww[i, j], wake_vortex_panels[j], wake_vortex_panels[i]
        ) for i in 1:length(wake_vortex_panels), j in 1:length(wake_vortex_panels)
    ]

    # axial components
    vx_ww = [
        A_ww[i, j][1] for i in 1:length(wake_vortex_panels),
        j in 1:length(wake_vortex_panels)
    ]

    # radial components
    vr_ww = [
        A_ww[i, j][2] for i in 1:length(wake_vortex_panels),
        j in 1:length(wake_vortex_panels)
    ]

    ##### ----- Induced Velcocities on Duct Wake ----- #####
    # - body to duct wake - #
    A_dwb = [
        assemble_induced_velocity_matrices(mesh_dwb[i, j], body_panels, duct_wake_panels)
        for i in 1:1, j in 1:1
    ]

    # axial components
    vx_dwb = [A_dwb[i, j][1] for i in 1:1, j in 1:1]

    # radial components
    vr_dwb = [A_dwb[i, j][2] for i in 1:1, j in 1:1]

    # - rotor to duct wake - #
    A_dwr = [
        assemble_induced_velocity_matrices(
            mesh_dwr[i, j], rotor_source_panels[j], duct_wake_panels
        ) for i in 1:1, j in 1:length(rotor_source_panels)
    ]

    # axial components
    vx_dwr = [A_dwr[i, j][1] for i in 1:1, j in 1:length(rotor_source_panels)]

    # radial components
    vr_dwr = [A_dwr[i, j][2] for i in 1:1, j in 1:length(rotor_source_panels)]

    # - wake to duct wake - #
    A_dww = [
        assemble_induced_velocity_matrices(
            mesh_dww[i, j], wake_vortex_panels[j], duct_wake_panels
        ) for i in 1:1, j in 1:length(wake_vortex_panels)
    ]

    # axial components
    vx_dww = [A_dww[i, j][1] for i in 1:1, j in 1:length(wake_vortex_panels)]

    # radial components
    vr_dww = [A_dww[i, j][2] for i in 1:1, j in 1:length(wake_vortex_panels)]

    ##### ----- Induced Velcocities on Hub Wake ----- #####
    # - body to hub wake - #
    A_hwb = [
        assemble_induced_velocity_matrices(mesh_hwb[i, j], body_panels, hub_wake_panels) for
        i in 1:1, j in 1:1
    ]

    # axial components
    vx_hwb = [A_hwb[i, j][1] for i in 1:1, j in 1:1]

    # radial components
    vr_hwb = [A_hwb[i, j][2] for i in 1:1, j in 1:1]

    # - rotor to hub wake - #
    A_hwr = [
        assemble_induced_velocity_matrices(
            mesh_hwr[i, j], rotor_source_panels[j], hub_wake_panels
        ) for i in 1:1, j in 1:length(rotor_source_panels)
    ]

    # axial components
    vx_hwr = [A_hwr[i, j][1] for i in 1:1, j in 1:length(rotor_source_panels)]

    # radial components
    vr_hwr = [A_hwr[i, j][2] for i in 1:1, j in 1:length(rotor_source_panels)]

    # - wake to hub wake - #
    A_hww = [
        assemble_induced_velocity_matrices(
            mesh_hww[i, j], wake_vortex_panels[j], hub_wake_panels
        ) for i in 1:1, j in 1:length(wake_vortex_panels)
    ]

    # axial components
    vx_hww = [A_hww[i, j][1] for i in 1:1, j in 1:length(wake_vortex_panels)]

    # radial components
    vr_hww = [A_hww[i, j][2] for i in 1:1, j in 1:length(wake_vortex_panels)]

    ## -- Miscellaneous Values for Indexing -- ##

    # - Get rotor panel edges and centers - #
    rotor_panel_edges = [rgrid[rotor_indices_in_wake[i], 1:length(rpe)] for i in 1:nrotor]
    rotor_panel_centers = [rotor_source_panels[i].panel_center[:, 2] for i in 1:nrotor]

    rotor_panel_edges = reduce(hcat, rotor_panel_edges)
    rotor_panel_centers = reduce(hcat, rotor_panel_centers)

    # get the total number of vortex panels on the bodies
    num_body_panels = length(b_bf)

    return (;
        converged=[false],
        #freestream
        freestream,
        #reference for post process
        reference_parameters,
        # body_geometry, # body geometry
        # - rotors
        blade_elements, # blade elements
        num_rotors=nrotor,
        rotor_panel_edges,
        rotor_panel_centers,
        # panels
        rotor_indices_in_wake,
        rotor_indices_on_duct,
        rotor_indices_on_hub,
        num_wake_x_panels=num_wake_x_panels,
        num_body_panels,
        ductTE_index=tip_gaps[1] == 0.0 ? ductTE_index : nothing,
        hubTE_index=!nohub ? hubTE_index : nothing,
        ductwakeidx=ductwakeidx,
        hubwakeidx=hubwakeidx,
        # body_panels, # body paneling
        # rotor_source_panels, # rotor paneling
        # wake_vortex_panels, # wake paneling
        # - unit induced velocities (INCLUDING PANEL LENGTH)
        A_bb_raw, #for initialization
        A_bb, # body to body
        b_bf, # freestream contribution to body boundary conditions
        A_br, # rotor to body (total)
        A_bw, # wake to body (total)
        A_rb, # body to rotor
        A_wb, # body to wake
        # A_wba, # body to "wake"
        A_wr, # rotor to wake
        A_ww, # wake to wake
        A_dwb, # body to duct wake
        A_dwr, # rotor to duct wake
        A_dww, # wake to duct wake
        A_hwb, # body to hub wake
        A_hwr, # rotor to hub wake
        A_hww, # wake to hub wake
        vx_rb, # body to rotor (x-direction)
        vr_rb, # body to rotor (r-direction)
        vx_rr, # rotor to rotor (x-direction)
        vr_rr, # rotor to rotor ( r-direction)
        vx_rw, # wake to rotor (x-direction)
        vr_rw, # wake to rotor ( r-direction)
        vx_wb, # body to wake (x-direction)
        vr_wb, # body to wake ( r-direction)
        # vx_wba, # body to "wake" (x-direction)
        # vr_wba, # body to "wake" ( r-direction)
        vx_wr, # rotor to wake (x-direction)
        vr_wr, # rotor to wake ( r-direction)
        vx_ww, # wake to wake (x-direction)
        vr_ww, # wake to wake ( r-direction)
        vx_dwb, # body to duct wake (x-direction)
        vr_dwb, # body to duct wake ( r-direction)
        vx_dwr, # rotor to duct wake (x-direction)
        vr_dwr, # rotor to duct wake ( r-direction)
        vx_dww, # wake to duct wake (x-direction)
        vr_dww, # wake to duct wake ( r-direction)
        vx_hwb, # body to hub wake (x-direction)
        vr_hwb, # body to hub wake ( r-direction)
        vx_hwr, # rotor to hub wake (x-direction)
        vr_hwr, # rotor to hub wake ( r-direction)
        vx_hww, # wake to hub wake (x-direction)
        vr_hww, # wake to hub wake ( r-direction)
        # kutta condition indices (body of rev=false for duct, true for hub)
        kutta_idxs=get_kutta_indices([false; true], mesh_bb),
        # operating conditions
        Vinf=freestream.Vinf, # freestream parameters
        # - Debugging/Plotting
        duct_coordinates=(noduct ? nothing : t_duct_coordinates),
        hub_coordinates=(nohub ? nothing : rp_hub_coordinates),
        isduct=!noduct,
        ishub=!nohub,
        body_panels,
        rotor_source_panels,
        wakexgrid=xgrid[:, 1:length(rpe)],
        wakergrid=rgrid[:, 1:length(rpe)],
        wake_vortex_panels,
        duct_wake_panels,
        hub_wake_panels,
        # wake_affect_panels,
        mesh_bb,
        mesh_rb,
        mesh_br,
        mesh_bw,
        mesh_rw,
    )
end

"""
    initialize_states(inputs)

Calculate an initial guess for the state variables
"""
function initialize_states(inputs)

    # - Initialize body vortex strengths (rotor-off linear problem) - #
    gamb = solve_body_system(inputs.A_bb_raw, inputs.b_bf, inputs.kutta_idxs) # get circulation strengths from solving body to body problem

    # - Initialize blade circulation and source strengths (assume open rotor) - #

    # get floating point type
    TF = promote_type(
        eltype(inputs.blade_elements[1].chords),
        eltype(inputs.blade_elements[1].twists),
        eltype(inputs.blade_elements[1].Omega),
        eltype(gamb),
    )

    # get problem dimensions (number of radial stations x number of rotors)
    nr = length(inputs.blade_elements[1].rbe)
    nrotor = length(inputs.blade_elements)
    nxwake = inputs.num_wake_x_panels

    # initialize outputs
    Gamr = zeros(TF, nr, nrotor)
    sigr = zeros(TF, nr, nrotor)
    gamw = zeros(TF, nr + 1, nxwake)
    # wake_vortex_strengths = repeat(gamw; inner=(1, inputs.num_wake_x_panels))
    vx_rotor, vr_rotor, vtheta_rotor = calculate_induced_velocities_on_rotors(
        inputs.blade_elements,
        Gamr,
        inputs.vx_rw,
        inputs.vr_rw,
        gamw,
        inputs.vx_rr,
        inputs.vr_rr,
        sigr,
        inputs.vx_rb,
        inputs.vr_rb,
        gamb,
    )
    # println("init vxrotor: ", vx_rotor)
    # println("init vrrotor: ", vr_rotor)

    Wx_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor = reframe_rotor_velocities(
        vx_rotor,
        vr_rotor,
        vtheta_rotor,
        inputs.Vinf,
        inputs.blade_elements[1].Omega,
        inputs.rotor_panel_centers,
    )

    # initialize circulation and source panel strengths
    calculate_gamma_sigma!(
        Gamr,
        sigr,
        inputs.blade_elements,
        Wm_rotor,
        Wtheta_rotor,
        Wmag_rotor,
        inputs.freestream,
    )

    # - update wake strengths - #
    calculate_wake_vortex_strengths!(
        Gamr, gamw, similar(sigr) .= 0.0, similar(gamb) .= 0.0, inputs; debug=false
    )

    # - Combine initial states into one vector - #
    states = vcat(
        gamb,               # body vortex panel strengths
        reduce(vcat, gamw), # wake vortex sheet strengths
        reduce(vcat, Gamr), # rotor circulation strengths
        reduce(vcat, sigr), # rotor source panel strengths
    )

    return states
end
