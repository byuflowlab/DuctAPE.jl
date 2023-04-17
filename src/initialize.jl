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
    freestream;
    finterp=fm.akima,
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

    # - Methods - #
    # method for the body to body problem
    # body_method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true])

    # method for vortex panels
    # vortex_method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true])

    # method for source panels
    # source_method = ff.AxisymmetricProblem(Source(Constant()), Neumann(), [true])

    #------------------------------------#
    # Discretize Wake and Repanel Bodies #
    #------------------------------------#
    if hub_coordinates == nothing
        nohub = true
        hub_coordinates = [
            minimum(duct_coordinates[:, 1]) 0.0
            maximum(duct_coordinates[:, 1]) 0.0
        ]
    else
        nohub = false
    end

    # - Discretize Wake x-coordinates - #
    # also returns indices of rotor locations in the wake
    xwake, rotor_indices = discretize_wake(
        duct_coordinates,
        hub_coordinates,
        rotor_parameters.xrotor,
        paneling_constants.wake_length,
        paneling_constants.npanels,
    )

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

    #don't include the dummy sheet used for tip gap if tipgap is present
    if tip_gaps[1] == 0.0
        xwakegrid = view(xgrid, :, :)
        rwakegrid = view(rgrid, :, :)
    else
        xwakegrid = view(xgrid, :, 1:(length(rwake) - 1))
        rwakegrid = view(rgrid, :, 1:(length(rwake) - 1))
    end

    # generate wake sheet paneling
    wake_vortex_panels = generate_wake_panels(xwakegrid, rwakegrid)

    #------------------------------------------#
    # Generate Rotor Panels and Blade Elements #
    #------------------------------------------#

    # rotor source panel objects
    rotor_source_panels = [
        generate_rotor_panels(rotor_parameters[i].xrotor, rwakegrid[rotor_indices[i], :])
        for i in 1:nrotor
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
        ) for i in 1:1, j in 1:length(rotor_source_panels)
    ]

    # - wake to body - #
    A_bw = [
        assemble_induced_velocity_on_body_matrix(
            mesh_bw[i, j], [wake_vortex_panels[j]], body_panels; singularity="vortex"
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

    ## -- Miscellaneous Values for Indexing -- ##

    # - Get rotor panel edges and centers - #
    rotor_panel_edges = [rwakegrid[rotor_indices[i], :] for i in 1:nrotor]
    rotor_panel_centers = [rotor_source_panels[i].panel_center[:, 2] for i in 1:nrotor]

    rotor_panel_edges = reduce(hcat, rotor_panel_edges)
    rotor_panel_centers = reduce(hcat, rotor_panel_centers)

    # get the total number of vortex panels on the bodies
    num_body_panels = length(b_bf)

    return (;
        converged=[false],
        # body_geometry, # body geometry
        # - rotors
        blade_elements, # blade elements
        num_rotors=nrotor,
        rotor_panel_edges,
        rotor_panel_centers,
        # panels
        rotor_indices,
        num_wake_x_panels=length(xwake) - 1,
        num_body_panels,
        # body_panels, # body paneling
        # rotor_source_panels, # rotor paneling
        # wake_vortex_panels, # wake paneling
        # - unit induced velocities (INCLUDING PANEL LENGTH)
        A_bb, # body to body
        b_bf, # freestream contribution to body boundary conditions
        A_br, # rotor to body (total)
        A_bw, # wake to body (total)
        vx_rb, # body to rotor (x-direction)
        vr_rb, # body to rotor (r-direction)
        vx_rr, # rotor to rotor (x-direction)
        vr_rr, # rotor to rotor ( r-direction)
        vx_rw, # wake to rotor (x-direction)
        vr_rw, # wake to rotor ( r-direction)
        # kutta condition indices
        kutta_idxs=get_kutta_indices([false; true], mesh_bb),
        # operating conditions
        Vinf=freestream.Vinf, # freestream parameters
        # - Debugging/Plotting
        t_duct_coordinates,
        rp_hub_coordinates,
        body_panels,
        rotor_source_panels,
        wake_vortex_panels,
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
    A = inputs.A_bb # AIC matrix for body to body problem
    b = inputs.b_bf # right hand side for body to body problem
    gamb = solve_body_system(A, b, inputs.kutta_idxs) # get circulation strengths from solving body to body problem

    # - Initialize blade circulation and source strengths (assume open rotor) - #
    # use rotor rotation for tangential velocity of each blade section
    Wθ_init = -inputs.rotor_panel_centers .* inputs.blade_elements.Omega'
    # use freestream magnitude as meridional velocity at each blade section
    Wm_init = similar(Wθ_init) .= inputs.Vinf
    # magnitude is simply freestream and rotation
    W_init = sqrt.(Wθ_init .^ 2 .+ Wm_init .^ 2)
    # initialize circulation and source panel strengths
    Gamr, sigr = calculate_gamma_sigma(inputs.blade_elements, Wm_init, Wθ_init, W_init)

    # - Initialize wake vortex strengths - #
    gamw = initialize_wake_vortex_strengths(
        inputs.Vinf,
        Gamr,
        inputs.blade_elements.Omega,
        inputs.blade_elements.B,
        inputs.rotor_panel_edges,
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
