"""
    initialize_parameters(duct_coordinates, hub_coordinates, rotor_parameters, freestream; kwargs...)

Initializes the geometry, panels, and aerodynamic influence coefficient matrices.

# Arguments:
- `duct_coordinates::Array{Float64,2}` : duct coordinates, starting from trailing edge, going clockwise
- `hub_coordinates::Array{Float64,2}` : hub coordinates, starting from leading edge
- `rotor_parameters`: named tuple of rotor parameters
- `freestream`: freestream parameters

# Keyword Arguments
- `wake_length=1.0` : non-dimensional length (based on maximum duct chord) that the wake
    extends past the furthest trailing edge.
- `nwake_sheets=length(rotor_parameters[1].rblade)`: Number of radial stations to use when defining
    the wake
- `finterp=FLOWMath.akima`: Method used to interpolate the duct and hub geometries.
- `xwake`: May be used to define a custom set of x-coordinates for the wake.
- `nhub`: Number of panels between the hub leading edge and the first rotor.
- `nduct_inner`: Number of panels between the duct leading edge and the first rotor.
- `nduct_outer`: Number of panels on the duct outer surface.

NOTE: The `nwake_sheets`, `xwake`, `nhub`, `nduct_inner`, and `nduct_outer` arguments should
always be provided when performing gradient-based optimizatoin in order to ensure that the
resulting paneling is consistent between optimization iterations.
"""
function initialize_parameters(
    duct_coordinates,
    hub_coordinates,
    rotor_parameters,
    freestream;
    wake_length=1.0,
    nwake_sheets=length(rotor_parameters[1].rblade),
    finterp=FLOWMath.akima,
    xwake=default_discretization(
        duct_coordinates, hub_coordinates, rotor_parameters, wake_length, nwake_sheets
    ),
    nhub=nothing,
    nduct_inner=nothing,
    nduct_outer=nothing,
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

    # update the body paneling to match the wake discretization
    updated_duct_coordinates, updated_hub_coordinates = update_body_geometry(
        duct_coordinates,
        hub_coordinates,
        xwake,
        nhub,
        nduct_inner,
        nduct_outer;
        finterp=finterp,
    )

    # generate body paneling
    body_panels = generate_body_panels(updated_duct_coordinates, updated_hub_coordinates)

    # --- First Rotor --- #

    # initialize blade elements for first rotor
    first_blade_elements = generate_blade_elements(
        rotor_parameters[1].B,
        rotor_parameters[1].Omega,
        rotor_parameters[1].xrotor,
        rotor_parameters[1].rblade,
        rotor_parameters[1].chords,
        rotor_parameters[1].twists,
        rotor_parameters[1].solidity,
        rotor_parameters[1].airfoils,
        updated_duct_coordinates,
        updated_hub_coordinates,
        rwake,
    )

    # initialize panels for first rotor
    first_rotor_panels = generate_rotor_panels(first_blade_elements, source_method)
    first_dummy_rotor_panels = generate_dummy_rotor_panels(
        first_blade_elements, source_method
    )

    # --- Wake Grid --- #

    # Generate a streamline-aligned wake grid
    xgrid, rgrid, rotor_indices = generate_wake_grid(xwake, rwake, xrotors, wake_length)

    # construct wake panels
    wake_panels = generate_wake_panels(xgrid, rgrid; method=vortex_method)

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
            rotor_parameters[i].Omega,
            rotor_parameters[i].xrotor,
            rotor_parameters[i].rblade,
            rotor_parameters[i].chords,
            rotor_parameters[i].twists,
            rotor_parameters[i].solidity,
            rotor_parameters[i].airfoils,
            updated_duct_coordinates,
            updated_hub_coordinates,
            view(rgrid, rotor_indices[i], :),
        )

        # update aft rotor panels
        rotor_panels[i] = generate_rotor_panels(blade_elements[i], source_method)
        dummy_rotor_panels[i] = generate_dummy_rotor_panels(
            blade_elements[i], source_method
        )
    end

    # --- Meshes --- #

    # body to body
    body_mesh = ff.generate_mesh(body_method, body_panels)

    # body to rotor
    mesh_rtb = generate_one_way_mesh.(Ref(body_panels), dummy_rotor_panels)

    # body to wake
    mesh_btw = generate_one_way_mesh.(Ref(body_panels), wake_panels)

    # rotor to body
    mesh_btr = generate_one_way_mesh.(rotor_panels, Ref(body_panels); singularity="source")

    # rotor to rotor
    mesh_rr =
        generate_one_way_mesh.(rotor_panels, dummy_rotor_panels'; singularity="source")

    # rotor to wake
    mesh_rtw = generate_one_way_mesh.(rotor_panels, wake_panels'; singularity="source")

    # wake to body
    mesh_wtb = generate_one_way_mesh.(wake_panels, Ref(body_panels))

    # wake to rotor
    mesh_wtr = generate_one_way_mesh.(wake_panels, dummy_rotor_panels')

    # wake to wake
    mesh_ww = generate_one_way_mesh.(wake_panels, wake_panels')

    # --- Coefficient Matrices --- #

    body_system = ff.generate_inviscid_system(body_method, body_panels, body_mesh)

    # body to body
    A_bb = body_system.A

    # freestream to body
    b_bf = body_system.b .* freestream.Vinf

    # body to rotor
    A_btr =
        assemble_one_way_coefficient_matrix.(mesh_rtb, Ref(body_panels), dummy_rotor_panels)
    vx_rb = getindex.(A_btr, 1)
    vr_rb = getindex.(A_btr, 2)

    # rotor to body
    A_rtb =
        assemble_one_way_coefficient_matrix.(
            mesh_btr, rotor_panels, Ref(body_panels); singularity="source"
        )
    vx_br = getindex.(A_rtb, 1)
    vr_br = getindex.(A_rtb, 2)

    # rotor to rotor
    A_rr =
        assemble_one_way_coefficient_matrix.(
            mesh_rr, rotor_panels, dummy_rotor_panels'; singularity="source"
        )
    vx_rr = getindex.(A_rr, 1)
    vr_rr = getindex.(A_rr, 2)

    # rotor to wake
    A_rtw =
        assemble_one_way_coefficient_matrix.(
            mesh_rtw, rotor_panels, wake_panels'; singularity="source"
        )
    vx_wr = getindex.(A_rtw, 1)
    vr_wr = getindex.(A_rtw, 2)

    # wake to body
    A_wtb = assemble_one_way_coefficient_matrix.(mesh_wtb, wake_panels, Ref(body_panels))
    vx_bw = getindex.(A_wtb, 1)
    vr_bw = getindex.(A_wtb, 2)

    # wake to rotor
    A_wtr = assemble_one_way_coefficient_matrix.(mesh_wtr, wake_panels, dummy_rotor_panels')
    vx_rw = getindex.(A_wtr, 1)
    vr_rw = getindex.(A_wtb, 2)

    return (
        # body
        body_geometry=body_geometry, # body geometry
        # rotors
        blade_elements=blade_elements, # blade elements
        rotor_indices=rotor_indices, # rotor locations
        # panels
        body_panels=body_panels, # body paneling
        rotor_panels=rotor_panels, # rotor paneling
        wake_panels=wake_panels, # wake paneling
        # aerodynamic influence coefficients
        A_bb=A_bb, # body to body
        b_bf=b_bf, # freestream contribution to body boundary conditions
        vx_rb=vx_rb, # body to rotor (x-direction)
        vr_rb=vr_rb, # body to rotor (r-direction)
        vx_wb=vx_wb, # body to wake (x-direction)
        vr_wb=vr_wb, # body to wake ( r-direction)
        vx_br=vx_br, # rotor to body (x-direction)
        vr_br=vr_br, # rotor to body ( r-direction)
        vx_rr=vx_rr, # rotor to rotor (x-direction)
        vr_rr=vr_rr, # rotor to rotor ( r-direction)
        vx_wr=vx_wr, # rotor to wake (x-direction)
        vr_wr=vr_wr, # rotor to wake ( r-direction)
        vx_bw=vx_bw, # wake to body (x-direction)
        vr_bw=vr_bw, # wake to body ( r-direction)
        vx_rw=vx_rw, # wake to rotor (x-direction)
        vr_rw=vr_rw, # wake to rotor ( r-direction)
        vx_ww=vx_ww, # wake to wake (x-direction)
        vr_ww=vr_ww, # wake to wake ( r-direction)
        # operating conditions
        freestream=freestream, # freestream parameters
    )
end
