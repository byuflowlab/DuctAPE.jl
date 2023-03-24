
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
- `nwake=length(rotor_parameters[1].rblade)`: Number of radial stations to use when defining
    the wake
- `finterp=FLOWMath.akima`: Method used to interpolate the duct and hub geometries.
- `xwake`: May be used to define a custom set of x-coordinates for the wake.
- `nhub`: Number of panels between the hub leading edge and the first rotor.
- `nduct_inner`: Number of panels between the duct leading edge and the first rotor.
- `nduct_outer`: Number of panels on the duct outer surface.

NOTE: The `nwake`, `xwake`, `nhub`, `nduct_inner`, and `nduct_outer` arguments should
always be provided when performing gradient-based optimizatoin in order to ensure that the
resulting paneling is consistent between optimization iterations.
"""
function initialize_parameters(duct_coordinates, hub_coordinates, rotor_parameters, freestream;
    wake_length=1.0, nwake = length(rotor_parameters[1].rblade), finterp = FLOWMath.akima,
    xwake = default_discretization(duct_coordinates, hub_coordinates, rotor_parameters,
        wake_length, nwake),
    nhub = nothing, nduct_inner = nothing, nduct_outer = nothing)

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
        duct_coordinates, hub_coordinates, xwake, nhub, nduct_inner, nduct_outer;
        finterp=finterp)

    # generate body paneling
    body_panels = generate_body_panels(updated_duct_coordinates, updated_hub_coordinates)

    # --- First Rotor --- #

    # initialize blade elements for first rotor
    first_blade_elements = generate_blade_elements(
        rotor_parameters[1].B,
        rotor_parameters[1].omega,
        rotor_parameters[1].xrotor,
        rotor_parameters[1].rblade,
        rotor_parameters[1].chords,
        rotor_parameters[1].twists,
        rotor_parameters[1].solidity,
        rotor_parameters[1].airfoils,
        updated_duct_coordinates,
        updated_hub_coordinates,
        rwake)

    # initialize panels for first rotor
    first_rotor_panels = generate_rotor_panels(first_blade_elements, source_method)
    first_dummy_rotor_panels = generate_dummy_rotor_panels(first_blade_elements, source_method)

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
            rotor_parameters[i].omega,
            rotor_parameters[i].xrotor,
            rotor_parameters[i].rblade,
            rotor_parameters[i].chords,
            rotor_parameters[i].twists,
            rotor_parameters[i].solidity,
            rotor_parameters[i].airfoils,
            updated_duct_coordinates,
            updated_hub_coordinates,
            view(rgrid, rotor_indices[i], :)
            )

        # update aft rotor panels
        rotor_panels[i] = generate_rotor_panels(blade_elements[i], source_method)
        dummy_rotor_panels[i] = generate_dummy_rotor_panels(blade_elements[i], source_method)

    end

    # --- Meshes --- #

    # body to body
    body_mesh = ff.generate_mesh(body_method, body_panels)

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

"""
    discretize_wake(duct_coordinates, hub_coordinates, rotor_parameters, wake_length, nwake)

Calculate wake x-coordinates.
"""
function discretize_wake(duct_coordinates, hub_coordinates, rotor_parameters, wake_length, nwake)

    # extract rotor locations
    xrotors = getproperty.(rotor_parameters, :xrotor)

    # extract duct leading and trailing edge locations
    xduct_le = minimum(duct_coordinates[:, 1])
    xduct_te = maximum(duct_coordinates[:, 1])

    # extract hub leading and trailing edge location
    xhub_le = minimum(hub_coordinates[:, 1])
    xhub_te = maximum(hub_coordinates[:, 1])

    # calculate duct chord
    duct_chord = max(xduct_te, xhub_te) - min(xduct_le, xhub_le)

    # dimensionalize wake_length
    wake_length = duct_chord * wake_length

    # ensure rotors are ordered properly
    @assert issorted(xrotors)

    # make sure that all rotors are inside the duct
    for xrotor in xrotors
        @assert xrotor > xduct_le "Rotor is in front of duct leading edge."
        @assert xrotor < xduct_te "Rotor is behind duct trailing edge."
        @assert xrotor > xhub_le "Rotor is in front of hub leading edge."
        @assert xrotor < xhub_te "Rotor is behind hub trailing edge."
    end

    # combine all discrete locations into one ordered array
    if isapprox(xhub_te, xduct_te)
        xd = vcat(xrotors, xhub_te, xhub_te + wake_length)
    elseif xduct_te < xhub_te
        xd = vcat(xrotors, xduct_te, xhub_te, xhub_te + wake_length)
    else
        xd = vcat(xrotors, xhub_te, xduct_te, xduct_te + wake_length)
    end

    # find (approximate) hub and tip locations
    _, leidx = findmin(view(duct_coordinates, :, 1))
    _, ihub = findmin(x -> abs(x - xrotors[1]), view(hub_coordinates, :, 1))
    _, iduct = findmin(x -> abs(x - xrotors[1]), view(duct_coordinates, 1:leidx, 1))
    Rhub = hub_coordinates[ihub, 2]
    Rtip = duct_coordinates[iduct, 2]

    # calculate number of panels per unit length
    panel_density = nwake / (Rtip - Rhub)

    # calculate number of panels for each discrete section
    npanels = [ceil(Int, (xd[i+1] - xd[i]) * panel_density) for i = 1:length(xd)-1]

    # calculate indices for the start of each discrete section
    indices = cumsum(vcat(1, npanels))

    # construct x-discretization
    xwake = similar(xd, sum(npanels)+1)
    for i = 1:length(xd)-1
        xrange = range(xd[i], xd[i+1], length=npanels[i]+1)
        xwake[indices[i] : indices[i+1]] .= xrange
    end

    # return dimensionalized wake x-coordinates
    return xwake
end