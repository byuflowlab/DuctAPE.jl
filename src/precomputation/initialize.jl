"""
    generate_geometry(
    duct_coordinates,
    hub_coordinates,
    paneling_constants,
    rotorstator_parameters;
    kwargs...
    )

Generates repaneled geometry based on the input coordinates and parameters.

Note that the contents of `paneling_constants` will for the most part be required to remain constant across optimization iterations if gradient-based optimization is used.

# Arguments:
- `duct_coordinates::Array{Float64,2}` : duct coordinates, starting from trailing edge, going clockwise
- `hub_coordinates::Array{Float64,2}` : hub coordinates, starting from leading edge
- `paneling_constants.NTuple' : Named tuple containing parameters necessary for repaneling including:
    - `:npanels::Vector{Int}` : Vector of the number of panels to place in the axial direction aft of the first rotor and between discrete system locations such as rotor locations, solid body trailing edges, and the end of the wake.  For example, if there is a single rotor and the duct and center body trailing edges align, 'npanels' would be a 2-element vector containing 1) the number of panels from the rotor to the duct trailing edge and 2) the number of panels from the duct trailing edge to the end of the wake.
    - `:nhub_inlet::Int` : number of nodes to include from the center body leading edge to the foremost rotor location
    - `:nduct_inlet::Int` : number of nodes to include from the duct leading edge to the foremost rotor location
    - `:wake_length::Float` : non-dimensional (relative to duct chord length) distance to extend wake aft of duct trailing edge
    - `:nwake_sheets::Int` : number of wake sheets (1 more than the number of blade elements to use)
- `rotorstator_parameters::Vector{NTuple}`: Vector of named tuple of rotor parameters
    - 'rotorzloc::Float` : axial position of rotor
    - 'r::Vector{Float}` : non-dimensional radial positions of stations along the blade
    - 'chords::Vector{Float}` : dimensional chord distribution
    - 'twists::Vector{Float}` : dimensional twist distribution (radians)
    - 'airfoils::Vector{::AFType}` : Vector of airfoil objects.
    - 'Rtip::Float` : dimensional tip radius
    - 'Rhub::Float` : dimensional hub radius
    - 'tip_gap::Float` : gap from rotor tip to duct wall (not yet functional)
    - 'B::Int` : number of blades
    - 'Omega::Float` : rotor rotation rate in radians per second
    - 'fliplift::Bool` : boolean as to whether to flip the looked up cl for the blade element (may be useful for stator sections).


# Keyword Arguments
- `finterp=FLOWMath.akima`: Method used to interpolate the duct and hub geometries.
- `autoshiftduct::Bool=false' : Boolean for whether to automatically shift the duct to the correct radial position based on foremost rotor tip radius and tip gap.

# Returns:
- `system_geometry::NTuple` : Named tuple of geometry items required as an input to the `precomputed_inputs()` function:
    - `duct_coordinates::Array{Float64,2}` : duct coordinates, starting from trailing edge, going clockwise
    - `hub_coordinates::Array{Float64,2}` : hub coordinates, starting from leading edge
    - `Rtips::Vector{Float}` : values of rotor tip radii
    - `Rhubs::Vector{Float}` : values of rotor hub radii
    - `rpe::Vector{Float}` : rotor node locations
    - `grid::Array{Float,3}` : nodes locations for wake, [z-r dimension, axial direction, radial direction]
    - `zwake::Vector{Float}` : axial positions of repaneled geometry
    - `rwake::Vector{Float}` : radial positions of first wake nodes
    - `rotor_indices_in_wake::Vector{Int}` : indices in axial direction of where rotors lie in the wake
    - `ductTE_index::Int` : index in axial direction of where duct trailing edge hits in the wake
    - `hubTE_index::Int` : index in axial direction of where center body trailing edge hits in the wake
    - `nohub::Bool` : flag if hub coordinates = `nothing`
    - `noduct::Bool` : flag if duct coordinates = `nothing`
    - `rotoronly::Bool` : flag if no bodies are defined (not currently used)


"""
function generate_geometry(
    duct_coordinates,
    hub_coordinates,
    paneling_constants,
    rotorstator_parameters; #vector of named tuples
    finterp=fm.akima,
    autoshiftduct=false,
)

    ## -- Rename for Convenience -- ##
    # promoted type
    TF = promote_type(
        eltype(rotorstator_parameters[1].chords),
        eltype(rotorstator_parameters[1].twists),
        eltype(rotorstator_parameters[1].Omega),
        eltype(duct_coordinates),
        eltype(hub_coordinates),
    )

    # number of rotors
    num_rotors = length(rotorstator_parameters)

    # - Determine if duct or center body is not present - #
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

    #------------------------------------#
    # Discretize Wake and Repanel Bodies #
    #------------------------------------#

    # - Discretize Wake z-coordinates - #
    # also returns indices of rotor locations and duct and center body trailng edges in the wake
    zwake, rotor_indices_in_wake, ductTE_index, hubTE_index = discretize_wake(
        duct_coordinates,
        hub_coordinates,
        rotorstator_parameters.rotorzloc,
        paneling_constants.wake_length,
        paneling_constants.npanels,
    )

    # - Repanel Bodies - #
    rp_duct_coordinates, rp_hub_coordinates = repanel_bodies(
        duct_coordinates,
        hub_coordinates,
        zwake,
        paneling_constants.nhub_inlet,
        paneling_constants.nduct_inlet;
        finterp=finterp,
    )

    #-----------------------------------#
    # Position Duct and Get Rotor Radii #
    #-----------------------------------#
    # check that tip gap isn't too small
    # set vector of tip gaps to zero for initialization (any overrides to zero thus taken c of automatically)
    tip_gaps = zeros(eltype(rotorstator_parameters.tip_gap), num_rotors)
    if rotorstator_parameters[1].tip_gap != 0.0
        if rotorstator_parameters[1].tip_gap < 1e-4
            @warn "You have selected a tip gap for the foremost rotor that is smaller than 1e-4. Overriding to 0.0 to avoid near singularity issues."
        else
            tip_gaps[1] = rotorstator_parameters[1].tip_gap
        end
    end

    # can't have non-zero tip gaps for aft rotors
    for ir in 2:num_rotors
        if rotorstator_parameters[ir].tip_gap != 0.0
            @warn "DuctAPE does not currently have capabilities for adding tip gap to any but the foremost rotor. Overriding to 0.0."
        else
            tip_gaps[ir] = rotorstator_parameters[ir].tip_gap
        end
    end

    # if hub was nothing, set hub radius to dimensional inner rotor radius
    if nohub
        rp_hub_coordinates[:, 2] .=
            rotorstator_parameters[1].r[1] * rotorstator_parameters[1].Rtip
    end

    if autoshiftduct
        place_duct!(
            rp_duct_coordinates,
            rotorstator_parameters[1].Rtip,
            rotorstator_parameters[1].rotorzloc,
            tip_gaps[1],
        )
    end

    # - Get dimensional rotor blade hub and tip radii based on duct and center body geometry and rotor z locations - #
    Rtips, Rhubs = get_blade_ends_from_body_geometry(
        rp_duct_coordinates, rp_hub_coordinates, tip_gaps, rotorstator_parameters.rotorzloc
    )

    for i in 1:num_rotors
        @assert Rtips[i] > Rhubs[i] "Rotor #$i Tip Radius is set to be less than its Hub Radius."
    end

    #----------------------------------#
    # Generate Discretized Wake Sheets #
    #----------------------------------#
    # get discretization of wakes at leading rotor position

    for i in 1:num_rotors
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
    grid = initialize_wake_grid(rp_duct_coordinates, rp_hub_coordinates, zwake, rwake)

    # Relax "Grid"
    relax_grid!(grid; max_iterations=100, tol=1e-9, verbose=false)

    return (;
        duct_coordinates=rp_duct_coordinates,
        hub_coordinates=rp_hub_coordinates,
        Rtips,
        Rhubs,
        rpe,
        grid,
        zwake,
        rwake,
        rotor_indices_in_wake,
        ductTE_index,
        hubTE_index,
        nohub,
        noduct,
        rotoronly,
    )
end

"""
    precomputed_inputs(duct_coordinates, hub_coordinates, paneling_constants, rotorstator_parameters, freestream, reference_parameters; kwargs...)
    precomputed_inputs(system_geometry, rotorstator_parameters, freestream, reference_parameters)

Initializes the geometry, panels, and aerodynamic influence coefficient matrices.

# Arguments:
- `duct_coordinates::Array{Float64,2}` : duct coordinates, starting from trailing edge, going clockwise
- `hub_coordinates::Array{Float64,2}` : hub coordinates, starting from leading edge
- `paneling_constants.wake_length=1.0` : non-dimensional length (based on maximum duct chord) that the wake extends past the furthest trailing edge.
- `geometry::NTuple' : Named tuple containing the outputs of the generate_geometry() function
- `rotorstator_parameters`: named tuple of rotor parameters
- `freestream`: freestream parameters
- `reference_parameters`: reference parameters

# Keyword Arguments
- `finterp=FLOWMath.akima`: Method used to interpolate the duct and hub geometries.
- `autoshiftduct::Bool=false' : Boolean for whether to automatically shift the duct to the correct radial position based on foremost rotor tip radius and tip gap.

"""
function precomputed_inputs(
    duct_coordinates,
    hub_coordinates,
    paneling_constants,
    rotorstator_parameters, #vector of named tuples
    freestream,
    reference_parameters;
    finterp=fm.akima,
    autoshiftduct=false,
    debug=false,
)
    system_geometry = generate_geometry(
        duct_coordinates,
        hub_coordinates,
        paneling_constants,
        rotorstator_parameters; #vector of named tuples
        finterp=finterp,
        autoshiftduct=autoshiftduct,
    )

    return precomputed_inputs(
        system_geometry,
        rotorstator_parameters, #vector of named tuples
        freestream,
        reference_parameters;
        debug=debug,
    )
end

function precomputed_inputs(
    system_geometry,
    rotorstator_parameters, #vector of named tuples
    freestream,
    reference_parameters;
    debug=false,
)
    (;
        duct_coordinates,
        hub_coordinates,
        Rtips,
        Rhubs,
        rpe,
        grid,
        zwake,
        rwake,
        rotor_indices_in_wake,
        ductTE_index,
        hubTE_index,
        nohub,
        noduct,
        rotoronly,
    ) = system_geometry

    ## -- Rename for Convenience -- ##
    # promoted type
    TF = promote_type(
        eltype(rotorstator_parameters[1].chords),
        eltype(rotorstator_parameters[1].twists),
        eltype(rotorstator_parameters[1].Omega),
        eltype(duct_coordinates),
        eltype(hub_coordinates),
        eltype(freestream.Vinf),
    )

    # number of rotors
    num_rotors = length(rotorstator_parameters)

    #---------------------------------#
    #         GENERATE PANELS         #
    #---------------------------------#
    # generate body paneling
    if  nohub
        body_vortex_panels = generate_panels(duct_coordinates)
    else
        body_vortex_panels = generate_panels([duct_coordinates, hub_coordinates])
    end

    # generate wake sheet paneling
    wake_vortex_panels = generate_wake_panels(
        grid[1, :, 1:length(rpe)], grid[2, :, 1:length(rpe)]
    )

    # generate body wake panels for convenience (used in getting surface pressure of body wakes in post processing)
    duct_wake_panels = generate_panels(
        [grid[1, ductTE_index:end, end]'; grid[2, ductTE_index:end, end]']
    )
    hub_wake_panels = generate_panels(
        [grid[1, hubTE_index:end, 1]'; grid[2, hubTE_index:end, 1]']
    )

    # rotor source panel objects
    rotor_source_panels = [
        generate_rotor_panels(
            rotorstator_parameters[i].rotorzloc,
            grid[2, rotor_indices_in_wake[i], 1:length(rpe)],
        ) for i in 1:num_rotors
    ]

    #---------------------------------#
    #           BOOKKEEPING           #
    #---------------------------------#
    # book keep wake wall interfaces
    # !NOTE: assumes duct geometry given first in paneling.  could generalize this, but not worth it for now.
    if nohub
        rotor_indices_on_hub = nothing
        hubidsaftofrotors = nothing
    else
        rotor_indices_on_hub = [
            findlast(
                x -> x < rotorstator_parameters.rotorzloc[i],
                body_vortex_panels.controlpoint[1, :],
            ) for i in 1:length(rotorstator_parameters.rotorzloc)
        ]
        hwidraw = sort([body_vortex_panels.npanel; rotor_indices_on_hub])
        hubidsaftofrotors = reduce(
            vcat, [[(hwidraw[i] + 1):hwidraw[i + 1]] for i in 1:num_rotors]
        )
    end

    if noduct
        rotor_indices_on_duct = nothing
        ductidsaftofrotors = nothing
    else
        rotor_indices_on_duct = [
            findfirst(
                x -> x < rotorstator_parameters.rotorzloc[i],
                body_vortex_panels.controlpoint[1, :],
            ) - 1 for i in 1:length(rotorstator_parameters.rotorzloc)
        ]
        dwidraw = sort([0; rotor_indices_on_duct])
        ductidsaftofrotors = reverse(
            reduce(vcat, [[(dwidraw[i] + 1):dwidraw[i + 1]] for i in 1:num_rotors])
        )
    end

    # calculate radius dependent "constant" for wake strength calcualtion
    wakeK = get_wake_k(wake_vortex_panels)

    # Go through the wake panels and determine the index of the aftmost rotor infront and the blade node from which the wake strength is defined.
    rotorwakepanelid = ones(Int, wake_vortex_panels.totpanel, 2)
    num_wake_z_panels = length(zwake) - 1
    for i in 1:(length(rwake))
        rotorwakepanelid[(1 + (i - 1) * num_wake_z_panels):(i * num_wake_z_panels), 1] .= i
    end
    for (i, wn) in enumerate(eachcol(wake_vortex_panels.controlpoint))
        rotorwakepanelid[i, 2] = findlast(x -> x <= wn[1], rotorstator_parameters.rotorzloc)
    end

    # Go through the wake panels and determine the index of the aftmost rotor infront and the blade node from which the wake strength is defined.
    rotorwakeid = ones(Int, wake_vortex_panels.totnode, 2)
    num_wake_z_nodes = length(zwake)
    for i in 1:(rotorstator_parameters[1].nwake_sheets)
        rotorwakeid[(1 + (i - 1) * num_wake_z_nodes):(i * num_wake_z_nodes), 1] .= i
    end
    for (i, wn) in enumerate(eachcol(wake_vortex_panels.node))
        # TODO: DFDC geometry doesn't line up wake and rotor perfectly, so need a more robust option.
        # rotorwakeid[i, 2] = findlast(x -> x <= wn[1], rotorstator_parameters.rotorzloc)
        # TODO: current tests are passing, but look here if things break in the future.
        rotorwakeid[i, 2] = findmin(x -> abs(x - wn[1]), rotorstator_parameters.rotorzloc)[2]
    end

    # Go through the wake panels and determine the indices that have interfaces with the hub and wake
    hubwakeinterfaceid = 1:(hubTE_index - 1) #first rotor-wake-body interface is at index 1, this is also on the first wake sheet, so the hub trailing edge index in the zwake vector should be (or one away from, need to check) the last interface point
    ductwakeinterfaceid =
        num_wake_z_panels * (length(rwake) - 1) .+ (1:(ductTE_index - 1))

    #------------------------------------------#
    #          Generate Blade Elements         #
    #------------------------------------------#

    # rotor blade element objects
    blade_elements = [
        generate_blade_elements(
            rotorstator_parameters[i].B,
            rotorstator_parameters[i].Omega,
            rotorstator_parameters[i].rotorzloc,
            rotorstator_parameters[i].r,
            rotorstator_parameters[i].chords,
            rotorstator_parameters[i].twists,
            rotorstator_parameters[i].airfoils,
            Rtips[i],
            Rhubs[i],
            rotor_source_panels[i].controlpoint[2, :];
            fliplift=rotorstator_parameters[i].fliplift,
        ) for i in 1:num_rotors
    ]

    #---------------------------------#
    # Calculate Coefficient Matrices  #
    #---------------------------------#

    ##### ----- Induced Velcocities on Bodies ----- #####

    # - body to body - #
    # A_bb = doublet_panel_influence_matrix(body_vortex_panels.nodes, body_vortex_panels)
    # LHS = A_bb

    # -- Assemble LHS Matrix -- #
    # - Boundary on boundary influence coefficients - #
    AICn, AICt = vortex_aic_boundary_on_boundary(
        body_vortex_panels.controlpoint,
        body_vortex_panels.normal,
        body_vortex_panels.tangent,
        body_vortex_panels.node,
        body_vortex_panels.nodemap,
        body_vortex_panels.influence_length,
    )

    # - Boundary on internal psuedo control point influence coefficients - #
    AICpcp, unused = vortex_aic_boundary_on_field(
        body_vortex_panels.itcontrolpoint,
        body_vortex_panels.itnormal,
        body_vortex_panels.ittangent,
        body_vortex_panels.node,
        body_vortex_panels.nodemap,
        body_vortex_panels.influence_length,
    )

    # - Add Trailing Edge Gap Panel Influences - #
    add_te_gap_aic!(
        AICn,
        AICt,
        body_vortex_panels.controlpoint,
        body_vortex_panels.normal,
        body_vortex_panels.tangent,
        body_vortex_panels.tenode,
        body_vortex_panels.teinfluence_length,
        body_vortex_panels.tendotn,
        body_vortex_panels.tencrossn,
        body_vortex_panels.teadjnodeidxs,
    )

    add_te_gap_aic!(
        AICpcp,
        unused,
        body_vortex_panels.itcontrolpoint,
        body_vortex_panels.itnormal,
        body_vortex_panels.ittangent,
        body_vortex_panels.tenode,
        body_vortex_panels.teinfluence_length,
        body_vortex_panels.tendotn,
        body_vortex_panels.tencrossn,
        body_vortex_panels.teadjnodeidxs,
    )

    prelhs = assemble_lhs_matrix(AICn, AICpcp, body_vortex_panels; dummyval=1.0)

    # - LU Decomp - #
    LHS = lu!(prelhs, NoPivot(); check=false)
    # check if success
    lu_decomp_flag = issuccess(LHS)

    # - Freetream influence for RHS vector - #
    # Define freestream on panels
    Vs = freestream.Vinf * [1.0; 0.0] # axisymmetric, so no radial component

    vdnb = freestream_influence_vector(
        body_vortex_panels.normal, repeat(Vs; outer=(1, body_vortex_panels.totpanel))
    )
    vdnpcp = freestream_influence_vector(
        body_vortex_panels.itnormal,
        repeat(Vs; outer=(1, size(body_vortex_panels.itcontrolpoint, 2))),
    )

    # initial RHS includes only freestream values
    RHS = assemble_rhs_matrix(vdnb, vdnpcp, body_vortex_panels)

    # initial body strengths are the strengths without the rotor inductions
    gamb = LHS \ RHS

    # - rotor to body - #
    # preallocate AICnr and AICtr
    Cb = size(body_vortex_panels.controlpoint, 2)
    Nr = size(rotor_source_panels[1].node, 2)

    AICnr = zeros(TF, num_rotors, Cb, Nr)
    AICtr = zeros(TF, num_rotors, Cb, Nr)

    # loop through rotor objects
    for (j, rp) in enumerate(rotor_source_panels)
        source_aic!(
            view(AICnr, j, :, :),
            view(AICtr, j, :, :),
            body_vortex_panels.controlpoint,
            body_vortex_panels.normal,
            body_vortex_panels.tangent,
            rp.node,
            rp.nodemap,
            rp.influence_length,
        )
    end

    # - wake to body - #
    AICnw, AICtw = vortex_aic_boundary_on_field(
        body_vortex_panels.controlpoint,
        body_vortex_panels.normal,
        body_vortex_panels.tangent,
        wake_vortex_panels.node,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.influence_length,
    )

    # TODO: for the vz_xx's from here down, look into trying the @views macro to see if you can speed things up

    ##### ----- Induced Velcocities on Rotors ----- #####
    # - rotor to rotor - #

    v_rr = [
        induced_velocities_from_source_panels_on_points(
            rotor_source_panels[i].controlpoint,
            rotor_source_panels[j].node,
            rotor_source_panels[j].nodemap,
            rotor_source_panels[j].influence_length,
            ones(TF, 2, rotor_source_panels[j].totnode),
        ) for i in 1:length(rotor_source_panels), j in 1:length(rotor_source_panels)
    ]

    # axial components
    vz_rr = [
        v_rr[i, j][:, :, 1] for j in 1:length(rotor_source_panels),
        i in 1:length(rotor_source_panels)
    ]

    # radial components
    vr_rr = [
        v_rr[i, j][:, :, 2] for j in 1:length(rotor_source_panels),
        i in 1:length(rotor_source_panels)
    ]

    # - body to rotor - #

    v_rb = [
        induced_velocities_from_vortex_panels_on_points(
            rotor_source_panels[i].controlpoint,
            body_vortex_panels.node,
            body_vortex_panels.nodemap,
            body_vortex_panels.influence_length,
            ones(TF, 2, body_vortex_panels.totpanel),
        ) for i in 1:length(rotor_source_panels)
    ]

    # - body TE to rotor - #
    for i in 1:num_rotors
        induced_velocities_from_trailing_edge_gap_panel!(
            view(v_rb[i], :, :, :),
            rotor_source_panels[i].controlpoint,
            body_vortex_panels.tenode,
            body_vortex_panels.teinfluence_length,
            body_vortex_panels.tendotn,
            body_vortex_panels.tencrossn,
            body_vortex_panels.teadjnodeidxs,
        )
    end

    # axial components
    vz_rb = [v_rb[i][:, :, 1] for i in 1:length(rotor_source_panels)]

    # radial components
    vr_rb = [v_rb[i][:, :, 2] for i in 1:length(rotor_source_panels)]

    # - wake to rotor - #
    v_rw = [
        induced_velocities_from_vortex_panels_on_points(
            rotor_source_panels[i].controlpoint,
            wake_vortex_panels.node,
            wake_vortex_panels.nodemap,
            wake_vortex_panels.influence_length,
            ones(TF, 2, wake_vortex_panels.totpanel),
        ) for i in 1:length(rotor_source_panels)
    ]

    # axial components
    vz_rw = [v_rw[i][:, :, 1] for i in 1:length(rotor_source_panels)]

    # radial components
    vr_rw = [v_rw[i][:, :, 2] for i in 1:length(rotor_source_panels)]

    ##### ----- Induced Velocities on Wake ----- #####
    # - body to wake - #

    v_wb = induced_velocities_from_vortex_panels_on_points(
        wake_vortex_panels.controlpoint,
        body_vortex_panels.node,
        body_vortex_panels.nodemap,
        body_vortex_panels.influence_length,
        ones(TF, 2, body_vortex_panels.totpanel),
    )

    # - body TE to wake - #
    induced_velocities_from_trailing_edge_gap_panel!(
        view(v_wb, :, :, :),
        wake_vortex_panels.controlpoint,
        body_vortex_panels.tenode,
        body_vortex_panels.teinfluence_length,
        body_vortex_panels.tendotn,
        body_vortex_panels.tencrossn,
        body_vortex_panels.teadjnodeidxs,
    )

    # axial components
    vz_wb = v_wb[:, :, 1]

    # radial components
    vr_wb = v_wb[:, :, 2]

    # - rotor to wake - #
    v_wr = [
        induced_velocities_from_source_panels_on_points(
            wake_vortex_panels.controlpoint,
            rotor_source_panels[j].node,
            rotor_source_panels[j].nodemap,
            rotor_source_panels[j].influence_length,
            ones(TF, 2, rotor_source_panels[j].totpanel),
        ) for j in 1:length(rotor_source_panels)
    ]

    # axial components
    vz_wr = [v_wr[j][:, :, 1] for j in 1:length(rotor_source_panels)]

    # radial components
    vr_wr = [v_wr[j][:, :, 2] for j in 1:length(rotor_source_panels)]

    # - wake to wake - #
    v_ww = induced_velocities_from_vortex_panels_on_points(
        wake_vortex_panels.controlpoint,
        wake_vortex_panels.node,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.influence_length,
        ones(TF, 2, wake_vortex_panels.totpanel),
    )

    # axial components
    vz_ww = v_ww[:, :, 1]

    # radial components
    vr_ww = v_ww[:, :, 2]

    # TODO: you are here, but it seems like the rest of these are redundant and simply keeping track of indices in the wake would suffice.  Need to explore more before changing things though.

    # ##### ----- Induced Velcocities on Duct Wake ----- #####
    # # - body to duct wake - #
    # v_dwb = influencefromdoubletpanels(
    #     duct_wake_panels.controlpoint,
    #     body_vortex_panels.nodes,
    #     ones(TF, 2, body_vortex_panels.totpanel),
    # )

    # # axial components
    # vz_dwb = v_dwb[:, :, 1]

    # # - body TE to ductwake - #
    # v_dwbte = influencefromTE(
    #     duct_wake_panels.controlpoint,
    #     body_vortex_panels.TEnodes,
    #     ones(TF, 2, body_vortex_panels.totpanel),
    # )

    # # axial components
    # vz_dwbte = v_dwbte[:, :, 1]

    # # radial components
    # vr_dwbte = v_dwbte[:, :, 2]

    # # radial components
    # vr_dwb = v_dwb[:, :, 2]

    # # - rotor to duct wake - #
    # v_dwr = [
    #     influencefromsourcepanels(
    #         duct_wake_panels.controlpoint,
    #         rotor_source_panels[j].controlpoint,
    #         rotor_source_panels[j].len,
    #         ones(TF, 2, rotor_source_panels[j].npanels),
    #     ) for j in 1:length(rotor_source_panels)
    # ]

    # # axial components
    # vz_dwr = [v_dwr[j][:, :, 1] for j in 1:length(rotor_source_panels)]

    # # radial components
    # vr_dwr = [v_dwr[j][:, :, 2] for j in 1:length(rotor_source_panels)]

    # # - wake to duct wake - #
    # v_dww = influencefromvortexpanels(
    #     duct_wake_panels.controlpoint,
    #     wake_vortex_panels.controlpoint,
    #     wake_vortex_panels.len,
    #     ones(TF, 2, wake_vortex_panels.totpanel),
    # )

    # # axial components
    # vz_dww = v_dww[:, :, 1]

    # # radial components
    # vr_dww = v_dww[:, :, 2]

    # ##### ----- Induced Velcocities on Hub Wake ----- #####
    # # - body to duct wake - #
    # v_hwb = influencefromdoubletpanels(
    #     hub_wake_panels.controlpoint,
    #     body_vortex_panels.nodes,
    #     ones(TF, 2, body_vortex_panels.totpanel),
    # )

    # # axial components
    # vz_hwb = v_hwb[:, :, 1]

    # # radial components
    # vr_hwb = v_hwb[:, :, 2]

    # # - body TE to hubwake - #
    # v_hwbte = influencefromTE(
    #     hub_wake_panels.controlpoint,
    #     body_vortex_panels.TEnodes,
    #     ones(TF, 2, body_vortex_panels.totpanel),
    # )

    # # axial components
    # vz_hwbte = v_hwbte[:, :, 1]

    # # radial components
    # vr_hwbte = v_hwbte[:, :, 2]

    # # - rotor to duct wake - #
    # v_hwr = [
    #     influencefromsourcepanels(
    #         hub_wake_panels.controlpoint,
    #         rotor_source_panels[j].controlpoint,
    #         rotor_source_panels[j].len,
    #         ones(TF, 2, rotor_source_panels[j].npanels),
    #     ) for j in 1:length(rotor_source_panels)
    # ]

    # # axial components
    # vz_hwr = [v_hwr[j][:, :, 1] for j in 1:length(rotor_source_panels)]

    # # radial components
    # vr_hwr = [v_hwr[j][:, :, 2] for j in 1:length(rotor_source_panels)]

    # # - wake to duct wake - #
    # v_hww = influencefromvortexpanels(
    #     hub_wake_panels.controlpoint,
    #     wake_vortex_panels.controlpoint,
    #     wake_vortex_panels.len,
    #     ones(TF, 2, wake_vortex_panels.totpanel),
    # )

    # # axial components
    # vz_hww = v_hww[:, :, 1]

    # # radial components
    # vr_hww = v_hww[:, :, 2]

    #---------------------------------#
    #     Misc Values for Indexing    #
    #---------------------------------#
    # - Get rotor panel edges and centers - #
    rotor_panel_edges = [
        grid[2, rotor_indices_in_wake[i], 1:length(rpe)] for i in 1:num_rotors
    ]
    rotor_panel_centers = [rotor_source_panels[i].controlpoint[2, :] for i in 1:num_rotors]

    rotor_panel_edges = reduce(vcat, rotor_panel_edges)
    rotor_panel_centers = reduce(vcat, rotor_panel_centers)

    # get the total number of vortex panels on the bodies
    num_body_panels = body_vortex_panels.totpanel

    return (;
        converged=[false],
        Vconv=[0.0],
        iterations=[0],
        lu_decomp_flag,
        #freestream
        freestream,
        #reference for post process
        reference_parameters,
        # - Panels - #
        body_vortex_panels,
        rotor_source_panels,
        wake_vortex_panels,
        wakeK, #constant based on radial location that is used to define the wake panel strength
        duct_wake_panels,
        hub_wake_panels,
        # - rotors - #
        blade_elements, # blade elements
        num_rotors,
        rotor_panel_edges,
        rotor_panel_centers,
        # - Book Keeping - #
        num_body_panels,
        # ductTE_index=tip_gaps[1] == 0.0 ? ductTE_index : nothing,
        # hubTE_index=!nohub ? hubTE_index : nothing,
        ductidsaftofrotors,
        hubidsaftofrotors,
        ductwakeinterfaceid, # wake panel indices that lie on top of duct wall
        hubwakeinterfaceid, # wake panel indices that lie on top of hub wall
        rotorwakeid, # [rotor panel edge index, and closest forward rotor id] for each wake panel
        rotorwakepanelid, # [rotor panel index, and closest forward rotor id] for each wake panel
        num_wake_z_nodes, # number of nodes in axial direction of wake sheet
        num_wake_z_panels, # number of panels in each wake sheet
        # - Linear System - #
        # body_system_matrices, # includes the various LHS and RHS matrics and vectors for solving the linear system for the body
        # - Influence Matrices - #
        A_bb=LHS, # body to body
        AICt,
        b_bf=RHS, # freestream contribution to body boundary conditions
        RHS = similar(RHS).=0.0,
        gamb, # Body strengths
        A_br=AICnr, # rotor to body (total)
        AICtr,
        A_bw=AICnw, # wake to body (total)
        AICtw,
        v_rb, # body to rotor
        v_wb, # body to wake
        v_wr, # rotor to wake
        v_ww, # wake to wake
        # v_dwb, # body to duct wake
        # v_dwr, # rotor to duct wake
        # v_dww, # wake to duct wake
        # v_hwb, # body to hub wake
        # v_hwr, # rotor to hub wake
        # v_hww, # wake to hub wake
        # TODO: are both the separated out and full versions needed?  why not just one or the other?
        vz_rb, # body to rotor (x-direction)
        vr_rb, # body to rotor (r-direction)
        vz_rr, # rotor to rotor (x-direction)
        vr_rr, # rotor to rotor ( r-direction)
        vz_rw, # wake to rotor (x-direction)
        vr_rw, # wake to rotor ( r-direction)
        vz_wb, # body to wake (x-direction)
        vr_wb, # body to wake ( r-direction)
        vz_wr, # rotor to wake (x-direction)
        vr_wr, # rotor to wake ( r-direction)
        vz_ww, # wake to wake (x-direction)
        vr_ww, # wake to wake ( r-direction)
        # vz_dwb, # body to duct wake (x-direction)
        # vr_dwb, # body to duct wake ( r-direction)
        # vz_dwr, # rotor to duct wake (x-direction)
        # vr_dwr, # rotor to duct wake ( r-direction)
        # vz_dww, # wake to duct wake (x-direction)
        # vr_dww, # wake to duct wake ( r-direction)
        # vz_hwb, # body to hub wake (x-direction)
        # vr_hwb, # body to hub wake ( r-direction)
        # vz_hwr, # rotor to hub wake (x-direction)
        # vr_hwr, # rotor to hub wake ( r-direction)
        # vz_hww, # wake to hub wake (x-direction)
        # vr_hww, # wake to hub wake ( r-direction)
        # operating conditions
        Vinf=freestream.Vinf, # freestream parameters
        # - Debugging/Plotting
        duct_coordinates,
        hub_coordinates,
        isduct=!noduct,
        ishub=!nohub,
        grid=grid[1, :, 1:length(rpe)], #TODO: what is this used for, and why is it not the whole grid?
        TF, #floating point type to pass around
    )
end


function initialize_rotorwake_aero(inputs)

    # promoted type
    TF =  inputs.TF

    # Initialize
    Gamr = zeros(TF, length(inputs.blade_elements[1].rbe),length(inputs.blade_elements))
    sigr = zeros(TF, length(inputs.blade_elements[1].rbe)+1,length(inputs.blade_elements))
    gamw = zeros(TF, inputs.wake_vortex_panels.totnode)


 initialize_rotorwake_aero!(Gamr, sigr, gamw, inputs)

 return Gamr, sigr, gamw
 end


"""
"""
function initialize_rotorwake_aero!(Gamr, sigr, gamw, inputs)

    # - Rename for convenience - #
    freestream = inputs.freestream
    rotor_panel_centers = inputs.rotor_panel_centers

    # - inialize for re-use - #
    Vinf = similar(Gamr) .= inputs.freestream.Vinf # if desired, you can add the body induced velocities at each rotor control point to this array.
    vzind = zeros(eltype(Gamr), size(rotor_panel_centers, 1))
    vthetaind = zeros(eltype(Gamr), size(rotor_panel_centers, 1))
    Wm_dist = similar(Gamr, size(sigr,1)) .= 0
    Wm_wake = similar(Gamr, inputs.wake_vortex_panels.totpanel) .= 0

    for (irotor, (G, s, rpc, be, V)) in enumerate(
        zip(
            eachcol(Gamr),
            eachcol(sigr),
            eachcol(rotor_panel_centers),
            inputs.blade_elements,
            eachcol(Vinf),
        ),
    )

        # - Setup and Run CCBlade - #
        # define rotor, do not apply any corrections (including a tip correction)
        rotor = c4b.Rotor(be.Rhub, be.Rtip, be.B; tip=nothing)

        # define rotor sections
        sections = c4b.Section.(rpc, be.chords, be.twists, be.inner_airfoil)

        # define operating points using induced velocity from rotors ahead of this one
        op = [
            c4b.OperatingPoint(
                V[ir] .+ vz, # axial velocity V is freestream + body induced (maybe), vz is induced by rotor(s) ahead
                be.Omega * rpc[ir] .- vt, # tangential velocity #TODO: should this be + or - the induced velocity of the rotor ahead?
                freestream.rhoinf,
                0.0, #pitch is zero
                freestream.muinf,
                freestream.asound,
               ) for (ir, (vz, vt)) in enumerate(zip(vzind, vthetaind))
        ]

        # solve CCBlade problem for this rotor
        out = c4b.solve.(Ref(rotor), sections, op)

        # update induced velocities using far-field velocities
        vzind .+= out.u * 2.0
        vthetaind .+= out.v * 2.0

        # Get circulation starting point for this rotor
        G .= 0.5 .* out.cl .* out.W .* be.chords

        # Get source strength starting point (note we need the values at the nodes not centers, so average the values and use the end values on the end points)
        s[1] = @. 0.5 * out.cd[1] * out.W[1] * be.chords[1]
        @. s[2:(end - 1)] =
            (0.5 * out.cd[2:end] * out.W[2:end] * be.chords[2:end] +
            0.5 * out.cd[1:(end - 1)] * out.W[1:(end - 1)] * be.chords[1:(end - 1)])/2.0
        s[end] = @. 0.5 * out.cd[end] * out.W[end] * be.chords[end]

        # Get velocity distribution between this rotor and the next (same deal as with source panels, since wakes extend from source panel endpoints, we need to average velocities and use the ends for endpoints)
        Wm_dist[1]  =  sqrt((V[1] + vzind[1])^2 + vthetaind[1]^2)
        Wm_dist[2:(end - 1)] = (sqrt.((V[2:end] .+ vzind[2:end]).^2 .+ vthetaind[2:end].^2)
       .+sqrt.((V[2:end] .+ vzind[1:(end - 1)]).^2 .+ vthetaind[1:(end - 1)].^2))/2.0
        Wm_dist[end] = sqrt((V[end] + vzind[end])^2 + vthetaind[end]^2)

        inputs.Vconv[1] = sum(Wm_dist)/length(Wm_dist)

        # populate this section of the wake average velocities
        for (wid, wmap) in enumerate(eachrow(inputs.rotorwakepanelid))
            if wmap[2] >= irotor && wmap[2] < irotor + 1
                Wm_wake[wid] = Wm_dist[wmap[1]]
            end
        end
    end

    # - initialize wake strengths - #
    calculate_wake_vortex_strengths!(gamw, Gamr, Wm_wake, inputs)

    return Gamr, sigr, gamw
end
