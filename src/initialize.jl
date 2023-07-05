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

    ## -- Rename for Convenience -- ##
    # promoted type
    TF = promote_type(
        eltype(rotor_parameters[1].chords),
        eltype(rotor_parameters[1].twists),
        eltype(rotor_parameters[1].Omega),
        eltype(duct_coordinates),
        eltype(hub_coordinates),
        eltype(freestream.Vinf),
    )

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

    # get prescribed panel to be near duct leading edge
    _, leidx = findmin(t_duct_coordinates[:, 1])
    # prescribedpanels = [(leidx, 0.0)]
    prescribedpanels = [(1, 0.0)]

    # generate body paneling
    if nohub
        body_doublet_panels = generate_panels(t_duct_coordinates)
    else
        body_doublet_panels = generate_panels([t_duct_coordinates, rp_hub_coordinates])
    end

    # book keep wake wall interfaces
    # !NOTE: assumes duct geometry given first in paneling.  could generalize this, but not worth it for now.
    if nohub
        rotor_indices_on_hub = nothing
        hubidsaftofrotors = nothing
    else
        rotor_indices_on_hub = [
            findlast(
                x -> x < rotor_parameters.xrotor[i], body_doublet_panels.controlpoint[:, 1]
            ) for i in 1:length(rotor_parameters.xrotor)
        ]
        hwidraw = sort([body_doublet_panels.npanels; rotor_indices_on_hub])
        hubidsaftofrotors = reduce(
            vcat, [[(hwidraw[i] + 1):hwidraw[i + 1]] for i in 1:nrotor]
        )
    end

    if noduct
        rotor_indices_on_duct = nothing
        ductidsaftofrotors = nothing
    else
        rotor_indices_on_duct = [
            findfirst(
                x -> x < rotor_parameters.xrotor[i], body_doublet_panels.controlpoint[:, 1]
            ) - 1 for i in 1:length(rotor_parameters.xrotor)
        ]
        dwidraw = sort([0; rotor_indices_on_duct])
        ductidsaftofrotors = reverse(
            reduce(vcat, [[(dwidraw[i] + 1):dwidraw[i + 1]] for i in 1:nrotor])
        )
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

    # calculate radius dependent "constant" for wake strength calcualtion
    wakeK = get_wake_k(wake_vortex_panels)

    # Go through the wake panels and determine the index of the aftmost rotor infront and the blade node from which the wake strength is defined.
    rotorwakeid = ones(Int, paneling_constants.nwake_sheets, num_wake_x_panels, 2)
    rotorwakeid[:, :, 1] .= range(1, paneling_constants.nwake_sheets)
    for i in 1:num_wake_x_panels
        rotorwakeid[:, i, 2] .= findlast(
            x -> x < wake_vortex_panels.controlpoint[i, 1], rotor_parameters.xrotor
        )
    end
    rotorwakeid = reshape(
        rotorwakeid, (paneling_constants.nwake_sheets * num_wake_x_panels, 2)
    )

    # Go through the wake panels and determine the indices that have interfaces with the hub and wake
    hubwakeinterfaceid = 1:(hubTE_index - 1) #first rotor-wake-body interface is at index 1, this is also on the first wake sheet, so the hub trailing edge index in the xwake vector should be (or one away from, need to check) the last interface point
    ductwakeinterfaceid =
        num_wake_x_panels * (paneling_constants.nwake_sheets - 1) .+ (1:(ductTE_index - 1))

    # generate body wake panels for convenience (used in getting surface pressure of body wakes in post processing)
    duct_wake_panels = generate_panels(
        [xgrid[ductTE_index:end, end] rgrid[ductTE_index:end, end]]
    )
    hub_wake_panels = generate_panels([xgrid[hubTE_index:end, 1] rgrid[hubTE_index:end, 1]])

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
            rotor_source_panels[i].controlpoint[:, 2],
        ) for i in 1:nrotor
    ]

    #---------------------------------#
    # Calculate Coefficient Matrices  #
    #---------------------------------#

    ##### ----- Induced Velcocities on Bodies ----- #####

    # - body to body - #
    A_bb = doublet_panel_influence_matrix(body_doublet_panels.nodes, body_doublet_panels)

    # apply Kutta condition
    body_lhs_kutta!(A_bb, body_doublet_panels)

    # right hand side from freestream
    Vinfmat = repeat([freestream.Vinf 0], body_doublet_panels.npanels)
    b_bf = freestream_influence_vector(body_doublet_panels.normal, Vinfmat)

    # Set up matrices for body linear system
    Nsys = size(A_bb, 1)
    Npp = length(prescribedpanels)
    body_system_matrices = (;
        cache=false
        # LHS=A_bb,
        # LHSred=zeros(TF, Nsys, Nsys - Npp),
        # LHSlsq=zeros(TF, Nsys - Npp, Nsys - Npp),
        # RHS=zeros(TF, Nsys),
        # RHSvinf=b_bf,
        # RHSlsq=zeros(TF, Nsys - Npp),
        # prescribedpanels,
    )

    # - rotor to body - #
    A_br = [
        source_panel_influence_matrix(rotor_source_panels[j], body_doublet_panels) for
        j in 1:length(rotor_source_panels)
    ]

    # - wake to body - #
    A_bw = vortex_panel_influence_matrix(wake_vortex_panels, body_doublet_panels)

    ##### ----- Induced Velcocities on Rotors ----- #####
    # - rotor to rotor - #

    v_rr = [
        influencefromsourcepanels(
            rotor_source_panels[i].controlpoint,
            rotor_source_panels[j].controlpoint,
            rotor_source_panels[j].len,
            ones(TF, rotor_source_panels[j].npanels),
        ) for i in 1:length(rotor_source_panels), j in 1:length(rotor_source_panels)
    ]

    # axial components
    vx_rr = [
        v_rr[i, j][:, :, 1] for i in 1:length(rotor_source_panels),
        j in 1:length(rotor_source_panels)
    ]

    # radial components
    vr_rr = [
        v_rr[i, j][:, :, 2] for i in 1:length(rotor_source_panels),
        j in 1:length(rotor_source_panels)
    ]

    # - body to rotor - #

    v_rb = [
        influencefromdoubletpanels(
            rotor_source_panels[i].controlpoint,
            body_doublet_panels.nodes,
            ones(TF, body_doublet_panels.npanels),
        ) for i in 1:length(rotor_source_panels)
    ]

    # axial components
    vx_rb = [v_rb[i][:, :, 1] for i in 1:length(rotor_source_panels)]

    # radial components
    vr_rb = [v_rb[i][:, :, 2] for i in 1:length(rotor_source_panels)]

    # - body TE to rotor - #
    v_rbte = [
        influencefromTE(
            rotor_source_panels[i].controlpoint,
            body_doublet_panels.endpoints,
            body_doublet_panels.endpointidxs,
            ones(TF, body_doublet_panels.npanels),
        ) for i in 1:length(rotor_source_panels)
    ]

    # axial components
    vx_rbte = [v_rbte[i][:, :, 1] for i in 1:length(rotor_source_panels)]

    # radial components
    vr_rbte = [v_rbte[i][:, :, 2] for i in 1:length(rotor_source_panels)]

    # - wake to rotor - #
    v_rw = [
        influencefromvortexpanels(
            rotor_source_panels[i].controlpoint,
            wake_vortex_panels.controlpoint,
            wake_vortex_panels.len,
            ones(TF, wake_vortex_panels.npanels),
        ) for i in 1:length(rotor_source_panels), j in 1:length(wake_vortex_panels)
    ]

    # axial components
    vx_rw = [v_rw[i][:, :, 1] for i in 1:length(rotor_source_panels)]

    # radial components
    vr_rw = [v_rw[i][:, :, 2] for i in 1:length(rotor_source_panels)]

    ##### ----- Induced Velocities on Wake ----- #####
    # - body to wake - #

    v_wb = influencefromdoubletpanels(
        wake_vortex_panels.controlpoint,
        body_doublet_panels.nodes,
        ones(TF, body_doublet_panels.npanels),
    )

    # axial components
    vx_wb = v_wb[:, :, 1]

    # radial components
    vr_wb = v_wb[:, :, 2]

    # - body TE to wake - #
    v_wbte = influencefromTE(
        wake_vortex_panels.controlpoint,
        body_doublet_panels.endpoints,
        body_doublet_panels.endpointidxs,
        ones(TF, body_doublet_panels.npanels),
    )

    # axial components
    vx_wbte = v_wbte[:, :, 1]

    # radial components
    vr_wbte = v_wbte[:, :, 2]

    # - rotor to wake - #
    v_wr = [
        influencefromsourcepanels(
            wake_vortex_panels.controlpoint,
            rotor_source_panels[j].controlpoint,
            rotor_source_panels[j].len,
            ones(TF, rotor_source_panels[j].npanels),
        ) for j in 1:length(rotor_source_panels)
    ]

    # axial components
    vx_wr = [v_wr[j][:, :, 1] for j in 1:length(rotor_source_panels)]

    # radial components
    vr_wr = [v_wr[j][:, :, 2] for j in 1:length(rotor_source_panels)]

    # - wake to wake - #
    v_ww = influencefromvortexpanels(
        wake_vortex_panels.controlpoint,
        wake_vortex_panels.controlpoint,
        wake_vortex_panels.len,
        ones(TF, wake_vortex_panels.npanels),
    )

    # axial components
    vx_ww = v_ww[:, :, 1]

    # radial components
    vr_ww = v_ww[:, :, 2]

    ##### ----- Induced Velcocities on Duct Wake ----- #####
    # - body to duct wake - #
    A_dwb =
        v_dwb = influencefromdoubletpanels(
            duct_wake_panels.controlpoint,
            body_doublet_panels.nodes,
            ones(TF, body_doublet_panels.npanels),
        )

    # axial components
    vx_dwb = v_dwb[:, :, 1]

    # - body TE to ductwake - #
    v_dwbte = influencefromTE(
        duct_wake_panels.controlpoint,
        body_doublet_panels.endpoints,
        body_doublet_panels.endpointidxs,
        ones(TF, body_doublet_panels.npanels),
    )

    # axial components
    vx_dwbte = v_dwbte[:, :, 1]

    # radial components
    vr_dwbte = v_dwbte[:, :, 2]

    # radial components
    vr_dwb = v_dwb[:, :, 2]

    # - rotor to duct wake - #
    v_dwr = [
        influencefromsourcepanels(
            duct_wake_panels.controlpoint,
            rotor_source_panels[j].controlpoint,
            rotor_source_panels[j].len,
            ones(TF, rotor_source_panels[j].npanels),
        ) for j in 1:length(rotor_source_panels)
    ]

    # axial components
    vx_dwr = [v_dwr[j][:, :, 1] for j in 1:length(rotor_source_panels)]

    # radial components
    vr_dwr = [v_dwr[j][:, :, 2] for j in 1:length(rotor_source_panels)]

    # - wake to duct wake - #
    v_dww = influencefromvortexpanels(
        duct_wake_panels.controlpoint,
        wake_vortex_panels.controlpoint,
        wake_vortex_panels.len,
        ones(TF, wake_vortex_panels.npanels),
    )

    # axial components
    vx_dww = v_dww[:, :, 1]

    # radial components
    vr_dww = v_dww[:, :, 2]

    ##### ----- Induced Velcocities on Hub Wake ----- #####
    # - body to duct wake - #
    v_hwb = influencefromdoubletpanels(
        hub_wake_panels.controlpoint,
        body_doublet_panels.nodes,
        ones(TF, body_doublet_panels.npanels),
    )

    # axial components
    vx_hwb = v_hwb[:, :, 1]

    # radial components
    vr_hwb = v_hwb[:, :, 2]

    # - body TE to hubwake - #
    v_hwbte = influencefromTE(
        hub_wake_panels.controlpoint,
        body_doublet_panels.endpoints,
        body_doublet_panels.endpointidxs,
        ones(TF, body_doublet_panels.npanels),
    )

    # axial components
    vx_hwbte = v_hwbte[:, :, 1]

    # radial components
    vr_hwbte = v_hwbte[:, :, 2]

    # - rotor to duct wake - #
    v_hwr = [
        influencefromsourcepanels(
            hub_wake_panels.controlpoint,
            rotor_source_panels[j].controlpoint,
            rotor_source_panels[j].len,
            ones(TF, rotor_source_panels[j].npanels),
        ) for j in 1:length(rotor_source_panels)
    ]

    # axial components
    vx_hwr = [v_hwr[j][:, :, 1] for j in 1:length(rotor_source_panels)]

    # radial components
    vr_hwr = [v_hwr[j][:, :, 2] for j in 1:length(rotor_source_panels)]

    # - wake to duct wake - #
    v_hww = influencefromvortexpanels(
        hub_wake_panels.controlpoint,
        wake_vortex_panels.controlpoint,
        wake_vortex_panels.len,
        ones(TF, wake_vortex_panels.npanels),
    )

    # axial components
    vx_hww = v_hww[:, :, 1]

    # radial components
    vr_hww = v_hww[:, :, 2]

    #---------------------------------#
    #     Misc Values for Indexing    #
    #---------------------------------#
    # - Get rotor panel edges and centers - #
    rotor_panel_edges = [rgrid[rotor_indices_in_wake[i], 1:length(rpe)] for i in 1:nrotor]
    rotor_panel_centers = [rotor_source_panels[i].controlpoint[:, 2] for i in 1:nrotor]

    rotor_panel_edges = reduce(hcat, rotor_panel_edges)
    rotor_panel_centers = reduce(hcat, rotor_panel_centers)

    # get the total number of vortex panels on the bodies
    num_body_panels = body_doublet_panels.npanels

    return (;
        converged=[false],
        #freestream
        freestream,
        #reference for post process
        reference_parameters,
        # - Panels - #
        body_doublet_panels,
        rotor_source_panels,
        wake_vortex_panels,
        wakeK, #constant based on radial location that is used to define the wake panel strength
        duct_wake_panels,
        hub_wake_panels,
        # - rotors - #
        blade_elements, # blade elements
        num_rotors=nrotor,
        rotor_panel_edges,
        rotor_panel_centers,
        # - Book Keeping - #
        num_wake_x_panels, # number of wake panels in the axial direction
        num_body_panels,
        # ductTE_index=tip_gaps[1] == 0.0 ? ductTE_index : nothing,
        # hubTE_index=!nohub ? hubTE_index : nothing,
        ductidsaftofrotors,
        hubidsaftofrotors,
        ductwakeinterfaceid, # wake panel indices that lie on top of duct wall
        hubwakeinterfaceid, # wake panel indices that lie on top of hub wall
        rotorwakeid, # [rotor panel edge index, and closest forward rotor id] for each wake panel
        # - Linear System - #
        body_system_matrices, # includes the various LHS and RHS matrics and vectors for solving the linear system for the body
        # - Influence Matrices - #
        A_bb, # body to body
        b_bf, # freestream contribution to body boundary conditions
        prescribedpanels, # prescribed panels
        A_br, # rotor to body (total)
        A_bw, # wake to body (total)
        v_rb, # body to rotor
        v_rbte,
        v_wb, # body to wake
        v_wbte,
        v_wr, # rotor to wake
        v_ww, # wake to wake
        v_dwb, # body to duct wake
        v_dwr, # rotor to duct wake
        v_dww, # wake to duct wake
        v_hwb, # body to hub wake
        v_hwr, # rotor to hub wake
        v_hww, # wake to hub wake
        vx_rb, # body to rotor (x-direction)
        vr_rb, # body to rotor (r-direction)
        vx_rbte, # bodyTE to rotor (x-direction)
        vr_rbte, # bodyTE to rotor (r-direction)
        vx_rr, # rotor to rotor (x-direction)
        vr_rr, # rotor to rotor ( r-direction)
        vx_rw, # wake to rotor (x-direction)
        vr_rw, # wake to rotor ( r-direction)
        vx_wb, # body to wake (x-direction)
        vr_wb, # body to wake ( r-direction)
        vx_wbte, # bodyTE to wake (x-direction)
        vr_wbte, # bodyTE to wake ( r-direction)
        vx_wr, # rotor to wake (x-direction)
        vr_wr, # rotor to wake ( r-direction)
        vx_ww, # wake to wake (x-direction)
        vr_ww, # wake to wake ( r-direction)
        vx_dwb, # body to duct wake (x-direction)
        vr_dwb, # body to duct wake ( r-direction)
        vx_dwbte, # bodyTE to duct wake (x-direction)
        vr_dwbte, # bodyTE to duct wake ( r-direction)
        vx_dwr, # rotor to duct wake (x-direction)
        vr_dwr, # rotor to duct wake ( r-direction)
        vx_dww, # wake to duct wake (x-direction)
        vr_dww, # wake to duct wake ( r-direction)
        vx_hwb, # body to hub wake (x-direction)
        vr_hwb, # body to hub wake ( r-direction)
        vx_hwbte, # bodyTE to hub wake (x-direction)
        vr_hwbte, # bodyTE to hub wake ( r-direction)
        vx_hwr, # rotor to hub wake (x-direction)
        vr_hwr, # rotor to hub wake ( r-direction)
        vx_hww, # wake to hub wake (x-direction)
        vr_hww, # wake to hub wake ( r-direction)
        # operating conditions
        Vinf=freestream.Vinf, # freestream parameters
        # - Debugging/Plotting
        duct_coordinates=(noduct ? nothing : t_duct_coordinates),
        hub_coordinates=(nohub ? nothing : rp_hub_coordinates),
        isduct=!noduct,
        ishub=!nohub,
        wakexgrid=xgrid[:, 1:length(rpe)],
        wakergrid=rgrid[:, 1:length(rpe)],
    )
end

"""
    initialize_states(inputs)

Calculate an initial guess for the state variables
"""
function initialize_states(inputs)

    # get floating point type
    TF = promote_type(
        eltype(inputs.blade_elements[1].chords),
        eltype(inputs.blade_elements[1].twists),
        eltype(inputs.blade_elements[1].Omega),
        eltype(inputs.duct_coordinates),
    )

    # - Initialize body vortex strengths (rotor-off linear problem) - #
    # set up right hand side
    RHS = update_RHS(
        inputs.b_bf,
        inputs.A_bw,
        zeros(TF, inputs.wake_vortex_panels.npanels),
        inputs.A_br,
        zeros(TF, inputs.rotor_source_panels[1].npanels, 2),
    )

    # solve body-only problem
    mub = solve_body_strengths(inputs.A_bb, RHS, inputs.prescribedpanels)

    # - Initialize blade circulation and source strengths (assume open rotor) - #

    # get problem dimensions (number of radial stations x number of rotors)
    nr = length(inputs.blade_elements[1].rbe)
    nrotor = length(inputs.blade_elements)
    nwake = inputs.wake_vortex_panels.npanels

    # initialize outputs
    Gamr = zeros(TF, nr, nrotor)
    sigr = zeros(TF, nr, nrotor)
    gamw = zeros(TF, nwake)

    # initialize velocities on rotor blade elements
    _, _, _, _, Wtheta_rotor, Wm_rotor, Wmag_rotor = calculate_rotor_velocities(
        Gamr, gamw, sigr, mub, inputs
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
    Wm_wake = calculate_wake_velocities(gamw, sigr, mub, inputs)
    calculate_wake_vortex_strengths!(gamw, Gamr, Wm_wake, inputs)

    # - Combine initial states into one vector - #
    states = vcat(
        mub,                # body vortex panel strengths
        gamw,               # wake vortex sheet strengths
        reduce(vcat, Gamr), # rotor circulation strengths
        reduce(vcat, sigr), # rotor source panel strengths
    )

    return states
end
