"""
"""
function reinterpolate_geometry(
    problem_dimensions,
    duct_coordinates,
    centerbody_coordinates,
    rotorstator_parameters,
    paneling_constants;
    autoshiftduct=true,
    grid_solver_options=GridSolverOptions(),
    finterp=FLOWMath.akima,
    verbose=false,
    silence_warnings=true,
)

    # - Get Problem Dimensions - #
    (;
        nrotor, # number of rotors
        ndn,    # number of duct nodes
        ncbn,   # number of centerbody nodes
        nws,    # number of wake sheets (also rotor nodes)
        nwsn,   # number of nodes in each wake sheet
    ) = problem_dimensions

    # - Promote Type - #
    TF = promote_type(
        eltype(duct_coordinates),
        eltype(centerbody_coordinates),
        eltype(rotorstator_parameters.r),
        eltype(rotorstator_parameters.Rhub),
        eltype(rotorstator_parameters.Rtip),
        eltype(rotorstator_parameters.rotorzloc),
    )

    wake_grid = zeros(TF, 2, nwsn, nws)
    rp_duct_coordinates = zeros(TF, 2, ndn)
    rp_centerbody_coordinates = zeros(TF, 2, ncbn)
    rotor_indices_in_wake = ones(Int, nrotor)

    reinterpolate_geometry!(
        @view(wake_grid[:, :, :]),
        @view(rp_duct_coordinates[:, :]),
        @view(rp_centerbody_coordinates[:, :]),
        @view(rotor_indices_in_wake[:]),
        duct_coordinates,
        centerbody_coordinates,
        rotorstator_parameters,
        paneling_constants;
        autoshiftduct=autoshiftduct,
        grid_solver_options=grid_solver_options,
        finterp=finterp,
        verbose=verbose,
        silence_warnings=silence_warnings,
    )

    return wake_grid, rp_duct_coordinates, rp_centerbody_coordinates, rotor_indices_in_wake
end

"""
"""
function reinterpolate_geometry!(
    wake_grid,
    rp_duct_coordinates,
    rp_centerbody_coordinates,
    rotor_indices_in_wake,
    duct_coordinates,
    centerbody_coordinates,
    rotorstator_parameters,
    paneling_constants;
    autoshiftduct=true,
    grid_solver_options=GridSolverOptions(),
    finterp=FLOWMath.akima,
    verbose=false,
    silence_warnings=true,
)

    ##### ----- Extract Tuples ----- #####
    (; B, Rhub, Rtip, tip_gap, r, chords, twists, rotorzloc, airfoils, fliplift) =
        rotorstator_parameters
    (; npanels, ncenterbody_inlet, nduct_inlet, wake_length, nwake_sheets, dte_minus_cbte) =
        paneling_constants

    ##### ----- Re-interpolate bodies and rotors ----- #####

    # - Discretize Wake z-coordinates - #
    # also returns indices of rotor locations and duct and center body trailng edges in the wake
    zwake, rotor_indices_in_wake[:] = discretize_wake(
        duct_coordinates,
        centerbody_coordinates,
        rotorzloc, # rotor axial locations
        wake_length,
        npanels,
        dte_minus_cbte,
    )

    # - Re-interpolate Bodies - #
    reinterpolate_bodies!(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        duct_coordinates,
        centerbody_coordinates,
        zwake,
        ncenterbody_inlet,
        nduct_inlet;
        finterp=finterp,
    )

    # - Fix user errors in aft rotor tip gap definitions - #
    # for irotor in 2:num_rotors #TODO: use this version once tip gaps are functional on the first rotor
    for (irotor, tg) in enumerate(tip_gap)
        if tg != 0.0
            if !silence_warnings
                # @warn "DuctAPE does not currently have capabilities for adding tip gap to any but the foremost rotor. OverWRITING to 0.0."
                @warn "DuctAPE does not currently have capabilities for adding tip gaps to rotors. OverWRITING to 0.0."
            end
            tg = 0.0
        end
    end

    # - Move duct to correct position if user didn't provide coordintes with radial placement - #
    if autoshiftduct
        place_duct!(rp_duct_coordinates, Rtip[1], rotorzloc[1], tip_gap[1])
    end

    # - Fix any user errors in rotor radius definitons - #
    get_blade_ends_from_body_geometry!(
        Rtip,
        Rhub,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        tip_gap,
        rotorzloc;
        silence_warnings=silence_warnings,
    )

    ##### ----- Generate Wake Grid ----- #####

    generate_wake_grid!(
        wake_grid,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        Rhub[1],
        Rtip[1],
        tip_gap[1],
        zwake;
        grid_solver_options=grid_solver_options,
        verbose=verbose,
        silence_warnings=silence_warnings,
    )

    return rp_duct_coordinates, rp_centerbody_coordinates, wake_grid, rotor_indices_in_wake
end

"""
"""
function generate_all_panels(
    rp_duct_coordinates,
    rp_centerbody_coordinates,
    nwake_sheets,
    rotor_indices_in_wake,
    rotorzloc,
    wake_grid;
    itcpshift=0.05,
    axistol=1e-15,
    tegaptol=1e1 * eps(),
    silence_warnings=true,
)

    ##### ----- Fill Panel Objects ----- #####
    # - Body Panels - #
    body_vortex_panels = generate_panels([rp_duct_coordinates, rp_centerbody_coordinates])

    # - Rotor Panels - #
    #TODO: test this function
    rotor_source_panels = generate_rotor_panels(
        rotorzloc, wake_grid, rotor_indices_in_wake, nwake_sheets
    )

    # - Wake Panels - #
    # TODO: test this function
    wake_vortex_panels = generate_wake_panels(wake_grid[:, :, 1:nwake_sheets])

    #TODO; what other panels are actually needed? do you need the body wake panels or no?
    return (; body_vortex_panels, rotor_source_panels, wake_vortex_panels)
end

"""
"""
function generate_all_panels!(
    panels,
    wake_grid,
    rp_duct_coordinates,
    rp_centerbody_coordinates,
    rotor_indices_in_wake,
    rotorzloc,
    nwake_sheets;
    itcpshift=0.05,
    axistol=1e-15,
    tegaptol=1e1 * eps(),
    silence_warnings=true,
)

    # - Extract Tuples - #
    (; body_vortex_panels, rotor_source_panels, wake_vortex_panels) = panels

    ##### ----- Fill Panel Objects ----- #####
    # - Body Panels - #
    # TODO: test this function
    generate_panels!(body_vortex_panels, [rp_duct_coordinates, rp_centerbody_coordinates])

    # - Rotor Panels - #
    #TODO: test this function
    generate_rotor_panels!(
        rotor_source_panels, rotorzloc, wake_grid, rotor_indices_in_wake, nwake_sheets
    )

    # - Wake Panels - #
    # TODO: test this function
    generate_wake_panels!(wake_vortex_panels, wake_grid[:, :, 1:nwake_sheets])

    #TODO; what other panels are actually needed? do you need the body wake panels or no?
    return (; body_vortex_panels, rotor_source_panels, wake_vortex_panels)
end

"""
"""
function calculate_unit_induced_velocities(problem_dimensions, panels, integration_options)
    (;
        nrotor, # number of rotors
        nwn,    # number of wake nodes
        nwp,    # number of wake panels
        nbn,    # number of body nodes
        nbp,    # number of body panels
        nws,    # number of wake sheets (also rotor nodes)
        nbe,    # number of blade elements (also rotor panels)
    ) = problem_dimensions

    TF = promote_type(
        eltype(panels.body_vortex_panels.controlpoint),
        eltype(panels.rotor_source_panels.controlpoint),
        eltype(panels.wake_vortex_panels.controlpoint),
    )

    ivr = (;
        v_rr=zeros(TF, nbe * nrotor, nws * nrotor, 2),
        v_rw=zeros(TF, nbe * nrotor, nwn, 2),
        v_rb=zeros(TF, nbe * nrotor, nbn, 2),
    )

    ivw = (;
        v_wr=zeros(TF, nwp, nws * nrotor, 2),
        v_ww=zeros(TF, nwp, nwn, 2),
        v_wb=zeros(TF, nwp, nbn, 2),
    )

    ivb = (;
        v_br=zeros(TF, nbp, nws * nrotor, 2),
        v_bw=zeros(TF, nbp, nwn, 2),
        v_bb=zeros(TF, nbp, nbn, 2),
    )

    return calculate_unit_induced_velocities!(ivr, ivw, ivb, panels, integration_options)
end

"""
"""
function calculate_unit_induced_velocities!(ivr, ivw, ivb, panels, integration_options)
    # - Reset Tuples - #
    reset_containers!(ivr)
    reset_containers!(ivw)
    reset_containers!(ivb)

    # - Extract Tuples - #
    # Extract induced velocities on rotor
    (; v_rr, v_rw, v_rb) = ivr

    # Extract induced velocities on wake
    (; v_wr, v_ww, v_wb) = ivw

    # Extract induced velocities on body
    (; v_br, v_bw, v_bb) = ivb

    TF = promote_type(eltype.([ivr...])..., eltype.([ivw...])..., eltype.([ivb...])...)

    # extract panels
    (; body_vortex_panels, rotor_source_panels, wake_vortex_panels) = panels

    ##### ----- Velocities on Bodies ----- #####
    # - Bodies on Bodies - #
    # body panels on body panels
    induced_velocities_from_vortex_panels_on_points!(
        v_bb,
        body_vortex_panels.controlpoint,
        body_vortex_panels.node,
        body_vortex_panels.nodemap,
        body_vortex_panels.influence_length,
        #TODO: set things up differently so that you don't have to provide this allocated vector
        ones(TF, 2, Int(body_vortex_panels.totpanel[])),
        integration_options,
    )

    # Add influence of body trailing edge gap "panels"
    induced_velocities_from_trailing_edge_gap_panel!(
        v_bb,
        body_vortex_panels.controlpoint,
        body_vortex_panels.tenode,
        body_vortex_panels.teinfluence_length,
        body_vortex_panels.tendotn,
        body_vortex_panels.tencrossn,
        body_vortex_panels.teadjnodeidxs,
        integration_options,
    )

    # - Rotors on Bodies - #
    # rotor panels to body panels
    induced_velocities_from_source_panels_on_points!(
        v_br,
        body_vortex_panels.controlpoint,
        rotor_source_panels.node,
        rotor_source_panels.nodemap,
        rotor_source_panels.influence_length,
        ones(TF, 2, Int(rotor_source_panels.totnode[])),
        integration_options,
    )

    # - Wake on Bodies - #
    # wake panels to body panels
    induced_velocities_from_vortex_panels_on_points!(
        v_bw,
        body_vortex_panels.controlpoint,
        wake_vortex_panels.node,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.influence_length,
        ones(TF, 2, Int(wake_vortex_panels.totpanel[])),
        integration_options,
    )

    # wake "TE panels" to body panels
    induced_velocities_from_trailing_edge_gap_panel!(
        v_bw,
        body_vortex_panels.controlpoint,
        wake_vortex_panels.tenode,
        wake_vortex_panels.teinfluence_length,
        wake_vortex_panels.tendotn,
        wake_vortex_panels.tencrossn,
        wake_vortex_panels.teadjnodeidxs,
        integration_options;
        wake=true,
    )

    ##### ----- Velocities on Rotors ----- #####
    # - Rotors on Rotors - #
    induced_velocities_from_source_panels_on_points!(
        v_rr,
        rotor_source_panels.controlpoint,
        rotor_source_panels.node,
        rotor_source_panels.nodemap,
        rotor_source_panels.influence_length,
        ones(TF, 2, Int(rotor_source_panels.totnode[])),
        integration_options,
    )

    # - Bodies on Rotors - #
    # body panels on rotor panels
    induced_velocities_from_vortex_panels_on_points!(
        v_rb,
        rotor_source_panels.controlpoint,
        body_vortex_panels.node,
        body_vortex_panels.nodemap,
        body_vortex_panels.influence_length,
        ones(TF, 2, Int(body_vortex_panels.totpanel[])),
        integration_options,
    )

    # add influence from body trailing edge gap "panels"
    induced_velocities_from_trailing_edge_gap_panel!(
        v_rb,
        rotor_source_panels.controlpoint,
        body_vortex_panels.tenode,
        body_vortex_panels.teinfluence_length,
        body_vortex_panels.tendotn,
        body_vortex_panels.tencrossn,
        body_vortex_panels.teadjnodeidxs,
        integration_options,
    )

    # - Wake on Rotors - #
    # wake panels on rotor panels
    induced_velocities_from_vortex_panels_on_points!(
        v_rw,
        rotor_source_panels.controlpoint,
        wake_vortex_panels.node,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.influence_length,
        ones(TF, 2, Int(wake_vortex_panels.totpanel[])),
        integration_options,
    )

    # add influence from wake "trailing edge panels"
    induced_velocities_from_trailing_edge_gap_panel!(
        v_rw,
        rotor_source_panels.controlpoint,
        wake_vortex_panels.tenode,
        wake_vortex_panels.teinfluence_length,
        wake_vortex_panels.tendotn,
        wake_vortex_panels.tencrossn,
        wake_vortex_panels.teadjnodeidxs,
        integration_options;
        wake=true,
    )

    ##### ----- Velocities on Wakes ----- #####
    # - Bodies on Wakes - #
    # body panels to wake panels
    induced_velocities_from_vortex_panels_on_points!(
        v_wb,
        wake_vortex_panels.controlpoint,
        body_vortex_panels.node,
        body_vortex_panels.nodemap,
        body_vortex_panels.influence_length,
        ones(TF, 2, Int(body_vortex_panels.totpanel[])),
        integration_options,
    )

    # add influence from body trailing edge gap "panels"
    induced_velocities_from_trailing_edge_gap_panel!(
        v_wb,
        wake_vortex_panels.controlpoint,
        body_vortex_panels.tenode,
        body_vortex_panels.teinfluence_length,
        # -body_vortex_panels.teinfluence_length,
        body_vortex_panels.tendotn,
        body_vortex_panels.tencrossn,
        body_vortex_panels.teadjnodeidxs,
        integration_options,
    )

    # - Rotors to Wakes - #
    induced_velocities_from_source_panels_on_points!(
        v_wr,
        wake_vortex_panels.controlpoint,
        rotor_source_panels.node,
        rotor_source_panels.nodemap,
        rotor_source_panels.influence_length,
        ones(TF, 2, Int(rotor_source_panels.totpanel[])),
        integration_options,
    )

    # - Wake on Wake - #
    # wake panels on wake panels
    induced_velocities_from_vortex_panels_on_points!(
        v_ww,
        wake_vortex_panels.controlpoint,
        wake_vortex_panels.node,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.influence_length,
        ones(TF, 2, Int(wake_vortex_panels.totpanel[])),
        integration_options,
    )

    # add influence from wake "trailing edge panels"
    induced_velocities_from_trailing_edge_gap_panel!(
        v_ww,
        wake_vortex_panels.controlpoint,
        wake_vortex_panels.tenode,
        wake_vortex_panels.teinfluence_length,
        wake_vortex_panels.tendotn,
        wake_vortex_panels.tencrossn,
        wake_vortex_panels.teadjnodeidxs,
        integration_options;
        wake=true,
    )

    return ivr, ivw, ivb
end

"""
"""
function initialize_linear_system(
    ivb,
    body_vortex_panels,
    rotor_source_panels,
    wake_vortex_panels,
    Vinf,
    integration_options,
)

    # velocities on body
    (; v_br, v_bw, v_bb) = ivb

    # body panels
    (;
        npanel,
        nnode,
        totpanel,
        totnode,
        prescribednodeidxs,
        nodemap,
        node,
        tenode,
        influence_length,
        teinfluence_length,
        controlpoint,
        itcontrolpoint,
        normal,
        itnormal,
        tangent,
        tendotn,
        tencrossn,
        teadjnodeidxs,
    ) = body_vortex_panels

    ##### ----- Generate LHS ----- #####

    AICn = calculate_normal_velocity(v_bb, normal)

    # Boundary on internal psuedo control point influence coefficients
    AICpcp = vortex_aic_boundary_on_field(
        itcontrolpoint, itnormal, node, nodemap, influence_length, integration_options
    )

    # Add Trailing Edge Gap Panel Influences to internal pseudo control point
    add_te_gap_aic!(
        AICpcp,
        itcontrolpoint,
        itnormal,
        tenode,
        teinfluence_length,
        tendotn,
        tencrossn,
        teadjnodeidxs,
        integration_options,
    )

    # Assemble Raw LHS Matrix into A_bb
    A_bb = assemble_lhs_matrix(
        AICn, AICpcp, npanel, nnode, totpanel, totnode, prescribednodeidxs;
    )

    # - LU Decomposition - #
    A_bb_LU = lu!(copy(A_bb), NoPivot(); check=false)

    # check if success
    lu_decomp_flag = [eltype(v_bb)(issuccess(A_bb_LU))]

    # - Freestream RHS - #
    # vinf dot normal body
    vinfvec = [Vinf; 0.0]
    vdnb = [dot(vinfvec, nhat) for nhat in eachcol(normal)]
    vdnpcp = [dot(vinfvec, nhat) for nhat in eachcol(itnormal)]
    b_bf = assemble_rhs_matrix(
        vdnb, vdnpcp, npanel, nnode, totpanel, totnode, prescribednodeidxs
    )

    ##### ----- Rotor AIC ----- #####
    A_br = calculate_normal_velocity(v_br, normal)

    A_pr = source_aic(
        itcontrolpoint,
        itnormal,
        rotor_source_panels.node,
        rotor_source_panels.nodemap,
        rotor_source_panels.influence_length,
        integration_options,
    )

    ##### ----- Wake AIC ----- #####
    # - wake panels to body panels - #
    A_bw = calculate_normal_velocity(v_bw, normal)

    # - wake panels on internal psuedo control point influence coefficients - #
    A_pw = vortex_aic_boundary_on_field(
        itcontrolpoint,
        itnormal,
        wake_vortex_panels.node,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.influence_length,
        integration_options,
    )

    # add contributions from wake "trailing edge panels" on pseudo control point
    add_te_gap_aic!(
        A_pw,
        itcontrolpoint,
        itnormal,
        wake_vortex_panels.tenode,
        wake_vortex_panels.teinfluence_length,
        wake_vortex_panels.tendotn,
        wake_vortex_panels.tencrossn,
        wake_vortex_panels.teadjnodeidxs,
        integration_options;
        wake=true,
    )

    return (; A_bb, b_bf, A_br, A_pr, A_bw, A_pw), A_bb_LU, lu_decomp_flag
end

"""
"""
function initialize_linear_system!(
    linsys,
    ivb,
    body_vortex_panels,
    rotor_source_panels,
    wake_vortex_panels,
    Vinf,
    intermediate_containers,
    integration_options,
)

    # - Clear Containers - #
    reset_containers!(intermediate_containers)
    reset_containers!(linsys)

    # - Extract Tuples - #

    # containers for intermediate calcs
    (; AICn, AICpcp, vdnb, vdnpcp) = intermediate_containers

    # linear system
    (; A_bb, b_bf, A_br, A_pr, A_bw, A_pw) = linsys

    # velocities on body
    (; v_br, v_bw, v_bb) = ivb

    # body panels
    (;
        npanel,
        nnode,
        totpanel,
        totnode,
        prescribednodeidxs,
        nodemap,
        node,
        tenode,
        influence_length,
        teinfluence_length,
        controlpoint,
        itcontrolpoint,
        normal,
        itnormal,
        tangent,
        tendotn,
        tencrossn,
        teadjnodeidxs,
    ) = body_vortex_panels

    ##### ----- Generate LHS ----- #####

    # TODO officially test this function
    calculate_normal_velocity!(AICn, v_bb, normal)

    # Boundary on internal psuedo control point influence coefficients
    vortex_aic_boundary_on_field!(
        AICpcp,
        itcontrolpoint,
        itnormal,
        node,
        nodemap,
        influence_length,
        integration_options,
    )

    # Add Trailing Edge Gap Panel Influences to internal pseudo control point
    add_te_gap_aic!(
        AICpcp,
        itcontrolpoint,
        itnormal,
        tenode,
        teinfluence_length,
        tendotn,
        tencrossn,
        teadjnodeidxs,
        integration_options,
    )

    # Assemble Raw LHS Matrix into A_bb
    assemble_lhs_matrix!(
        A_bb,
        AICn,
        AICpcp,
        npanel,
        nnode,
        totpanel,
        totnode,
        prescribednodeidxs;
        dummyval=1.0,
    )

    # - LU Decomposition - #
    # TODO: you aren't going to be able to put the factorization inside the cache. it needs to always be a float. so you need to set it up like the airfoils, where it's allocated and passed in along side linsys.
    A_bb_LU = factorize_LHS(A_bb)
    lu_decomp_flag = eltype(A_bb)(issuccess(A_bb_LU))

    # - Freestream RHS - #
    vinfvec = [Vinf; 0.0]
    vdnb[:] .= [dot(vinfvec, nhat) for nhat in eachcol(normal)]
    vdnpcp[:] .= [dot(vinfvec, nhat) for nhat in eachcol(itnormal)]
    assemble_rhs_matrix!(
        b_bf, vdnb, vdnpcp, npanel, nnode, totpanel, totnode, prescribednodeidxs
    )

    ##### ----- Rotor AIC ----- #####
    calculate_normal_velocity!(A_br, v_br, normal)

    source_aic!(
        A_pr,
        itcontrolpoint,
        itnormal,
        rotor_source_panels.node,
        rotor_source_panels.nodemap,
        rotor_source_panels.influence_length,
        integration_options,
    )

    ##### ----- Wake AIC ----- #####
    # - wake panels to body panels - #
    calculate_normal_velocity!(A_bw, v_bw, normal)

    # - wake panels on internal psuedo control point influence coefficients - #
    vortex_aic_boundary_on_field!(
        A_pw,
        itcontrolpoint,
        itnormal,
        wake_vortex_panels.node,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.influence_length,
        integration_options,
    )

    # add contributions from wake "trailing edge panels" on pseudo control point
    add_te_gap_aic!(
        A_pw,
        itcontrolpoint,
        itnormal,
        wake_vortex_panels.tenode,
        wake_vortex_panels.teinfluence_length,
        wake_vortex_panels.tendotn,
        wake_vortex_panels.tencrossn,
        wake_vortex_panels.teadjnodeidxs,
        integration_options;
        wake=true,
    )

    return A_bb_LU, lu_decomp_flag
end

"""
wnm = wake_vortex_panels.nodemap
wenids = wake_vortex_panels.endnodeidxs
nwn =  problem_dimensions.nwn
nwsn = problem_dimensions.nwsn
nbn = problem_dimensions.nbn
ndp = body_vortex_panels.npanel[1]
riiw = rotor_indices_in_wake

    wake_panel_sheet_be_map = ones(Int, wake_vortex_panels.totnode, 2)
    num_wake_z_nodes = length(zwake)
    for i in 1:(rotorstator_parameters[1].nwake_sheets)
        wake_panel_sheet_be_map[(1 + (i - 1) * num_wake_z_nodes):(i * num_wake_z_nodes), 1] .= i
    end
    for (i, wn) in enumerate(eachcol(wake_vortex_panels.node))
        # TODO: DFDC geometry doesn't line up wake and rotor perfectly, so need a more robust option.
        # wake_panel_sheet_be_map[i, 2] = findlast(x -> x <= wn[1], rotorstator_parameters.rotorzloc)
        # TODO: current tests are passing, but look here if things break in the future.
        wake_panel_sheet_be_map[i, 2] = findmin(x -> abs(x - wn[1]), rotorstator_parameters.rotorzloc)[2]
    end
"""
function set_index_maps(
    npanels,
    ncenterbody_inlet,
    nwake_sheets,
    dte_minus_cbte,
    wnm,
    wenids,
    nwp,
    nwsp,
    nbn,
    ndp,
    riiw,
    nrotor,
)

    # - Quick Access Maps - #
    #=
      NOTE: used in various locations throughout
    =#
    # from wake_vortex_panels, map of wake nodes for each panel
    wake_nodemap = wnm
    # from wake_vortex_panels, get the node indices at the beginning and end of each wake sheet
    wake_endnodeidxs = wenids
    # from body_vortex_panels, the total number of body nodes
    body_totnodes = nbn
    # indices along each wake sheet at which the rotor(s) are placed.
    rotor_indices_in_wake = riiw

    # - Wake-Body Interface Node Maps - #
    #=
      NOTE: used for setting wake strengths along wake-body interfaces (node maps), or setting the proper jump velocity across the body boundaries (panel maps)
    =#
    # indices of wake nodes interfacing with the centerbody wall
    if iszero(dte_minus_cbte) || dte_minus_cbte < 0
        index_of_cb_te_along_wake_sheet = sum(npanels[1:(end - 1)]) + 1
    else
        index_of_cb_te_along_wake_sheet = sum(npanels[1:(end - 2)]) + 1
    end
    wake_node_ids_along_centerbody_wake_interface = collect(
        range(1, index_of_cb_te_along_wake_sheet; step=1)
    )

    # indices of wake nodes interfacing with the casing wall
    if iszero(dte_minus_cbte) || dte_minus_cbte > 0
        index_of_duct_te_along_wake_sheet = sum(npanels[1:(end - 1)]) + 1
        index_of_wake_node_at_duct_te = wake_endnodeidxs[end] - (npanels[end])
    else
        index_of_duct_te_along_wake_sheet = sum(npanels[1:(end - 2)]) + 1
        index_of_wake_node_at_duct_te = wake_endnodeidxs[end] - sum(npanels[(end - 1):end])
    end

    wake_node_ids_along_casing_wake_interface = collect(
        range(
            index_of_wake_node_at_duct_te - index_of_duct_te_along_wake_sheet + 1,
            index_of_wake_node_at_duct_te;
            step=1,
        ),
    )

    # indices of wake panels interfacing with centerbody wall
    if iszero(dte_minus_cbte) || dte_minus_cbte < 0
        cb_te_id = sum(npanels[1:(end - 1)])
    else
        cb_te_id = sum(npanels[1:(end - 2)])
    end
    wake_panel_ids_along_centerbody_wake_interface = collect(range(1, cb_te_id; step=1))

    # indices of wake panels interfacing with casing wall
    if iszero(dte_minus_cbte) || dte_minus_cbte > 0
        duct_te_id = sum(npanels[1:(end - 1)])
        ductteinwake = sum(npanels) * nwake_sheets - npanels[end]
    else
        duct_te_id = sum(npanels[1:(end - 2)])
        ductteinwake = sum(npanels) * nwake_sheets - npanels[end] - 1
    end

    wake_panel_ids_along_casing_wake_interface = collect(
        range(ductteinwake - duct_te_id + 1, ductteinwake; step=1)
    )

    # indices of duct panels interfacing with wake
    if dte_minus_cbte >= 0
        duct_panel_ids_along_centerbody_wake_interface = collect(
            sum([npanels[i] for i in 1:(length(npanels) - 1)]):-1:1
        )
    else
        dte_minus_cbte < 0
        duct_panel_ids_along_centerbody_wake_interface = collect(
            sum([npanels[i] for i in 1:(length(npanels) - 2)]):-1:1
        )
    end

    # indices of centerbody panels interfacing with wake
    if dte_minus_cbte <= 0
        centerbody_panel_ids_along_centerbody_wake_interface = collect(
            (ncenterbody_inlet + 1):(ncenterbody_inlet .+ sum([
                npanels[i] for i in 1:(length(npanels) - 1)
            ])),
        )
    else
        centerbody_panel_ids_along_centerbody_wake_interface = collect(
            (ncenterbody_inlet + 1):(ncenterbody_inlet .+ sum([
                npanels[i] for i in 1:(length(npanels) - 2)
            ])),
        )
    end

    # - Indices of the first body wall panels after a rotor - #
    #=
      NOTE: used in post-processing for augmenting the surface pressure due to pressure rise across rotors.
    =#
    if dte_minus_cbte < 0
        id_of_first_casing_panel_aft_of_each_rotor = cumsum([
            npanels[i] for i in (length(npanels) - 2):-1:1
        ])
    elseif dte_minus_cbte > 0
        id_of_first_casing_panel_aft_of_each_rotor = cumsum([
            npanels[i] for i in (length(npanels) - 1):-1:1
        ])[2:end]
    else
        id_of_first_casing_panel_aft_of_each_rotor = cumsum([
            npanels[i] for i in (length(npanels) - 1):-1:1
        ])
    end

    if iszero(dte_minus_cbte)
        id_of_first_centerbody_panel_aft_of_each_rotor = [
            ncenterbody_inlet + 1
            ncenterbody_inlet .+ cumsum([npanels[i] for i in 1:(length(npanels) - 2)])
        ]
    else
        id_of_first_centerbody_panel_aft_of_each_rotor = [
            ncenterbody_inlet + 1
            ncenterbody_inlet .+ cumsum([npanels[i] for i in 1:(length(npanels) - 3)])
        ]
    end
    id_of_first_centerbody_panel_aft_of_each_rotor .+= ndp #add on number of duct panels

    # - Map of wake panel index to the wake sheet on which it resides and the last rotor ahead of the panel - #
    wake_panel_sheet_be_map = ones(Int, nwp, 2)
    for i in 1:nwake_sheets
        wake_panel_sheet_be_map[(1 + (i - 1) * nwsp):(i * nwsp), 1] .= i
        for (ir, r) in enumerate(
            eachrow(@view(wake_panel_sheet_be_map[(1 + (i - 1) * nwsp):(i * nwsp), :]))
        )
            r[2] = min(nrotor, searchsortedlast(rotor_indices_in_wake, ir))
        end
    end

    # - Map of wake node index to the wake sheet on which it resides and the last rotor ahead of the node - #
    nwsn = nwsp + 1
    wake_node_sheet_be_map = ones(Int, Int(wenids[end]), 2)
    for i in 1:nwake_sheets
        wake_node_sheet_be_map[(1 + (i - 1) * nwsn):(i * nwsn), 1] .= i
        for (ir, r) in enumerate(
            eachrow(@view(wake_node_sheet_be_map[(1 + (i - 1) * nwsn):(i * nwsn), :]))
        )
            r[2] = min(nrotor, searchsortedlast(rotor_indices_in_wake, ir))
        end
    end

    return (;
        wake_nodemap,
        wake_endnodeidxs,
        wake_panel_sheet_be_map,
        wake_node_sheet_be_map,
        wake_node_ids_along_casing_wake_interface,
        wake_node_ids_along_centerbody_wake_interface,
        wake_panel_ids_along_casing_wake_interface,
        wake_panel_ids_along_centerbody_wake_interface,
        duct_panel_ids_along_centerbody_wake_interface,
        centerbody_panel_ids_along_centerbody_wake_interface,
        id_of_first_casing_panel_aft_of_each_rotor,
        id_of_first_centerbody_panel_aft_of_each_rotor,
        rotor_indices_in_wake,
        body_totnodes,
    )
end

"""
"""
function precompute_parameters(
    propulsor;
    grid_solver_options=GridSolverOptions(),
    integration_options=IntegrationOptions(),
    autoshiftduct=true,
    itcpshift=0.05,
    axistol=1e-15,
    tegaptol=1e1 * eps(),
    finterp=FLOWMath.akima,
    silence_warnings=true,
    verbose=false,
)

    # - Extract propulsor - #
    (;
        duct_coordinates,
        centerbody_coordinates,
        rotorstator_parameters,
        paneling_constants,
        operating_point,
    ) = propulsor

    problem_dimensions = get_problem_dimensions(paneling_constants)

    # - Reinterpolate Geometry and Generate Wake Grid - #
    wake_grid, rp_duct_coordinates, rp_centerbody_coordinates, rotor_indices_in_wake = reinterpolate_geometry(
        problem_dimensions,
        duct_coordinates,
        centerbody_coordinates,
        rotorstator_parameters,
        paneling_constants;
        grid_solver_options=grid_solver_options,
        autoshiftduct=autoshiftduct,
        finterp=finterp,
        verbose=verbose,
        silence_warnings=silence_warnings,
    )

    # shift duct into position
    if autoshiftduct
        place_duct!(
            rp_duct_coordinates,
            rotorstator_parameters.Rtip[1],
            rotorstator_parameters.rotorzloc[1],
            rotorstator_parameters.tip_gap[1],
        )
    end

    # Get actual Blade end positions from body geometry so there's not odd overlaps
    Rtips, Rhubs = get_blade_ends_from_body_geometry(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        rotorstator_parameters.tip_gap,
        rotorstator_parameters.rotorzloc,
    )

    return precompute_parameters(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        wake_grid,
        rotor_indices_in_wake,
        Rtips,
        Rhubs,
        rotorstator_parameters,
        paneling_constants,
        operating_point,
        integration_options,
        problem_dimensions;
        itcpshift=itcpshift,
        axistol=axistol,
        tegaptol=tegaptol,
        silence_warnings=silence_warnings,
        verbose=verbose,
    )
end

"""
"""
function precompute_parameters(
    rp_duct_coordinates,
    rp_centerbody_coordinates,
    wake_grid,
    rotor_indices_in_wake,
    Rtips,
    Rhubs,
    rotorstator_parameters,
    paneling_constants,
    operating_point,
    integration_options,
    problem_dimensions=nothing;
    itcpshift=0.05,
    axistol=1e-15,
    tegaptol=1e1 * eps(),
    silence_warnings=true,
    verbose=false,
)

    # - Panel Everything - #
    body_vortex_panels, rotor_source_panels, wake_vortex_panels = generate_all_panels(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        paneling_constants.nwake_sheets,
        rotor_indices_in_wake,
        rotorstator_parameters.rotorzloc,
        wake_grid;
        itcpshift=itcpshift,
        axistol=axistol,
        tegaptol=tegaptol,
        silence_warnings=silence_warnings,
    )

    # - Get problem dimensions if not already done - #
    if isnothing(problem_dimensions)
        problem_dimensions = get_problem_dimensions(
            body_vortex_panels, rotor_source_panels, wake_vortex_panels
        )
    end

    # - Compute Influence Matrices - #
    ivr, ivw, ivb = calculate_unit_induced_velocities(
        problem_dimensions,
        (; body_vortex_panels, rotor_source_panels, wake_vortex_panels),
        integration_options,
    )

    # - Set up Linear System - #
    linsys, A_bb_LU, lu_decomp_flag = initialize_linear_system(
        ivb,
        body_vortex_panels,
        rotor_source_panels,
        wake_vortex_panels,
        operating_point.Vinf,
        integration_options,
    )

    # - Interpolate Blade Elements - #
    blade_elements, airfoils = interpolate_blade_elements(
        rotorstator_parameters,
        Rtips,
        Rhubs,
        rotor_source_panels.controlpoint[2, :],
        problem_dimensions.nbe,
    )

    # - Get geometry-based constants for wake node strength calculations - #
    wakeK = get_wake_k(wake_vortex_panels.node[2, :], problem_dimensions.nwn)

    # - Save all the index mapping (bookkeeping) - #
    idmaps = set_index_maps(
        paneling_constants.npanels,
        paneling_constants.ncenterbody_inlet,
        paneling_constants.nwake_sheets,
        paneling_constants.dte_minus_cbte,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.endnodeidxs,
        problem_dimensions.nwp,
        problem_dimensions.nwsp,
        problem_dimensions.nbn,
        problem_dimensions.ndn - 1, #number of duct panels
        rotor_indices_in_wake,
        problem_dimensions.nrotor,
    )

    return ivr,
    ivw,
    ivb,
    linsys,
    A_bb_LU,
    blade_elements,
    airfoils,
    wakeK,
    idmaps,
    (; body_vortex_panels, rotor_source_panels, wake_vortex_panels),
    problem_dimensions
end

"""
"""
function precompute_parameters!(
    ivr,
    ivw,
    blade_element_cache,
    linsys,
    wakeK,
    propulsor,
    prepost_containers,
    problem_dimensions;
    grid_solver_options=GridSolverOptions(),
    integration_options=IntegrationOptions(),
    autoshiftduct=true,
    itcpshift=0.05,
    axistol=1e-15,
    tegaptol=1e1 * eps(),
    finterp=fm.akima,
    silence_warnings=true,
    verbose=false,
)

    # - unpack propulsor - #
    (;
        duct_coordinates,
        centerbody_coordinates,
        rotorstator_parameters,
        paneling_constants,
        operating_point,
    ) = propulsor

    # - Unpack preprocess containers - #
    (;
     wake_grid,
     rp_duct_coordinates,
     rp_centerbody_coordinates,
     rotor_indices_in_wake,
    ) = prepost_containers

    # - Reinterpolate Geometry and Generate Wake Grid - #
    # note: currently has 3526 allocations and takes 103.21 ms
    reinterpolate_geometry!(
        wake_grid,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        rotor_indices_in_wake,
        duct_coordinates,
        centerbody_coordinates,
        rotorstator_parameters,
        paneling_constants;
        autoshiftduct=autoshiftduct,
        grid_solver_options=grid_solver_options,
        finterp=finterp,
        verbose=verbose,
        silence_warnings=silence_warnings,
    )

    # shift duct into position
    if autoshiftduct
        place_duct!(
            rp_duct_coordinates,
            rotorstator_parameters.Rtip[1],
            rotorstator_parameters.rotorzloc[1],
            rotorstator_parameters.tip_gap[1],
        )
    end

    return precompute_parameters!(
        ivr,
        ivw,
        blade_element_cache,
        linsys,
        wakeK,
        wake_grid,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        rotor_indices_in_wake,
        rotorstator_parameters,
        paneling_constants,
        operating_point,
        prepost_containers,
        problem_dimensions;
        integration_options=integration_options,
        itcpshift=itcpshift,
        axistol=axistol,
        tegaptol=tegaptol,
        finterp=finterp,
        silence_warnings=silence_warnings,
        verbose=verbose,
    )
end

"""
"""
function precompute_parameters!(
    ivr,
    ivw,
    blade_element_cache,
    linsys,
    wakeK,
    wake_grid,
    rp_duct_coordinates,
    rp_centerbody_coordinates,
    rotor_indices_in_wake,
    rotorstator_parameters,
    paneling_constants,
    operating_point,
    prepost_containers,
    problem_dimensions=nothing;
    integration_options=IntegrationOptions(),
    itcpshift=0.05,
    axistol=1e-15,
    tegaptol=1e1 * eps(),
    finterp=fm.akima,
    silence_warnings=true,
    verbose=false,
)

    # - Reset Caches - #
    reset_containers!(ivr)
    reset_containers!(ivw)
    reset_containers!(blade_element_cache)
    reset_containers!(linsys)
    reset_containers!(wakeK)
    reset_containers!(prepost_containers; exception_keys=[:wake_grid,
     :rp_duct_coordinates,
     :rp_centerbody_coordinates,
     :rotor_indices_in_wake])

    # - Get Floating Point Type - #
    TF = promote_type(
        eltype(rp_duct_coordinates),
        eltype(rp_centerbody_coordinates),
        eltype(operating_point.Vinf),
        eltype(operating_point.Omega),
        eltype(operating_point.rhoinf),
        eltype(operating_point.muinf),
        eltype(operating_point.asound),
        eltype(rotorstator_parameters.B),
        eltype(rotorstator_parameters.Rhub),
        eltype(rotorstator_parameters.Rtip),
        eltype(rotorstator_parameters.rotorzloc),
        eltype(rotorstator_parameters.chords),
        eltype(rotorstator_parameters.twists),
    )

    # - Unpack preprocess containers - #
    (;
        panels,
        ivb,
        AICn,
        AICpcp,
        vdnb,
        vdnpcp,
    ) = prepost_containers

    # - Panel Everything - #
    generate_all_panels!(
        panels,
        wake_grid,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        rotor_indices_in_wake,
        rotorstator_parameters.rotorzloc,
        paneling_constants.nwake_sheets;
        itcpshift=itcpshift,
        axistol=axistol,
        tegaptol=tegaptol,
        silence_warnings=silence_warnings,
    )

    # - Get problem dimensions if not already done - #
    if isnothing(problem_dimensions)
        problem_dimensions = get_problem_dimensions(
            panels.body_vortex_panels, panels.rotor_source_panels, panels.wake_vortex_panels
        )
    end

    # - Compute Influence Matrices - #
    calculate_unit_induced_velocities!(ivr, ivw, ivb, panels, integration_options)

    # - Set up Linear System - #
    A_bb_LU, lu_decomp_flag = initialize_linear_system!(
        linsys,
        ivb,
        panels.body_vortex_panels,
        panels.rotor_source_panels,
        panels.wake_vortex_panels,
        operating_point.Vinf[1],
        (; AICn, AICpcp, vdnb, vdnpcp),
        integration_options,
    )

    # - Interpolate Blade Elements - #
    airfoils = interpolate_blade_elements!(
        blade_element_cache,
        rotorstator_parameters,
        panels.rotor_source_panels.controlpoint[2, :],
        problem_dimensions.nbe,
    )

    # - Get geometry-based constants for wake node strength calculations - #
    get_wake_k!(wakeK, panels.wake_vortex_panels.node[2, :])

    # - Save all the index mapping (bookkeeping) - #
    idmaps = set_index_maps(
        paneling_constants.npanels,
        paneling_constants.ncenterbody_inlet,
        paneling_constants.nwake_sheets,
        paneling_constants.dte_minus_cbte,
        Int.(panels.wake_vortex_panels.nodemap),
        Int.(panels.wake_vortex_panels.endnodeidxs),
        problem_dimensions.nwp,
        problem_dimensions.nwsp,
        problem_dimensions.nbn,
        problem_dimensions.ndn - 1, #number of duct panels
        rotor_indices_in_wake,
        problem_dimensions.nrotor,
    )

    return ivb, A_bb_LU, lu_decomp_flag, airfoils, idmaps, problem_dimensions
end
