"""
"""
function reinterpolate_geometry(
    problem_dimensions,
    duct_coordinates,
    centerbody_coordinates,
    rotorstator_parameters,
    paneling_constants;
    autoshiftduct=true,
    wake_nlsolve_ftol=1e-14,
    wake_max_iter=100,
    max_wake_relax_iter=3,
    wake_relax_tol=1e-14,
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
        wake_nlsolve_ftol=wake_nlsolve_ftol,
        wake_max_iter=wake_max_iter,
        max_wake_relax_iter=max_wake_relax_iter,
        wake_relax_tol=wake_relax_tol,
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
    wake_nlsolve_ftol=1e-14,
    wake_max_iter=100,
    max_wake_relax_iter=3,
    wake_relax_tol=1e-14,
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
        wake_nlsolve_ftol=wake_nlsolve_ftol,
        wake_max_iter=wake_max_iter,
        max_wake_relax_iter=max_wake_relax_iter,
        wake_relax_tol=wake_relax_tol,
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
    idmaps,
    rp_duct_coordinates,
    rp_centerbody_coordinates,
    rotorstator_parameters,
    paneling_constants,
    wake_grid;
    itcpshift=0.05,
    axistol=1e-15,
    tegaptol=1e1 * eps(),
    silence_warnings=true,
)

    # - Extract Tuples - #
    (; body_vortex_panels, rotor_source_panels, wake_vortex_panels) = panels
    (; npanels, ncenterbody_inlet, nduct_inlet, wake_length, nwake_sheets, dte_minus_cbte) =
        paneling_constants

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
function calculate_unit_induced_velocities(problem_dimensions, panels)
    (;
         nrotor,    # number of rotors
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
        v_rr=zeros(TF, nbe*nrotor, nws*nrotor, 2),
        v_rw=zeros(TF, nbe*nrotor, nwn, 2),
        v_rb=zeros(TF, nbe*nrotor, nbn, 2),
    )

    ivw = (;
        v_wr=zeros(TF, nwp, nws*nrotor, 2),
        v_ww=zeros(TF, nwp, nwn, 2),
        v_wb=zeros(TF, nwp, nbn, 2),
    )

    ivb = (;
        v_br=zeros(TF, nbp, nws*nrotor, 2),
        v_bw=zeros(TF, nbp, nwn, 2),
        v_bb=zeros(TF, nbp, nbn, 2),
    )

    return calculate_unit_induced_velocities!(ivr, ivw, ivb, panels)
end

"""
"""
function calculate_unit_induced_velocities!(ivr, ivw, ivb, panels)
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
        view(v_bb, :, :, :),
        body_vortex_panels.controlpoint,
        body_vortex_panels.node,
        body_vortex_panels.nodemap,
        body_vortex_panels.influence_length,
        ones(TF, 2, body_vortex_panels.totpanel[1]),
    )

    # Add influence of body trailing edge gap "panels"
    induced_velocities_from_trailing_edge_gap_panel!(
        view(v_bb, :, :, :),
        body_vortex_panels.controlpoint,
        body_vortex_panels.tenode,
        body_vortex_panels.teinfluence_length,
        body_vortex_panels.tendotn,
        body_vortex_panels.tencrossn,
        body_vortex_panels.teadjnodeidxs,
    )

    # - Rotors on Bodies - #
    # rotor panels to body panels
    induced_velocities_from_source_panels_on_points!(
        view(v_br, :, :, :),
        body_vortex_panels.controlpoint,
        rotor_source_panels.node,
        rotor_source_panels.nodemap,
        rotor_source_panels.influence_length,
        ones(TF, 2, rotor_source_panels.totnode[1]),
    )

    # - Wake on Bodies - #
    # wake panels to body panels
    induced_velocities_from_vortex_panels_on_points!(
        view(v_bw, :, :, :),
        body_vortex_panels.controlpoint,
        wake_vortex_panels.node,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.influence_length,
        ones(TF, 2, wake_vortex_panels.totpanel[1]),
    )

    # wake "TE panels" to body panels
    induced_velocities_from_trailing_edge_gap_panel!(
        view(v_bw, :, :, :),
        body_vortex_panels.controlpoint,
        wake_vortex_panels.tenode,
        wake_vortex_panels.teinfluence_length,
        wake_vortex_panels.tendotn,
        wake_vortex_panels.tencrossn,
        wake_vortex_panels.teadjnodeidxs,
    )

    ##### ----- Velocities on Rotors ----- #####
    # - Rotors on Rotors - #
    induced_velocities_from_source_panels_on_points!(
        view(v_rr, :, :, :),
        rotor_source_panels.controlpoint,
        rotor_source_panels.node,
        rotor_source_panels.nodemap,
        rotor_source_panels.influence_length,
        ones(TF, 2, rotor_source_panels.totnode[1]),
    )

    # - Bodies on Rotors - #
    # body panels on rotor panels
    induced_velocities_from_vortex_panels_on_points!(
        view(v_rb, :, :, :),
        rotor_source_panels.controlpoint,
        body_vortex_panels.node,
        body_vortex_panels.nodemap,
        body_vortex_panels.influence_length,
        ones(TF, 2, body_vortex_panels.totpanel[1]),
    )

    # add influence from body trailing edge gap "panels"
    induced_velocities_from_trailing_edge_gap_panel!(
        view(v_rb, :, :, :),
        rotor_source_panels.controlpoint,
        body_vortex_panels.tenode,
        body_vortex_panels.teinfluence_length,
        # -body_vortex_panels.teinfluence_length,
        body_vortex_panels.tendotn,
        body_vortex_panels.tencrossn,
        body_vortex_panels.teadjnodeidxs,
    )

    # - Wake on Rotors - #
    # wake panels on rotor panels
    induced_velocities_from_vortex_panels_on_points!(
        view(v_rw, :, :, :),
        rotor_source_panels.controlpoint,
        wake_vortex_panels.node,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.influence_length,
        ones(TF, 2, wake_vortex_panels.totpanel[1]),
    )

    # add influence from wake "trailing edge panels"
    induced_velocities_from_trailing_edge_gap_panel!(
        view(v_rw, :, :, :),
        rotor_source_panels.controlpoint,
        wake_vortex_panels.tenode,
        wake_vortex_panels.teinfluence_length,
        wake_vortex_panels.tendotn,
        wake_vortex_panels.tencrossn,
        wake_vortex_panels.teadjnodeidxs,
    )

    ##### ----- Velocities on Wakes ----- #####
    # - Bodies on Wakes - #
    # body panels to wake panels
    induced_velocities_from_vortex_panels_on_points!(
        view(v_wb, :, :, :),
        wake_vortex_panels.controlpoint,
        body_vortex_panels.node,
        body_vortex_panels.nodemap,
        body_vortex_panels.influence_length,
        ones(TF, 2, body_vortex_panels.totpanel),
    )

    # add influence from body trailing edge gap "panels"
    induced_velocities_from_trailing_edge_gap_panel!(
        view(v_wb, :, :, :),
        wake_vortex_panels.controlpoint,
        body_vortex_panels.tenode,
        body_vortex_panels.teinfluence_length,
        # -body_vortex_panels.teinfluence_length,
        body_vortex_panels.tendotn,
        body_vortex_panels.tencrossn,
        body_vortex_panels.teadjnodeidxs,
    )

    # - Rotors to Wakes - #
    induced_velocities_from_source_panels_on_points!(
        view(v_wr, :, :, :),
        wake_vortex_panels.controlpoint,
        rotor_source_panels.node,
        rotor_source_panels.nodemap,
        rotor_source_panels.influence_length,
        ones(TF, 2, rotor_source_panels.totpanel),
    )

    # - Wake on Wake - #
    # wake panels on wake panels
    induced_velocities_from_vortex_panels_on_points!(
        view(v_ww, :, :, :),
        wake_vortex_panels.controlpoint,
        wake_vortex_panels.node,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.influence_length,
        ones(TF, 2, wake_vortex_panels.totpanel),
    )

    # add influence from wake "trailing edge panels"
    induced_velocities_from_trailing_edge_gap_panel!(
        view(v_ww, :, :, :),
        wake_vortex_panels.controlpoint,
        wake_vortex_panels.tenode,
        wake_vortex_panels.teinfluence_length,
        wake_vortex_panels.tendotn,
        wake_vortex_panels.tencrossn,
        wake_vortex_panels.teadjnodeidxs,
    )

    return ivr, ivw, ivb
end

"""
"""
function initialize_linear_system(
    ivb, body_vortex_panels, rotor_source_panels, wake_vortex_panels, Vinf
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
        itcontrolpoint, itnormal, node, nodemap, influence_length
    )

    # Add Trailing Edge Gap Panel Influences to panels
    add_te_gap_aic!(
        @view(AICn[:, :]),
        controlpoint,
        normal,
        tenode,
        teinfluence_length,
        tendotn,
        tencrossn,
        teadjnodeidxs,
    )

    # Add Trailing Edge Gap Panel Influences to internal pseudo control point
    add_te_gap_aic!(
        @view(AICpcp[:, :]),
        itcontrolpoint,
        itnormal,
        tenode,
        teinfluence_length,
        tendotn,
        tencrossn,
        teadjnodeidxs,
    )

    # Assemble Raw LHS Matrix into A_bb
    A_bb = assemble_lhs_matrix(
        AICn,
        AICpcp,
        npanel,
        nnode,
        totpanel,
        totnode,
        prescribednodeidxs;
    )

    # - LU Decomposition - #
    A_bb_LU = lu!(A_bb, NoPivot(); check=false)

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

    A_pr = induced_velocities_from_source_panels_on_points(
        itcontrolpoint,
        rotor_source_panels.node,
        rotor_source_panels.nodemap,
        rotor_source_panels.influence_length,
        ones(eltype(b_bf),size(v_br,2))
    )

    ##### ----- Wake AIC ----- #####
    # - wake panels to body panels - #
    A_bw = calculate_normal_velocity(v_bw, normal)

    # add contributions from wake "trailing edge panels"
    add_te_gap_aic!(
        @view(A_bw[:, :]),
        controlpoint,
        normal,
        wake_vortex_panels.tenode,
        wake_vortex_panels.teinfluence_length,
        wake_vortex_panels.tendotn,
        wake_vortex_panels.tencrossn,
        wake_vortex_panels.teadjnodeidxs;
        wake=true,
    )

    # - wake panels on internal psuedo control point influence coefficients - #
    A_pw = vortex_aic_boundary_on_field(
        itcontrolpoint,
        itnormal,
        wake_vortex_panels.node,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.influence_length,
    )

    # add contributions from wake "trailing edge panels" on pseudo control point
    add_te_gap_aic!(
        @view(A_pw[:, :]),
        itcontrolpoint,
        itnormal,
        wake_vortex_panels.tenode,
        wake_vortex_panels.teinfluence_length,
        wake_vortex_panels.tendotn,
        wake_vortex_panels.tencrossn,
        wake_vortex_panels.teadjnodeidxs;
        wake=true,
    )

    return (; A_bb, A_bb_LU, lu_decomp_flag, b_bf, A_br, A_pr, A_bw, A_pw)
end
"""
"""
function initialize_linear_system!(
    linsys, ivb, body_vortex_panels, rotor_source_panels, wake_vortex_panels, Vinf, intermediate_containers
)
    # - Extract Tuples - #

    # containers for intermediate calcs
    (; AICn, AICpcp, vdnb, vdnpcp) = intermediate_containers

    # linear system
    (; A_bb, A_bb_LU, lu_decomp_flag, b_bf, A_br, A_pr, A_bw, A_pw) = linsys


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
    calculate_normal_velocity!(@view(AICn[:, :]), v_bb, normal)

    # Boundary on internal psuedo control point influence coefficients
    vortex_aic_boundary_on_field!(
        @view(AICpcp[:, :]), itcontrolpoint, itnormal, node, nodemap, influence_length
    )

    # Add Trailing Edge Gap Panel Influences to panels
    add_te_gap_aic!(
        @view(AICn[:, :]),
        controlpoint,
        normal,
        tenode,
        teinfluence_length,
        tendotn,
        tencrossn,
        teadjnodeidxs,
    )

    # Add Trailing Edge Gap Panel Influences to internal pseudo control point
    add_te_gap_aic!(
        @view(AICpcp[:, :]),
        itcontrolpoint,
        itnormal,
        tenode,
        teinfluence_length,
        tendotn,
        tencrossn,
        teadjnodeidxs,
    )

    # Assemble Raw LHS Matrix into A_bb
    assemble_lhs_matrix!(
        @view(A_bb[:, :]),
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
    A_bb_LU.factors .= lu!(A_bb, NoPivot(); check=false).factors # need to do it this way if A_bb_LU contains a view into a PreallocationTools DiffCache this allocates quite a bit more than just overwriting it though.
    # A_bb_LU = lu!(prelhs, NoPivot(); check=false) #this will overwrite A_bb_LU, and not update anything it was previously pointing to

    # check if success
    lu_decomp_flag[1] = eltype(v_bb)(issuccess(LHS))

    # - Freestream RHS - #
    vinfvec = [Vinf; 0.0]
    vdnb[:] .= [dot(vinfvec, nhat) for nhat in eachcol(normal)]
    vdnpcp[:] .= [dot(vinfvec, nhat) for nhat in eachcol(itnormal)]
    assemble_rhs_matrix!(
    @view(b_bf[:]), vdnb, vdnpcp, npanel, nnode, totpanel, totnode, prescribednodeidxs
    )

    ##### ----- Rotor AIC ----- #####
    calculate_normal_velocity!(@view(A_br[:, :]), v_br, normal)
    vortex_aic_boundary_on_field!(
        @view(A_pr[:, :]),
        itcontrolpoint,
        rotor_source_panels.node,
        rotor_source_panels.nodemap,
        rotor_source_panels.influence_length,
        ones(eltype(b_bf),size(v_br,2)),
    )

    ##### ----- Wake AIC ----- #####
    # - wake panels to body panels - #
    calculate_normal_velocity!(@view(A_bw[:, :]), v_bw, normal)

    # add contributions from wake "trailing edge panels"
    add_te_gap_aic!(
        @view(A_bw[:, :]),
        controlpoint,
        normal,
        wake_vortex_panels.tenode,
        wake_vortex_panels.teinfluence_length,
        wake_vortex_panels.tendotn,
        wake_vortex_panels.tencrossn,
        wake_vortex_panels.teadjnodeidxs;
        wake=true,
    )

    # - wake panels on internal psuedo control point influence coefficients - #
    vortex_aic_boundary_on_field!(
        @view(A_pw[:, :]),
        itcontrolpoint,
        wake_vortex_panels.node,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.influence_length,
        ones(eltype(b_bf),size(v_bw,2)),
    )

    # add contributions from wake "trailing edge panels" on pseudo control point
    add_te_gap_aic!(
        @view(A_pw[:, :]),
        itcontrolpoint,
        itnormal,
        wake_vortex_panels.tenode,
        wake_vortex_panels.teinfluence_length,
        wake_vortex_panels.tendotn,
        wake_vortex_panels.tencrossn,
        wake_vortex_panels.teadjnodeidxs;
        wake=true,
    )

    return linsys
end

function interpolate_blade_elements(rsp, Rtips, Rhubs, rotor_panel_centers, nbe)

    nrotor = length(rsp.B)
    Rtip = Rtips
    Rhub = Rhubs
    B = rsp.B
    fliplift = rsp.fliplift
    chords = similar(rsp.chords, nbe,nrotor) .= 0
    twists = similar(rsp.twists, nbe,nrotor) .= 0
    stagger = similar(rsp.twists, nbe,nrotor) .= 0
    solidity = similar(rsp.chords, nbe,nrotor) .= 0
    outer_airfoil = similar(rsp.airfoils, nbe,nrotor)
    inner_airfoil = similar(rsp.airfoils, nbe,nrotor)
    inner_fraction = similar(rsp.r, nbe,nrotor) .= 0

    for irotor in 1:nrotor
        rpcs = rotor_panel_centers[(nbe * (irotor - 1) + 1):(nbe * irotor)]

        # dimensionalize the blade element radial positions
        rblade = FLOWMath.linear([0.0; 1.0], [0.0; Rtip[irotor]], rsp.r[:, irotor])

        # update chord lengths
        chords[:, irotor] .= FLOWMath.akima(rblade, rsp.chords[:,irotor], rpcs)

        # update twists
        twists[:, irotor] .= FLOWMath.akima(rblade, rsp.twists[:,irotor], rpcs)

        # update stagger
        stagger[:, irotor] .= get_stagger(twists[:,irotor])

        # update solidity
        solidity[:, irotor] .= get_local_solidity(B[irotor], chords[:, irotor], rpcs)

        for ir in 1:nbe
            # outer airfoil
            io = min(length(rblade), searchsortedfirst(rblade, rpcs[ir]))
            outer_airfoil[ir,irotor] = rsp.airfoils[io, irotor]

            # inner airfoil
            ii = max(1, io - 1)
            inner_airfoil[ir,irotor] = rsp.airfoils[ii, irotor]

            # fraction of inner airfoil's polars to use
            if rblade[io] == rblade[ii]
                inner_fraction[ir, irotor] = 1.0
            else
                inner_fraction[ir, irotor] = (rpcs[ir] - rblade[ii]) / (rblade[io] - rblade[ii])
            end

            # Check incorrect extrapolation
            if inner_fraction[ir,irotor] > 1.0
                inner_fraction[ir,irotor] = 1.0
            end
        end
    end

    return (;
        Rtip,
        Rhub,
        rotor_panel_centers=reshape(rotor_panel_centers, (nbe, nrotor)),
        B,
        fliplift,
        chords,
        twists,
        stagger,
        solidity,
        outer_airfoil,
        inner_airfoil,
        inner_fraction,
    )
end

function interpolate_blade_elements!(blade_elements, rsp, Rtips, Rhubs, rotor_panel_center)
    # - Extract Blade Elements - #
    (;
        B,
        Rhub,
        Rtip,
        rotor_panel_centers,
        chords,
        twists,
        stagger,
        solidity,
        inner_airfoil,
        outer_airfoil,
        inner_fraction,
        fliplift,
    ) = blade_elements

    Rtip .= Rtips
    Rhub .= Rhubs
    B .= rsp.B
    fliplift .= rsp.fliplift
    rotor_panel_centers .= rotor_panel_center

    for irotor in 1:length(B)

        # dimensionalize the blade element radial positions
        rblade = FLOWMath.linear([0.0; 1.0], [0.0; Rtip[irotor]], rsp.r)

        # update chord lengths
        chords[:, irotor] .= FLOWMath.akima(rblade, chords, rotor_panel_centers)

        # update twists
        twists[:, irotor] .= FLOWMath.akima(rblade, twists, rotor_panel_centers)

        # update stagger
        stagger[:, irotor] .= get_stagger(twists)

        # update solidity
        solidity[:, irotor] .= get_local_solidity(B, chords, rotor_panel_centers)

        # get bounding airfoil polars
        outer_airfoil[:, irotor] .= similar(airfoils, length(rotor_panel_centers))
        inner_airfoil[:, irotor] .= similar(airfoils, length(rotor_panel_centers))
        inner_fraction[:, irotor] .= similar(airfoils, TF, length(rotor_panel_centers))

        for ir in 1:(length(rotor_panel_centers))
            # outer airfoil
            io = min(length(rblade), searchsortedfirst(rblade, rotor_panel_centers[ir]))
            outer_airfoil[ir] = airfoils[io]

            # inner airfoil
            ii = max(1, io - 1)
            inner_airfoil[ir] = airfoils[ii]

            # fraction of inner airfoil's polars to use
            if rblade[io] == rblade[ii]
                inner_fraction[ir] = 1.0
            else
                inner_fraction[ir] =
                    (rotor_panel_centers[ir] - rblade[ii]) / (rblade[io] - rblade[ii])
            end

            # Check incorrect extrapolation
            if inner_fraction[ir] > 1.0
                inner_fraction[ir] = 1.0
            end
        end
    end

    return blade_elements
end

"""
wnm = wake_vortex_panels.nodemap
wenids = wake_vortex_panels.endnodeidxs
nwn =  problem_dimensions.nwn
nwsn = problem_dimensions.nwsn
nbn = problem_dimensions.nbn
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
    npanels, ncenterbody_inlet, nwake_sheets, dte_minus_cbte, wnm, wenids, nwp, nwsp, nbn, riiw, nrotor
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
        cb_te_id = sum(npanels[1:(end - 1)]) + 1
    else
        cb_te_id = sum(npanels[1:(end - 2)]) + 1
    end
    wake_node_ids_along_centerbody_wake_interface = collect(range(1, cb_te_id; step=1))

    # indices of wake nodes interfacing with the casing wall
    if iszero(dte_minus_cbte) || dte_minus_cbte > 0
        duct_te_id = sum(npanels[1:(end - 1)]) + 1
    else
        duct_te_id = sum(npanels[1:(end - 2)]) + 1
    end

    ductteinwake = (sum(npanels) + 1) * nwake_sheets - (npanels[end] + 1)

    wake_node_ids_along_casing_wake_interface = collect(
        range(ductteinwake - duct_te_id, ductteinwake; step=1)
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
    else
        duct_te_id = sum(npanels[1:(end - 2)])
    end

    ductteinwake = (sum(npanels)) * nwake_sheets - (npanels[end])

    wake_panel_ids_along_casing_wake_interface = collect(
        range(ductteinwake - duct_te_id, ductteinwake; step=1)
    )

    # indices of duct panels interfacing with wake
    if dte_minus_cbte <= 0
        duct_panel_ids_along_centerbody_wake_interface = collect(
            1:sum([npanels[i] for i in (length(npanels) - 1):-1:1])
        )
    else
        dte_minus_cbte > 0
        duct_panel_ids_along_centerbody_wake_interface = collect(
            1:sum([npanels[i] for i in (length(npanels) - 2):-1:1])
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
            npanels[i] for i in (length(npanels) - 1):-1:1
        ])[2:end]
    elseif dte_minus_cbte > 0
        id_of_first_casing_panel_aft_of_each_rotor = cumsum([
            npanels[i] for i in (length(npanels) - 2):-1:1
        ])
    else
        id_of_first_casing_panel_aft_of_each_rotor = cumsum([
            npanels[i] for i in (length(npanels) - 1):-1:1
        ])
    end

    if iszero(dte_minus_cbte)
        id_of_first_centerbody_panel_aft_of_each_rotor = [ncenterbody_inlet+1; ncenterbody_inlet .+ [npanels[i] for i in 1:length(npanels)-2]]
    else
        id_of_first_centerbody_panel_aft_of_each_rotor = [ncenterbody_inlet+1; ncenterbody_inlet .+ [npanels[i] for i in 1:length(npanels)-3]]
    end

    # - Map of wake panel index to the wake sheet on which it resides and the last rotor ahead of the node - #
    wake_panel_sheet_be_map = ones(Int, nwp, 2)
    for i in 1:nwake_sheets
        wake_panel_sheet_be_map[(1 + (i - 1) * nwsp):(i * nwsp), 1] .= i
        for (ir, r) in enumerate(
            eachrow(@view(wake_panel_sheet_be_map[(1 + (i - 1) * nwsp):(i * nwsp), :]))
        )
            r[2] = min(nrotor, searchsortedlast(rotor_indices_in_wake, ir))
        end
    end

    return (;
        wake_nodemap,
        wake_endnodeidxs,
        wake_panel_sheet_be_map,
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
wnm = wake_vortex_panels.nodemap
venids = wake_vortex_panels.endnodeidxs
nbn = body_vortex_panels.totnode
wpsbm = wake_panel_sheet_be_map
e_map
    wake_panel_sheet_be_map = ones(Int, wake_vortex_panels.totnode, 2)
    num_wake_z_nodes = length(zwake)
    for i in 1:(rotorstator_parameters[1].nwake_sheets)
        wake_panel_sheet_be_map[(1 + (i - 1) * num_wake_z_nodes):(i * num_wake_z_nodes), 1] .= i
    end
    for (i, wn) in enumerate(eachcol(wake_vortex_panels.node))
        # TODO: DFDC geometry doesn't line up wake and rotor perfectly, so need a more robust option.
        # wake_panel_sheet_be_map[i, 2] = findlast(x -> x <= wn[1], rotorstator_parameters.rotorzloc)
        # TODO: current tee_mapsts are passing, but look here if things break in the future.
        wake_panel_sheet_be_map[i, 2] = findmin(x -> abs(x - wn[1]), rotorstator_parameters.rotorzloc)[2]
    end
"""
function set_index_maps!(
    idmaps, npanels, nwake_sheets, ncenterbody_inlet, dte_minus_cbte, wnm, wenids, nbn, wpsbm
)
    # - Extract Index Maps - #
    (;
        wake_nodemap,
        wake_endnodeidxs,
        wake_panel_sheet_be_map,
        wake_node_ids_along_casing_wake_interface,
        wake_node_ids_along_centerbody_wake_interface,
        id_of_first_casing_panel_aft_of_each_rotor,
        id_of_first_centerbody_panel_aft_of_each_rotor,
        rotor_indices_in_wake,
        body_totnodes,
    ) = idmaps

    wake_nodemap .= wnm
    wake_endnodeidxs .= wenids
    body_totnodes .= nbn

    if iszero(dte_minus_cbte) || dte_minus_cbte < 0
        cb_te_id = sum(npanels[1:(end - 1)]) + 1
    else
        cb_te_id = sum(npanels[1:(end - 2)]) + 1
    end
    wake_node_ids_along_centerbody_wake_interface .= collect(range(1, cb_te_id; step=1))

    if iszero(dte_minus_cbte) || dte_minus_cbte > 0
        duct_te_id = sum(npanels[1:(end - 1)]) + 1
    else
        duct_te_id = sum(npanels[1:(end - 2)]) + 1
    end

    ductteinwake = (sum(npanels) + 1) * nwake_sheets - (npanels[end] + 1)

    wake_node_ids_along_casing_wake_interface .= collect(
        range(ductteinwake - duct_te_id, ductteinwake; step=1)
    )

    for i in 1:nwake_sheets
        wake_panel_sheet_be_map[(1 + (i - 1) * nwsn):(i * nwsn), 1] .= i
        for (ir,r) in enumerate(eachrow(@view(wake_panel_sheet_be_map[(1 + (i - 1) * nwsn):(i * nwsn), :])))
            r[2] = min(nrotor,searchsortedlast(rotor_indices_in_wake, ir))
        end
    end

    if dte_minus_cbte < 0
        id_of_first_casing_panel_aft_of_each_rotor = cumsum([npanels[i] for i in length(npanels)-1:-1:1])[2:end]
    elseif dte_minus_cbte > 0
        id_of_first_casing_panel_aft_of_each_rotor = cumsum([npanels[i] for i in length(npanels)-2:-1:1])
    else
        id_of_first_casing_panel_aft_of_each_rotor = cumsum([npanels[i] for i in length(npanels)-1:-1:1])
    end

    if iszero(dte_minus_cbte)
        id_of_first_centerbody_panel_aft_of_each_rotor = [ncenterbody_inlet+1; ncenterbody_inlet .+ [npanels[i] for i in 1:length(npanels)-2]]
    else
        id_of_first_centerbody_panel_aft_of_each_rotor = [ncenterbody_inlet+1; ncenterbody_inlet .+ [npanels[i] for i in 1:length(npanels)-3]]
    end





    return idmaps
end

"""
"""
function precompute_parameters_iad(
    propulsor;
    wake_nlsolve_ftol=1e-14,
    wake_max_iter=100,
    max_wake_relax_iter=3,
    wake_relax_tol=1e-14,
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
        duct_coordinates,       # Matrix
        centerbody_coordinates, # Matrix
        rotorstator_parameters, # TODO: make this a named tuple of vectors
        paneling_constants,     # NamedTuple of scalars and vectors of scalars
        operating_point,        # NamedTuple of vectors (TODO) of numbers
        reference_parameters,   # NamedTuple of vectors (TODO) of numbers
    ) = propulsor

    problem_dimensions = get_problem_dimensions(paneling_constants)

    # - Reinterpolate Geometry and Generate Wake Grid - #
    wake_grid, rp_duct_coordinates, rp_centerbody_coordinates, rotor_indices_in_wake = reinterpolate_geometry(
        problem_dimensions,
        duct_coordinates,
        centerbody_coordinates,
        rotorstator_parameters,
        paneling_constants;
        wake_nlsolve_ftol=wake_nlsolve_ftol,
        wake_max_iter=wake_max_iter,
        max_wake_relax_iter=max_wake_relax_iter,
        wake_relax_tol=wake_relax_tol,
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

    # - Compute Influence Matrices - #
    ivr, ivw, ivb = calculate_unit_induced_velocities(
        problem_dimensions, (; body_vortex_panels, rotor_source_panels, wake_vortex_panels)
    )

    # - Set up Linear System - #
    linsys = initialize_linear_system(
        ivb, body_vortex_panels, rotor_source_panels, wake_vortex_panels, operating_point.Vinf
    )

    # - Interpolate Blade Elements - #
    blade_elements = interpolate_blade_elements(
        rotorstator_parameters, Rtips, Rhubs, rotor_source_panels.controlpoint[2, :], problem_dimensions.nbe
    )

    # - Get geometry-based constants for wake node strength calculations - #
    wakeK = get_wake_k(wake_vortex_panels)

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
        rotor_indices_in_wake,
        problem_dimensions.nrotor,
    )

    return ivr, ivw, ivb, linsys, blade_elements, wakeK, idmaps, (; body_vortex_panels, rotor_source_panels, wake_vortex_panels), problem_dimensions
end

"""
"""
function precompute_parameters_iad!(
    ivr,
    ivw,
    ivb,
    linsys,
    blade_elements,
    idmaps,
    panels,
    propulsor,
    precomp_containers; # contains wake_grid and repaneled duct and centerbody coordinates
    wake_nlsolve_ftol=1e-14,
    wake_max_iter=100,
    max_wake_relax_iter=3,
    wake_relax_tol=1e-14,
    itcpshift=0.05,
    axistol=1e-15,
    tegaptol=1e1 * eps(),
    finterp=FLOWMath.akima,
    silence_warnings=true,
    verbose=false,
)

    # - Extract propulsor - #
    (;
        duct_coordinates,       # Matrix
        centerbody_coordinates, # Matrix
        rotorstator_parameters, #  NamedTuple of vectors of numbers
        paneling_constants,     # NamedTuple of scalars and vectors of scalars
        operating_point,        # NamedTuple of vectors of numbers
        reference_parameters,   # NamedTuple of vectors of numbers
    ) = propulsor

    # - Extract Precomp Containers - #
    (; panels, wake_grid, rp_duct_coordinates, rp_centerbody_coordinates, AICn, AICpcp) =
        precomp_containers

    # - Reinterpolate Geometry and Generate Wake Grid - #
    # TODO: test this function
    reinterpolate_geometry!(
        precomp_containers,
        duct_coordinates,
        centerbody_coordinates,
        rotorstator_parameters,
        paneling_constants,
        idmaps.rotor_indices_in_wake;
        wake_nlsolve_ftol=wake_nlsolve_ftol,
        wake_max_iter=wake_max_iter,
        max_wake_relax_iter=max_wake_relax_iter,
        wake_relax_tol=wake_relax_tol,
        finterp=finterp,
        silence_warnings=silence_warnings,
    )

    # - Panel Everything - #
    # TODO: test this function
    generate_all_panels!(
        panels,
        idmaps,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        rotorstator_parameters,
        paneling_constants,
        wake_grid;
        itcpshift=itcpshift,
        axistol=axistol,
        tegaptol=tegaptol,
        silence_warnings=silence_warnings,
    )

    # - Compute Influence Matrices - #
    # TODO: test this function
    calculate_unit_induced_velocities!(ivr, ivw, ivb, panels)

    # - Set up Linear System - #
    # TODO: test this function
    initialize_linear_system!(linsys, ivb, panels.body_vortex_panels, AICn, AICpcp)

    # - Interpolate Blade Elements - #
    # TODO: test this function
    interpolate_blade_elements!(blade_elements, panels.rotor_source_panels)

    # - Save all the index mapping (bookkeeping) - #
    # TODO: test this function
    set_index_maps!(idmaps, paneling_constants)

    # - Get geometry-based constants for wake node strength calculations - #
    # TODO: write in-place version of this function
    get_wake_k!(wakeK, wake_vortex_panels)

    return ivr, ivw, ivb, linsys, blade_elements, wakeK, idmaps, panels
end

"""
"""
function initialize_velocities(
    operating_point, blade_elements, linsys, ivr, ivw, body_totnodes, wake_panel_sheet_be_map
)

    ##### ----- Initialize ----- #####
    # - get type - #
    #=
      NOTE: use anything in the operating point, the wake on body AIC should cover any body and wake geometry changes, and the rotor-on-rotor velocity should cover any rotor changes.
    =#
    TF = promote_type(
          eltype(operating_point.Omega),
          eltype(operating_point.Vinf),
          eltype(operating_point.rhoinf),
          eltype(operating_point.muinf),
          eltype(operating_point.asound),
          eltype(linsys.A_bw),
          eltype(ivr.v_rr)
         )

    # - Initialize velocity intermediates and outputs - #
    # rename to only compute once
    nbe = size(blade_elements.rotor_panel_centers)
    nbe_col = size(blade_elements.rotor_panel_centers, 1)

    # outputs
    Cm_wake = zeros(TF, size(ivw.v_ww, 1))
    vz_rotor = zeros(TF, nbe)
    vtheta_rotor = zeros(TF, nbe)

    # intermediate values
    sigr = zeros(TF, nbe[1] + 1, nbe[2])
    Cm_wake_vec = zeros(TF, nbe_col + 1)
    vthetaind = zeros(TF, nbe_col)
    vzind = zeros(TF, nbe_col)
    vrind = zeros(TF, nbe_col)

    # Solve Linear System for gamb
    gamb = copy(linsys.b_bf)
    ImplicitAD.implicit_linear(linsys.A_bb, gamb; lsolve=ldiv!, Af=linsys.A_bb_LU)

    # - Get body-induced velocities on rotors - #
    vzb = zeros(TF, nbe)
    vrb = zeros(TF, nbe)
    for irotor in 1:length(operating_point.Omega)
        berange = nbe_col*(irotor-1)+1:nbe_col*irotor
        vzb[:, irotor] = ivr.v_rb[berange,:,1] * gamb[1:body_totnodes]
        vrb[:, irotor] = ivr.v_rb[berange,:,2] * gamb[1:body_totnodes]
    end

    ##### ----- Loop through rotors ----- #####
    for irotor in 1:length(operating_point.Omega)

        #remove body influece for previous rotor and add it for this rotor
        if irotor > 1
            vzind .-= vzb[:, irotor - 1]
            vrind .-= vrb[:, irotor - 1]
        end
        vzind .+= vzb[:, irotor]
        vrind .+= vrb[:, irotor]

        # - Setup and Run CCBlade - #
        # define rotor, do not apply any corrections (including a tip correction)
        rotor = c4b.Rotor(
            blade_elements.Rhub[irotor],
            blade_elements.Rtip[irotor],
            blade_elements.B[irotor];
            tip=nothing,
        )

        # define rotor sections
        sections =
            c4b.Section.(
                blade_elements.rotor_panel_centers[:, irotor],
                blade_elements.chords[:, irotor],
                blade_elements.twists[:, irotor],
                blade_elements.inner_airfoil[:, irotor],
            )

        # define operating points using induced velocity from rotors ahead of this one
        c4bop = [
            c4b.OperatingPoint(
              operating_point.Vinf[1] + vz, # axial velocity V is freestream, vz is induced by bodies and rotor(s) ahead
                operating_point.Omega[irotor] *
                blade_elements.rotor_panel_centers[ir, irotor] + vt, # tangential velocity
                operating_point.rhoinf[1],
                0.0, #pitch is zero
                operating_point.muinf[1],
                operating_point.asound[1],
            ) for (ir, (vz, vt)) in enumerate(zip(vzind, vthetaind))
        ]

        # solve CCBlade problem for this rotor
        out = c4b.solve.(Ref(rotor), sections, c4bop)

        ##### ----- Assign Initial vz_rotor, vtheta_rotor, and Cm_wake ----- #####

        # -  vz_rotor and V_theta rotor - #
        # self influence
        vz_rotor[:, irotor] .+= vzind .+ out.u
        vtheta_rotor[:, irotor] .+= vthetaind .+ out.v

        # - Get Cm_wake - #
        #=
          NOTE: we are going to estimate this by taking the velocities on the rotors (though using far field z terms) and applying them constantly straight back to the next rotor or end of wake.
        =#
        # Get source strengths
        #=
          NOTE: we need the values at the nodes not centers, so average the values and use the end values on the end points
        =#
        sigr[1] = @. blade_elements.B[irotor] / (4.0 * pi) *
            out.cd[1] *
            out.W[1] *
            blade_elements.chords[1, irotor]
        @. sigr[2:(end - 1)] =
            (
                blade_elements.B[irotor] / (4.0 * pi) *
                out.cd[2:end] *
                out.W[2:end] *
                blade_elements.chords[2:end, irotor] +
                blade_elements.B[irotor] / (4.0 * pi) *
                out.cd[1:(end - 1)] *
                out.W[1:(end - 1)] *
                blade_elements.chords[1:(end - 1), irotor]
            ) / 2.0
        sigr[end] = @. blade_elements.B[irotor] / (4.0 * pi) *
            out.cd[end] *
            out.W[end] *
            blade_elements.chords[end, irotor]

        # add influence of rotor radial induced velocity from self and rotors ahead
        for jrotor in 1:irotor
            berange = nbe_col*(irotor-1)+1:nbe_col*irotor
            berangej = (nbe_col+1)*(jrotor-1)+1:(nbe_col+1)*jrotor
            vrind .+= ivr.v_rr[berange, berangej, 2] * sigr[:, jrotor]
        end

        # add in axial and tangential influence aft of current rotor
        vzind .+= 2.0 * out.u
        vthetaind .-= 2.0 * out.v

        # since wakes extend from source panel endpoints, we need to average velocities and use the ends for endpoints
        Cm_wake_vec[1] = sqrt((operating_point.Vinf[1] + vzind[1])^2 + vrind[1]^2)
        Cm_wake_vec[2:(end - 1)] =
            (
                sqrt.((operating_point.Vinf[1] .+ vzind[2:end]) .^ 2 .+ vrind[2:end] .^ 2) .+
                sqrt.((operating_point.Vinf[1] .+ vzind[1:(end - 1)]) .^ 2 .+ vrind[1:(end - 1)] .^ 2)
            ) / 2.0
        Cm_wake_vec[end] = sqrt((operating_point.Vinf[1] + vzind[end])^2 + vrind[end]^2)

        # fill in the section of the wake aft of the current rotor and up to the next rotor (or end of wake)
        for (wid, wmap) in enumerate(eachrow(wake_panel_sheet_be_map))
            if wmap[2] >= irotor && wmap[2] < irotor + 1
                Cm_wake[wid] = Cm_wake_vec[wmap[1]]
            end
        end
    end # loop through rotors

    return vz_rotor, vtheta_rotor, Cm_wake
end

"""
TODO; move this to the post-processing file
"""
function postcompute_parameters_aid!(cache, propulsor)
    # - Re-compute the precomputed stuff - #
    # - Compute the rest of the stuff required for post-processing - #
    return nothing
end
