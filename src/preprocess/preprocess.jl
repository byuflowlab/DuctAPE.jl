"""
    reinterpolate_geometry(
        problem_dimensions,
        duct_coordinates,
        centerbody_coordinates,
        rotor,
        paneling_constants;
        autoshiftduct=true,
        grid_solver_options=GridSolverOptions(),
        finterp=(x,y,xp)->FLOWMath.akima(x,y,xp,2.0*eps(),eps()),
        verbose=false,
        silence_warnings=true,
    )

Re-interpolate the body geometry and return compatible body and way geometry.

# Arguments
- `problem_dimensions::ProblemDimensions` : A ProblemDimensions object
- `duct_coordinates::Matrix{Float}` : [z,r] coordinates of duct geometry
- `centerbody_coordinates::Matrix{Float}` : [z,r] coordinates of centerbody geometry
- `rotor::Rotor` : A Rotor object
- `paneling_constants::PanelingConstants` : A PanelingConstants object

# Keyword Arguments
- `autoshiftduct::Bool=true` : flag to shift duct geometry based on rotor tip radius
- `grid_solver_options::SolverOptionsType=GridSolverOptions()` : options for the wake grid position solver
- `finterp::Function=FLOWMath.akima` : interpolation method for re-interpolating body coordinates
- `verbose::Bool=false` : flag to print verbose statements
- `silence_warnings::Bool=true` : flag to silence warnings

# Returns
- `wake_grid::Array{Float}` : array containig the z and r elliptic grid points defning the wake geometry.
- `rp_duct_coordinates::Matrix{Float}` : matrix containing the re-paneled duct coordinates
- `rp_centerbody_coordinates::Matrix{Float}` : matrix containing the re-paneled centerbody coordinates
- `rotor_indices_in_wake::Vector{Int}` : vector containing the indices of where in the wake the rotors reside (used later to define the rotor panel edges).
"""
function reinterpolate_geometry(
    problem_dimensions,
    duct_coordinates,
    centerbody_coordinates,
    rotor,
    paneling_constants;
    autoshiftduct=true,
    grid_solver_options=GridSolverOptions(),
    finterp=(x,y,xp)->FLOWMath.akima(x,y,xp,2.0*eps(),eps()),
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
        eltype(rotor.r),
        eltype(rotor.Rhub),
        eltype(rotor.Rtip),
        eltype(rotor.rotorzloc),
    )

    wake_grid = zeros(TF, 2, nwsn, nws)
    rp_duct_coordinates = zeros(TF, 2, ndn)
    rp_centerbody_coordinates = zeros(TF, 2, ncbn)
    rotor_indices_in_wake = ones(Int, nrotor)
    blade_element_cache = (;Rtip=zeros(TF,nrotor), Rhub=zeros(TF,nrotor))

    reinterpolate_geometry!(
        wake_grid,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        rotor_indices_in_wake,
        duct_coordinates,
        centerbody_coordinates,
        rotor,
        blade_element_cache,
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
    reinterpolate_geometry!(
        wake_grid,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        rotor_indices_in_wake,
        duct_coordinates,
        centerbody_coordinates,
        rotor,
        blade_element_cache,
        paneling_constants;
        autoshiftduct=true,
        grid_solver_options=GridSolverOptions(),
        finterp=(x,y,xp)->FLOWMath.akima(x,y,xp,2.0*eps(),eps()),
        verbose=false,
        silence_warnings=true,
    )

In-place version of `reinterpolate_geometry`.
"""
function reinterpolate_geometry!(
    wake_grid,
    rp_duct_coordinates,
    rp_centerbody_coordinates,
    rotor_indices_in_wake,
    duct_coordinates,
    centerbody_coordinates,
    rotor,
    blade_element_cache,
    paneling_constants;
    autoshiftduct=true,
    grid_solver_options=GridSolverOptions(),
    finterp=(x,y,xp)->FLOWMath.akima(x,y,xp,2.0*eps(),eps()),
    verbose=false,
    silence_warnings=true,
    le_bracket=1,
)

    ##### ----- Extract Tuples ----- #####
    (;Rhub, Rtip) = blade_element_cache

    (; B, tip_gap, r, chords, twists, rotorzloc, airfoils, is_stator) =
        rotor

    @assert length(unique(rotorzloc)) == length(rotorzloc) "Cannot place rotors on top of eachother: rotorzloc = $rotorzloc"

    Rhub .= rotor.Rhub
    Rtip .= rotor.Rtip

    (; npanels, ncenterbody_inlet, nduct_inlet, wake_length, nwake_sheets, dte_minus_cbte) =
        paneling_constants

    ##### ----- Re-interpolate bodies and rotors ----- #####

    # - Discretize Wake z-coordinates - #
    # also returns indices of rotor locations and duct and center body trailng edges in the wake
    zwake, rotor_indices_in_wake[:], duct_le_coordinates = discretize_wake(
        duct_coordinates,
        centerbody_coordinates,
        rotorzloc, # rotor axial locations
        wake_length,
        npanels,
        dte_minus_cbte;
        le_bracket=le_bracket
    )

    # - Re-interpolate Bodies - #
    reinterpolate_bodies!(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        duct_coordinates,
        centerbody_coordinates,
        zwake,
        duct_le_coordinates,
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
    generate_all_panels(
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

Function that calls all of the various panel generation functions are returns a named tuple containing all the panels

# Arguments
- `rp_duct_coordinates::Matrix{Float}` : matrix containing the re-paneled duct coordinates
- `rp_centerbody_coordinates::Matrix{Float}` : matrix containing the re-paneled centerbody coordinates
- `nwake_sheets::Int` : number of wake sheets
- `rotor_indices_in_wake::Vector{Int}` : vector containing the indices of where in the wake the rotors reside (used later to define the rotor panel edges).
- `rotorzloc:Vector{Float}` : axial locations of rotor lifting lines (contained in Rotor)
- `wake_grid::Array{Float}` : array containig the z and r elliptic grid points defning the wake geometry.

# Keyword Arguments
- `itcpshift::Float=0.05` : value used in positioning the internal pseudo control point in the solid bodies. Default is DFDC hard-coded value.
- `axistol::Float=1e-15` : tolerance for how close to the axis of rotation to be considered on the axis.
- `tegaptol::Float=1e1 * eps()` : tolerance for how large of a trailing edge gap is considered a gap.
- `silence_warnings::Bool=true` : flag to silence warnings

# Returns
- `panels::NamedTuple` : A named tuple of named tuples containing paneling information, including:
  - `body_vortex_panels::NamedTuple`
  - `rotor_source_panels::NamedTuple`
  - `wake_vortex_panels::NamedTuple`
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
    generate_all_panels!(
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

In-place version of `generate_all_panels`.
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
    calculate_unit_induced_velocities(problem_dimensions, panels, integration_options)

Calculate all the unit-induced velocties of all panels on all control points

# Arguments
- `problem_dimensions::ProblemDimensions` : A ProblemDimensions object
- `panels::NamedTuple` : A named tuple containing all the paneling information
- `integration_options::IntegrationOptions` : Options used for integration of velocity kernals across panels

# Returns
- `ivr::NamedTuple` : A named tuple containing arrays of induced velocities on the rotors
- `ivw::NamedTuple` : A named tuple containing arrays of induced velocities on the wake
- `ivb::NamedTuple` : A named tuple containing arrays of induced velocities on the bodies
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
    calculate_unit_induced_velocities!(ivr, ivw, ivb, panels, integration_options)

In-place version of `calculate_unit_induced_velocities`.
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
        integration_options,
    )

    # add influence from body trailing edge gap "panels"
    induced_velocities_from_trailing_edge_gap_panel!(
        v_wb,
        wake_vortex_panels.controlpoint,
        body_vortex_panels.tenode,
        body_vortex_panels.teinfluence_length,
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
    initialize_linear_system(
        ivb,
        body_vortex_panels,
        rotor_source_panels,
        wake_vortex_panels,
        Vinf,
        integration_options,
    )

Set up the linear system used in the panel method solve.

# Arguments
- `ivb::NamedTuple` : the named tuple containing all the unit induced velocities on the bodies
- `body_vortex_panels::NamedTuple` : the named tuple containing the body vortex panel information
- `rotor_source_panels::NamedTuple` : the named tuple containing the rotor source panel information
- `wake_vortex_panels::NamedTuple` : the named tuple containing the wake vortex panel information
- `Vinf::Vector{Float}` : the one-element vector containing the Freestream velocity magnitude
- `integration_options::IntegrationOptions` : the integration options used in integrating the panel induced velocities

# Returns
- `linsys::NamedTuple` : A named tuple containing cacheable data for the linear system, including:
  - `A_bb::Array{Float}` : AIC (LHS) matrix for the panel method system
  - `b_bf::Array{Float}` : Initial system RHS vector based on freestrem magnitude
  - `A_br::Array{Float}` : Unit normal velocity from rotors onto body panels
  - `A_pr::Array{Float}` : Unit normal velocity from rotors onto body internal psuedo control points
  - `A_bw::Array{Float}` : Unit normal velocity from wake onto body panels
  - `A_pw::Array{Float}` : Unit normal velocity from wake onto body internal psuedo control points
- `A_bb_LU::LinearAlgebra.LU` : LinearAlgebra LU factorization of the LHS matrix
- `lu_decomp_flag::Vector{Bool}` : flag for whether factorization was successful
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
       AICn, AICpcp, npanel, nnode, totpanel[], totnode[], prescribednodeidxs
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
       vdnb, vdnpcp, npanel, nnode, totpanel[], totnode[], prescribednodeidxs
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
    initialize_linear_system!(
        linsys,
        ivb,
        body_vortex_panels,
        rotor_source_panels,
        wake_vortex_panels,
        Vinf,
        intermediate_containers,
        integration_options,
    )

In-place version of `initialize_linear_system`.
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
    set_index_maps(
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

Set values for index map to be used throughout solve and post-process.

# Arguments
- `npanels : paneling_constants.npanels`
- `ncenterbody_inlet : paneling_constants.ncenterbody_inlet`
- `nwake_sheets : paneling_constants.nwake_sheets`
- `dte_minus_cbte : paneling_constants.dte_minus_cbte`
- `wnm : wake_vortex_panels.nodemap`
- `wenids : wake_vortex_panels.endnodeidxs`
- `nwp :  problem_dimensions.nwp`
- `nwsp : problem_dimensions.nwsp`
- `nbn : problem_dimensions.nbn`
- `ndp : body_vortex_panels.npanel[1]`
- `riiw : rotor_indices_in_wake`
- `nrotor : problem_dimensions.nrotor`

# Returns
- `idmaps::NamedTuple` : A named tuple containing index mapping used in bookkeeping throughout solve and post-process
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
    precompute_parameters(
        ducted_rotor,
        operating_point;
        grid_solver_options=GridSolverOptions(),
        integration_options=IntegrationOptions(),
        autoshiftduct=true,
        itcpshift=0.05,
        axistol=1e-15,
        tegaptol=1e1 * eps(),
        finterp=(x,y,xp)->FLOWMath.akima(x,y,xp,2.0*eps(),eps()),
        silence_warnings=true,
        verbose=false,
    )

Out of place main pre-processing function that computes all the required parameters for the solve.

# Arguments
- `ducted_rotor::DuctedRotor` : A DuctedRotor object
- `operating_point::OperatingPoint` : A OperatingPoint object

# Keyword Arguments
- `grid_solver_options::GridSolverOptionsType=GridSolverOptions()` : A GridSolverOptionsType object
- `integration_options::IntegrationMethod=IntegrationOptions()` : An IntegrationMethod object
- `autoshiftduct::Bool=true` : flag to shift duct geometry based on rotor tip radius
- `itcpshift::Float=0.05` : value used in positioning the internal pseudo control point in the solid bodies. Default is DFDC hard-coded value.
- `axistol::Float=1e-15` : tolerance for how close to the axis of rotation to be considered on the axis.
- `tegaptol::Float=1e1 * eps()` : tolerance for how large of a trailing edge gap is considered a gap.
- `finterp::Function=FLOWMath.akima` : interpolation method for re-interpolating body coordinates
- `silence_warnings::Bool=true` : flag to silence warnings
- `verbose::Bool=false` : flag to print verbose statements

# Returns
- `ivr::NamedTuple` : A named tuple containing arrays of induced velocities on the rotors
- `ivw::NamedTuple` : A named tuple containing arrays of induced velocities on the wake
- `ivb::NamedTuple` : A named tuple containing arrays of induced velocities on the bodies
- `linsys::NamedTuple` : A named tuple containing cacheable data for the linear system, including:
  - `A_bb::Array{Float}` : AIC (LHS) matrix for the panel method system
  - `b_bf::Array{Float}` : Initial system RHS vector based on freestrem magnitude
  - `A_br::Array{Float}` : Unit normal velocity from rotors onto body panels
  - `A_pr::Array{Float}` : Unit normal velocity from rotors onto body internal psuedo control points
  - `A_bw::Array{Float}` : Unit normal velocity from wake onto body panels
  - `A_pw::Array{Float}` : Unit normal velocity from wake onto body internal psuedo control points
- `A_bb_LU::LinearAlgebra.LU` : LinearAlgebra LU factorization of the LHS matrix
- `lu_decomp_flag::Vector{Bool}` : flag for whether factorization was successful
- `blade_elements::NamedTuple` : A named tuple containing cacheable blade element information (see docs for `interpolate_blade_elements`)
- `airfoils::Vector{AFType}` : A matrix of airfoil types associated with each of the blade elements
- `wakeK::Matrix{Float}` : A matrix of precomputed geometric constants used in the calculation of the wake vortex strengths
- `idmaps::NamedTuple` : A named tuple containing index mapping used in bookkeeping throughout solve and post-process
- `panels::NamedTuple` : A named tuple of panel objects including:
  - `body_vortex_panels::NamedTuple` : the named tuple containing the body vortex panel information
  - `rotor_source_panels::NamedTuple` : the named tuple containing the rotor source panel information
  - `wake_vortex_panels::NamedTuple` : the named tuple containing the wake vortex panel information
- `problem_dimensions::ProblemDimensions` : A ProblemDimensions object
"""
function precompute_parameters(
    ducted_rotor,
    operating_point;
    grid_solver_options=GridSolverOptions(),
    integration_options=IntegrationOptions(),
    autoshiftduct=true,
    itcpshift=0.05,
    axistol=1e-15,
    tegaptol=1e1 * eps(),
    finterp=(x,y,xp)->FLOWMath.akima(x,y,xp,2.0*eps(),eps()),
    silence_warnings=true,
    verbose=false,
)

    # - Extract ducted_rotor - #
    (;
        duct_coordinates,
        centerbody_coordinates,
        rotor,
        paneling_constants,
        operating_point,
    ) = ducted_rotor

    problem_dimensions = get_problem_dimensions(paneling_constants)

    # - Reinterpolate Geometry and Generate Wake Grid - #
    wake_grid, rp_duct_coordinates, rp_centerbody_coordinates, rotor_indices_in_wake = reinterpolate_geometry(
        problem_dimensions,
        duct_coordinates,
        centerbody_coordinates,
        rotor,
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
            rotor.Rtip[1],
            rotor.rotorzloc[1],
            rotor.tip_gap[1],
        )
    end

    # Get actual Blade end positions from body geometry so there's not odd overlaps
    Rtips, Rhubs = get_blade_ends_from_body_geometry(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        rotor.tip_gap,
        rotor.rotorzloc,
    )

    return precompute_parameters(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        wake_grid,
        rotor_indices_in_wake,
        Rtips,
        Rhubs,
        rotor,
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
    precompute_parameters(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        wake_grid,
        rotor_indices_in_wake,
        Rtips,
        Rhubs,
        rotor,
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

An alternate version of precompute_parameters allowing for user defined geometry that does not go through a re-panling step (use with caution).

The first inputs are the outputs of the `reinterpolate_geometry` and `get_blade_ends_from_body_geometry` functions.
"""
function precompute_parameters(
    rp_duct_coordinates,
    rp_centerbody_coordinates,
    wake_grid,
    rotor_indices_in_wake,
    Rtips,
    Rhubs,
    rotor,
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
        rotor.rotorzloc,
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
        operating_point.Vinf[],
        integration_options,
    )

    # - Interpolate Blade Elements - #
    blade_elements, airfoils = interpolate_blade_elements(
        rotor,
        Rtips,
        Rhubs,
        rotor_source_panels.controlpoint[2, :],
        problem_dimensions.nbe,
    )

    # - Get geometry-based constants for wake node strength calculations - #
    wakeK = get_wake_k(wake_vortex_panels.node[2, :])

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
    lu_decomp_flag,
    blade_elements,
    airfoils,
    wakeK,
    idmaps,
    (; body_vortex_panels, rotor_source_panels, wake_vortex_panels),
    problem_dimensions
end


"""
    precompute_parameters!(
        ivr,
        ivw,
        blade_element_cache,
        linsys,
        wakeK,
        ducted_rotor,
        operating_point,
        prepost_containers,
        problem_dimensions;
        grid_solver_options=GridSolverOptions(),
        integration_options=IntegrationOptions(),
        autoshiftduct=true,
        itcpshift=0.05,
        axistol=1e-15,
        tegaptol=1e1 * eps(),
        finterp=(x,y,xp)->FLOWMath.akima(x,y,xp,2.0*eps(),eps()),
        silence_warnings=true,
        verbose=false,
    )

In-place version of `precompute_parameters`.
"""
function precompute_parameters!(
    ivr,
    ivw,
    blade_element_cache,
    linsys,
    wakeK,
    ducted_rotor,
    operating_point,
    prepost_containers,
    problem_dimensions;
    grid_solver_options=GridSolverOptions(),
    integration_options=IntegrationOptions(),
    autoshiftduct=true,
    itcpshift=0.05,
    axistol=1e-15,
    tegaptol=1e1 * eps(),
    finterp=(x,y,xp)->FLOWMath.akima(x,y,xp,2.0*eps(),eps()),
    silence_warnings=true,
    verbose=false,
)

    # - unpack ducted_rotor - #
    (;
        duct_coordinates,
        centerbody_coordinates,
        rotor,
        paneling_constants,
    ) = ducted_rotor

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
        rotor,
        blade_element_cache,
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
            rotor.Rtip[1],
            rotor.rotorzloc[1],
            rotor.tip_gap[1],
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
        rotor,
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
    precompute_parameters!(
        ivr,
        ivw,
        blade_element_cache,
        linsys,
        wakeK,
        wake_grid,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        rotor_indices_in_wake,
        rotor,
        paneling_constants,
        operating_point,
        prepost_containers,
        problem_dimensions=nothing;
        integration_options=IntegrationOptions(),
        itcpshift=0.05,
        axistol=1e-15,
        tegaptol=1e1 * eps(),
        finterp=(x,y,xp)->FLOWMath.akima(x,y,xp,2.0*eps(),eps()),
        silence_warnings=true,
        verbose=false,
    )

In-place version of the `precompute_parameters` function by-passing the geometry reinterpolateion. (Use with caution)
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
    rotor,
    paneling_constants,
    operating_point,
    prepost_containers,
    problem_dimensions=nothing;
    integration_options=IntegrationOptions(),
    itcpshift=0.05,
    axistol=1e-15,
    tegaptol=1e1 * eps(),
    finterp=(x,y,xp)->FLOWMath.akima(x,y,xp,2.0*eps(),eps()),
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
        eltype(rotor.B),
        eltype(rotor.Rhub),
        eltype(rotor.Rtip),
        eltype(rotor.rotorzloc),
        eltype(rotor.chords),
        eltype(rotor.twists),
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
        rotor.rotorzloc,
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
        rotor,
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

    return A_bb_LU, lu_decomp_flag, airfoils, idmaps, problem_dimensions
end
