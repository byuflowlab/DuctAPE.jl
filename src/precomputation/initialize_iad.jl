"""
"""
function reinterpolate_geometry!(
    precomp_containers,
    duct_coordinates,
    centerbody_coordinates,
    rotorstator_parameters,
    paneling_constants,
    rotor_indices_in_wake;
    max_wake_relax_iter=100,
    wake_relax_tol=1e-9,
    finterp=fm.akima,
    silence_warnings=true,
)

    ##### ----- Extract Tuples ----- #####
    (; B, Rhub, Rtip, tip_gap, r, chords, twists, rotorzloc, airfoils, fliplift) =
        rotorstator_parameters
    (; npanels, ncenterbody_inlet, nduct_inlet, wake_length, nwake_sheets, dte_minus_cbte) =
        paneling_constants
    (; wake_grid, rp_duct_coordinates, rp_centerbody_coordinates) = precomp_containers

    ##### ----- Re-interpolate bodies and rotors ----- #####

    # - Discretize Wake z-coordinates - #
    # also returns indices of rotor locations and duct and center body trailng edges in the wake
    # TODO: update test for this function
    zwake, rotor_indices_in_wake = discretize_wake(
        duct_coordinates,
        centerbody_coordinates,
        rotorzloc, # rotor axial locations
        wake_length,
        npanels,
        dte_minus_cbte;
    )

    # - Re-interpolate Bodies - #
    # TODO: update test for this function
    # repanel_bodies!(
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
    # TODO: re-test the out of place version of this function
    get_blade_ends_from_body_geometry!(
        Rtip,
        Rhub,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        tip_gaps,
        rotorzloc;
        silence_warnings=silence_warnings,
    )

    ##### ----- Generate Wake Grid ----- #####

    generate_wake_grid!(
        wake_grid,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        Rhub,
        Rtip,
        nwake_sheets;
        max_wake_relax_iter=max_wake_relax_iter,
        wake_relax_tol=wake_relax_tol,
        verbose=verbose,
        silence_warnings=silence_warnings,
    )

    return rp_duct_coordinates, rp_centerbody_coordinates, wake_grid, rotor_indices_in_wake
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
function calculate_unit_induced_velocities!(ivr, ivw, ivb, panels)
    # - Extract induced velocities on rotor - #
    (; v_rr, v_rr, v_rw, v_rw, v_rb, v_rb) = ivr

    # - Extract induced velocities on wake - #
    (; v_wr, v_wr, v_ww, v_ww, v_wb, v_wb) = ivw

    # - Extract induced velocities on body - #
    (; v_br, v_br, v_bw, v_bw, v_bb, v_bb) = ivb

    return nothing
end

"""
"""
function initialize_linear_system(linsys, ivb, body_vortex_panels)
    # - Extract Linear System - #
    (; A_bb, A_bb_LU, lu_decomp_flag, b_bf, A_br, A_pr, A_bw, A_pw) = linsys

    return linsys
end

"""
"""
function interpolate_blade_elements!(blade_elements, rotor_source_panels)
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

    return nothing
end

"""
"""
function set_index_maps!(idmaps, paneling_constants)
    # - Extract Index Maps - #
    (;) = idmaps
    return nothing
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
    max_wake_relax_iter=100,
    wake_relax_tol=1e-9,
    itcpshift=0.05,
    axistol=1e-15,
    tegaptol=1e1 * eps(),
    finterp=fm.akima,
    silence_warnings=true,
)

    # - Extract propulsor - #
    (;
        duct_coordinates, # Matrix
        centerbody_coordinates, # Matrix
        rotorstator_parameters, # Vector of NamedTuples of a bunch of stuff...
        paneling_constants, # NamedTuple of numbers and vectors of numbers
        operating_point, # NamedTuple of numbers
        reference_parameters, # NamedTuple of numbers
    ) = propulsor

    # - Reinterpolate Geometry and Generate Wake Grid - #
    # TODO: test this function
    reinterpolate_geometry!(
        precomp_containers,
        duct_coordinates,
        centerbody_coordinates,
        rotorstator_parameters,
        paneling_constants,
        idmaps.rotor_indices_in_wake;
        max_wake_relax_iter=max_wake_relax_iter,
        wake_relax_tol=wake_relax_tol,
        finterp=finterp,
        silence_warnings=silence_warnings,
    )

    # - Panel Everything - #
    # TODO: write a function that just does all the paneling
    generate_all_panels!(
        panels,
        idmaps,
        precomp_containers.rp_duct_coordinates,
        precomp_containers.rp_centerbody_coordinates,
        rotorstator_parameters,
        paneling_constants,
        precomp_containers.wake_grid;
        itcpshift=0.05,
        axistol=1e-15,
        tegaptol=1e1 * eps(),
        silence_warnings=silence_warnings,
    )

    # - Compute Influence Matrices - #
    # TODO: write a function that does all the influence matrix stuff
    calculate_unit_induced_velocities!(ivr, ivw, ivb, panels)

    # - Set up Linear System - #
    # TODO: write function that assembles the linear system
    initialize_linear_system!(linsys, ivb, panels.body_vortex_panels)

    # - Interpolate Blade Elements - #
    # TODO: write function that interpolates the blade elements
    interpolate_blade_elements!(blade_elements, panels.rotor_source_panels)

    # - Save all the index mapping (bookkeeping) - #
    # TODO: write a function for any additional index mapping that didn't get added in the above functions
    set_index_maps!(idmaps, paneling_constants)

    return ivr, ivw, ivb, linsys, blade_elements, idmaps, panels
end

"""
"""
function initialize_velocities(
    op, blade_elements, linsys, ivr, ivw, nbodynodes, rotorwakenodeid
)

    ##### ----- Initialize ----- #####
    # - get type - #
    # note: use anything in the operating point, the wake on body AIC should cover any body and wake geometry changes, and the rotor-on-rotor velocity should cover any rotor changes.
    TF = promote_type(eltype.([op...])..., eltype(linsys.A_bw), eltype(ivr.vz_rr))

    # - Initialize velocity intermediates and outputs - #
    # rename to only compute once
    nbe = size(blade_elements.rotor_panel_centers)
    nbe_col = size(blade_elements.rotor_panel_centers, 1)

    # outputs
    Cm_wake = zeros(TF, size(ivw.vz_ww, 1))
    Vz_rotor = zeros(TF, nbe)
    Vtheta_rotor = zeros(TF, nbe)

    # intermediate values
    sigr = zeros(TF, nbe[1] + 1, nbe[2])
    Cm_wake_vec = zeros(TF, nbe_col + 1)
    vthetaind = zeros(TF, nbe_col)
    vzind = zeros(TF, nbe_col)
    vrind = zeros(TF, nbe_col)

    # Solve Linear System for gamb
    gamb = copy(linsys.b_bf)
    iad.implicit_linear(linsys.A_bb, gamb; lsolve=ldiv!, Af=linsys.A_bb_LU)

    # - Get body-induced velocities on rotors - #
    vzb = zeros(TF, nbe)
    vrb = zeros(TF, nbe)
    for irotor in 1:length(blade_elements.Omega)
        vzb[:, irotor] = ivr.vz_rb[irotor] * gamb[1:nbodynodes]
        vrb[:, irotor] = ivr.vr_rb[irotor] * gamb[1:nbodynodes]
    end

    ##### ----- Loop through rotors ----- #####
    for irotor in 1:length(blade_elements.Omega)

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
                op.Vinf + vz, # axial velocity V is freestream, vz is induced by bodies and rotor(s) ahead
                blade_elements.Omega[irotor] *
                blade_elements.rotor_panel_centers[:, irotor] + vt, # tangential velocity
                op.rhoinf,
                0.0, #pitch is zero
                op.muinf,
                op.asound,
            ) for (ir, (vz, vt)) in enumerate(zip(vzind, vthetaind))
        ]

        # solve CCBlade problem for this rotor
        out = c4b.solve.(Ref(rotor), sections, c4bop)

        ##### ----- Assign Initial Vz_rotor, Vtheta_rotor, and Cm_wake ----- #####

        # -  Vz_rotor and V_theta rotor - #
        # self influence
        Vz_rotor[:, irotor] .+= vzind .+ out.u
        Vtheta_rotor[:, irotor] .+= vthetaind .+ out.v

        # - Get Cm_wake - #
        # note: we are going to estimate this by taking the velocities on the rotors (though using far field z terms) and applying them constantly straight back to the next rotor or end of wake.
        # Get source strengths
        # note: we need the values at the nodes not centers, so average the values and use the end values on the end points
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
            vrind .+= vr_rr[:, :, jrotor] * sigr[:, jrotor]
        end

        # add in axial and tangential influence aft of current rotor
        vzind .+= 2.0 * out.u
        vthetaind .-= 2.0 * out.v

        # note: since wakes extend from source panel endpoints, we need to average velocities and use the ends for endpoints
        Cm_wake_vec[1] = sqrt((op.Vinf + vzind[1])^2 + vrind[1]^2)
        Cm_wake_vec[2:(end - 1)] =
            (
                sqrt.((op.Vinf .+ vzind[2:end]) .^ 2 .+ vrind[2:end] .^ 2) .+
                sqrt.((op.Vinf .+ vzind[1:(end - 1)]) .^ 2 .+ vrind[1:(end - 1)] .^ 2)
            ) / 2.0
        Cm_wake_vec[end] = sqrt((op.Vinf + vzind[end])^2 + vrind[end]^2)

        # fill in the section of the wake aft of the current rotor and up to the next rotor (or end of wake)
        # TODO: need to change pre-computation to get rotorwakeNODEid rather than rotorwakePANELid. it's only used here as well.
        for (wid, wmap) in enumerate(eachrow(inputs.rotorwakenodeid))
            if wmap[2] >= irotor && wmap[2] < irotor + 1
                Cm_wake[wid] = Cm_wake_vec[wmap[1]]
            end
        end
    end # loop through rotors

    return Vz_rotor, Vtheta_rotor, Cm_wake
end

"""
TODO; move this to the post-processing file
"""
function postcompute_parameters_aid!(cache, propulsor)
    # - Re-compute the precomputed stuff - #
    # - Compute the rest of the stuff required for post-processing - #
    return nothing
end
