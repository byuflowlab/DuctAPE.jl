"""
"""
function initialize_velocities(
    operating_point,
    blade_elements,
    linsys,
    ivr,
    ivw,
    body_totnodes,
    wake_panel_sheet_be_map,
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
        eltype(ivr.v_rr),
        eltype(blade_elements.twists),
        eltype(blade_elements.chords),
        eltype(blade_elements.Rtip),
        eltype(blade_elements.Rhub),
    )

    # - Initialize velocity intermediates and outputs - #
    # rename to only compute once
    nbe = size(blade_elements.rotor_panel_centers)
    nbe = size(blade_elements.rotor_panel_centers, 1)

    # outputs
    vz_rotor = zeros(TF, nbe)
    vtheta_rotor = zeros(TF, nbe)
    Cm_wake = zeros(TF, size(ivw.v_ww, 1))

    return initialize_velocities!(
        vz_rotor,
        vtheta_rotor,
        Cm_wake,
        operating_point,
        blade_elements,
        linsys,
        ivr,
        ivw,
        body_totnodes,
        wake_panel_sheet_be_map,
    )
end

function initialize_velocities!(
    vz_rotor,
    vtheta_rotor,
    Cm_wake,
    operating_point,
    blade_elements,
    linsys,
    ivr,
    ivw,
    body_totnodes,
    wake_panel_sheet_be_map,
)

    # zero outputs:
    vz_rotor .= 0
    vtheta_rotor .= 0
    Cm_wake .= 0

    # - get floating point type - #
    TF = promote_type(
        eltype(operating_point.Omega),
        eltype(operating_point.Vinf),
        eltype(operating_point.rhoinf),
        eltype(operating_point.muinf),
        eltype(operating_point.asound),
        eltype(linsys.A_bw),
        eltype(ivr.v_rr),
        eltype(blade_elements.twists),
        eltype(blade_elements.chords),
        eltype(blade_elements.Rtip),
        eltype(blade_elements.Rhub),
    )

    # - Initialize intermediates - #
    # rename to only compute once
    nbe, nrotor = size(blade_elements.rotor_panel_centers)

    # intermediate values
    # TODO: put these in a precomp container cache eventually
    sigr = zeros(TF, nbe + 1, nrotor)
    Cm_wake_vec = zeros(TF, nbe + 1)
    vzind = zeros(TF, nbe)
    vrind = zeros(TF, nbe)
    vthetaind = zeros(TF, nbe)

    # Solve Linear System for gamb
    # TODO; consider having an option here where you can fill the rhs cache (which should be used here) based on the reference velocity to try and get a better starting point
    # #probably set that up in the precompute parameters function as this would be the first place that rhs vector would be seen.
    gamb = ImplicitAD.implicit_linear(
        linsys.A_bb, copy(linsys.b_bf); lsolve=ldiv!, Af=linsys.A_bb_LU
    )

    # - Get body-induced velocities on rotors - #
    vzb = zeros(TF, nbe, nrotor)
    vrb = zeros(TF, nbe, nrotor)
    for irotor in 1:length(operating_point.Omega)
        berange = (nbe * (irotor - 1) + 1):(nbe * irotor)
        vzb[:, irotor] = ivr.v_rb[berange, :, 1] * gamb[1:body_totnodes]
        vrb[:, irotor] = ivr.v_rb[berange, :, 2] * gamb[1:body_totnodes]
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
                operating_point.Vinf[] + vz, # axial velocity V is freestream, vz is induced by bodies and rotor(s) ahead
                operating_point.Omega[irotor] *
                blade_elements.rotor_panel_centers[ir, irotor] + vt, # tangential velocity
                operating_point.rhoinf[],
                0.0, #pitch is zero
                operating_point.muinf[],
                operating_point.asound[],
            ) for (ir, (vz, vt)) in enumerate(zip(vzind, vthetaind))
        ]

        # solve CCBlade problem for this rotor
        out = c4b.solve.(Ref(rotor), sections, c4bop)

        ##### ----- Assign Initial vz_rotor, vtheta_rotor, and Cm_wake ----- #####

        # -  vz_rotor and V_theta rotor - #
        # self influence
        vz_rotor[:, irotor] .+= vzind .+ getfield.(out, :u)
        vtheta_rotor[:, irotor] .+= vthetaind .+ getfield.(out, :v)

        # - Get Cm_wake - #
        #=
          NOTE: we are going to estimate this by taking the velocities on the rotors (though using far field z terms) and applying them constantly straight back to the next rotor or end of wake.
        =#
        # Get source strengths
        #=
          NOTE: we need the values at the nodes not centers, so average the values and use the end values on the end points
        =#
        sigr[1] = @. blade_elements.B[irotor] / (4.0 * pi) *
            getfield.(out, :cd)[1] *
            getfield.(out, :W)[1] *
            blade_elements.chords[1, irotor]
        @. sigr[2:(end - 1)] =
            (
                blade_elements.B[irotor] / (4.0 * pi) *
                getfield.(out, :cd)[2:end] *
                getfield.(out, :W)[2:end] *
                blade_elements.chords[2:end, irotor] +
                blade_elements.B[irotor] / (4.0 * pi) *
                getfield.(out, :cd)[1:(end - 1)] *
                getfield.(out, :W)[1:(end - 1)] *
                blade_elements.chords[1:(end - 1), irotor]
            ) / 2.0
        sigr[end] = @. blade_elements.B[irotor] / (4.0 * pi) *
            getfield.(out, :cd)[end] *
            getfield.(out, :W)[end] *
            blade_elements.chords[end, irotor]

        # add influence of rotor radial induced velocity from self and rotors ahead
        for jrotor in 1:irotor
            berange = (nbe * (irotor - 1) + 1):(nbe * irotor)
            berangej = ((nbe + 1) * (jrotor - 1) + 1):((nbe + 1) * jrotor)
            vrind .+= ivr.v_rr[berange, berangej, 2] * sigr[:, jrotor]
        end

        # add in axial and tangential influence aft of current rotor
        vzind .+= 2.0 * getfield.(out, :u)
        vthetaind .-= 2.0 * getfield.(out, :v)

        # since wakes extend from source panel endpoints, we need to average velocities and use the ends for endpoints
        Cm_wake_vec[1] = sqrt((operating_point.Vinf[1] + vzind[1])^2 + vrind[1]^2)
        Cm_wake_vec[2:(end - 1)] =
            (
                sqrt.(
                    (operating_point.Vinf[1] .+ vzind[2:end]) .^ 2 .+ vrind[2:end] .^ 2
                ) .+
                sqrt.(
                    (operating_point.Vinf[1] .+ vzind[1:(end - 1)]) .^ 2 .+
                    vrind[1:(end - 1)] .^ 2
                )
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

function initialize_strengths!(
    Gamr,
    sigr,
    gamw,
    operating_point,
    blade_elements,
    linsys,
    ivr,
    ivw,
    wakeK,
    body_totnodes,
    wake_nodemap,
    wake_endnodeidxs,
    wake_panel_sheet_be_map,
    wake_node_sheet_be_map,
    wake_node_ids_along_casing_wake_interface,
    wake_node_ids_along_centerbody_wake_interface,
)

    # zero outputs:
    Gamr .= 0
    sigr .= 0
    gamw .= 0

    # - get floating point type - #
    TF = promote_type(
        eltype(Gamr),
        eltype(sigr),
        eltype(gamw),
        eltype(operating_point.Omega),
        eltype(operating_point.Vinf),
        eltype(operating_point.rhoinf),
        eltype(operating_point.muinf),
        eltype(operating_point.asound),
        eltype(linsys.A_bw),
        eltype(ivr.v_rr),
        eltype(blade_elements.twists),
        eltype(blade_elements.chords),
        eltype(blade_elements.Rtip),
        eltype(blade_elements.Rhub),
    )

    # - Initialize intermediates - #
    # rename to only compute once
    nbe, nrotor = size(blade_elements.rotor_panel_centers)

    # intermediate values
    # TODO: put these in a precomp container cache eventually
    vz_rotor = zeros(TF, nbe, nrotor)
    vtheta_rotor = zeros(TF, nbe, nrotor)
    Cm_wake_vec = zeros(TF, nbe + 1)
    Cm_wake = zeros(TF, size(wake_panel_sheet_be_map, 1)) .= 0
    vthetaind = zeros(TF, nbe)
    vzind = zeros(TF, nbe)
    vrind = zeros(TF, nbe)

    # Solve Linear System for gamb
    # gamb = ImplicitAD.implicit_linear(
    #     linsys.A_bb, copy(linsys.b_bf); lsolve=ldiv!, Af=linsys.A_bb_LU
    # )
    gamb = zeros(size(ivr.v_rb, 2) + 2)

    # - Get body-induced velocities on rotors - #
    vzb = zeros(TF, nbe, nrotor)
    vrb = zeros(TF, nbe, nrotor)
    for irotor in 1:length(operating_point.Omega)
        berange = (nbe * (irotor - 1) + 1):(nbe * irotor)
        vzb[:, irotor] = ivr.v_rb[berange, :, 1] * gamb[1:body_totnodes]
        vrb[:, irotor] = ivr.v_rb[berange, :, 2] * gamb[1:body_totnodes]
    end

    ##### ----- Loop through rotors ----- #####
    for irotor in 1:nrotor

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
                operating_point.Vinf[] + vz, # axial velocity V is freestream, vz is induced by bodies and rotor(s) ahead
                operating_point.Omega[irotor] *
                blade_elements.rotor_panel_centers[ir, irotor] .- vt, # tangential velocity
                operating_point.rhoinf[],
                0.0, #pitch is zero
                operating_point.muinf[],
                operating_point.asound[],
            ) for (ir, (vz, vt)) in enumerate(zip(vzind, vthetaind))
        ]

        # solve CCBlade problem for this rotor
        out = c4b.solve.(Ref(rotor), sections, c4bop)

        ##### ----- Assign Initial Gamr, sigr, and gamw ----- #####
        # - Get Gamr - #
        Gamr[:, irotor] .=
            0.5 .* getfield.(out, :cl) .* getfield.(out, :W) .*
            blade_elements.chords[:, irotor]

        # # -  vz_rotor and V_theta rotor - #
        # # self influence
        # vz_rotor[:, irotor] .+= vzind .+ getfield.(out,:u
        # vtheta_rotor[:, irotor] .+= vthetaind .+ getfield.(out,:v

        # - Get Cm_wake - #
        #=
          NOTE: we are going to estimate this by taking the velocities on the rotors (though using far field z terms) and applying them constantly straight back to the next rotor or end of wake.
        =#
        # Get source strengths
        #=
          NOTE: we need the values at the nodes not centers, so average the values and use the end values on the end points
        =#
        sigr[1, irotor] = @. blade_elements.B[irotor] / (4.0 * pi) *
            getfield.(out, :cd)[1] *
            getfield.(out, :W)[1] *
            blade_elements.chords[1, irotor]
        @. sigr[2:(end - 1), irotor] =
            (
                blade_elements.B[irotor] / (4.0 * pi) *
                getfield.(out, :cd)[2:end] *
                getfield.(out, :W)[2:end] *
                blade_elements.chords[2:end, irotor] +
                blade_elements.B[irotor] / (4.0 * pi) *
                getfield.(out, :cd)[1:(end - 1)] *
                getfield.(out, :W)[1:(end - 1)] *
                blade_elements.chords[1:(end - 1), irotor]
            ) / 2.0
        sigr[end, irotor] = @. blade_elements.B[irotor] / (4.0 * pi) *
            getfield.(out, :cd)[end] *
            getfield.(out, :W)[end] *
            blade_elements.chords[end, irotor]

        # add influence of rotor radial induced velocity from self and rotors ahead
        for jrotor in 1:irotor
            berange = (nbe * (irotor - 1) + 1):(nbe * irotor)
            berangej = ((nbe + 1) * (jrotor - 1) + 1):((nbe + 1) * jrotor)
            vrind .+= ivr.v_rr[berange, berangej, 2] * sigr[:, jrotor]
        end

        # add in axial and tangential influence aft of current rotor
        vzind .+= 2.0 * getfield.(out, :u)
        vthetaind .-= 2.0 * getfield.(out, :v)

        # since wakes extend from source panel endpoints, we need to average velocities and use the ends for endpoints
        Cm_wake_vec[1] = sqrt((operating_point.Vinf[1] + vzind[1])^2 + vrind[1]^2)
        Cm_wake_vec[2:(end - 1)] =
            (
                sqrt.(
                    (operating_point.Vinf[1] .+ vzind[2:end]) .^ 2 .+ vrind[2:end] .^ 2
                ) .+
                sqrt.(
                    (operating_point.Vinf[1] .+ vzind[1:(end - 1)]) .^ 2 .+
                    vrind[1:(end - 1)] .^ 2
                )
            ) / 2.0
        Cm_wake_vec[end] = sqrt((operating_point.Vinf[1] + vzind[end])^2 + vrind[end]^2)

        # fill in the section of the wake aft of the current rotor and up to the next rotor (or end of wake)
        for (wid, wmap) in enumerate(eachrow(wake_panel_sheet_be_map))
            if wmap[2] >= irotor && wmap[2] < irotor + 1
                Cm_wake[wid] = Cm_wake_vec[wmap[1]]
            end
        end
    end # loop through rotors

    # - initialize wake strengths - #
    # TODO: these should be in solve_containers, but need to figure out how to organize that as an input in this case
    Gamma_tilde = zeros(TF, nbe, nrotor)
    H_tilde = zeros(TF, nbe, nrotor)
    deltaGamma2 = zeros(TF, nbe + 1, nrotor)
    deltaH = zeros(TF, nbe + 1, nrotor)
    Cm_avg = zeros(TF, size(gamw)) .= 0

    average_wake_velocities!(Cm_avg, Cm_wake, wake_nodemap, wake_endnodeidxs)

    # - Calculate Wake Panel Strengths - #
    # in-place solve for gamw,
    calculate_wake_vortex_strengths!(
        gamw,
        Gamma_tilde,
        H_tilde,
        deltaGamma2,
        deltaH,
        Gamr,
        Cm_avg,
        blade_elements.B,
        operating_point.Omega,
        wakeK,
        wake_node_sheet_be_map,
        wake_node_ids_along_casing_wake_interface,
        wake_node_ids_along_centerbody_wake_interface;
    )

    # Gamr struggles to converge if it's not initially positive...
    for g in eachindex(Gamr)
        if Gamr[g] < 0
            Gamr[g] = 0.05
        end
    end

    return Gamr, sigr, gamw
end
