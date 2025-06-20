"""
    initialize_velocities(
        solver_options::SolverOptionsType,
        operating_point,
        blade_elements,
        linsys,
        ivr,
        ivw,
        body_totnodes,
        wake_panel_sheet_be_map,
    )

Initialize velocity state variables.

# Arguments
- `solver_options::SolverOptionsType` : solver options type for dispatch
- `operating_point::OperatingPoint` : an OperatingPoint object
- `blade_elements::NamedTuple` : A named tuple containing the blade element geometry and airfoil information.
- `linsys::NamedTuple` : A named tuple containing the panel method linear system information.
- `ivr::NamedTuple` : A named tuple containing the unit induced velocities on the rotors
- `ivw::NamedTuple` : A named tuple containing the unit induced velocities on the wake
- `body_totnodes::Int` : the total number of panel nodes comprising the duct and center_body geometry
- `wake_panel_sheet_be_map::Matrix{Int}` : An index map from the wake panels to the nearest ahead rotor blade element along the wake sheets

# Returns
- `vz_rotor::Vector{Float}` : a vector of the velocity state variables associated with the rotor axially induced velocity
- `vtheta_rotor::Vector{Float}` : a vector of the velocity state variables associated with the rotor tangentially induced velocity
- `Cm_wake::Vector{Float}` : a vector of the velocity state variables associated with the wake control point meridional velocity
"""
function initialize_velocities(
    solver_options::TS,
    operating_point,
    blade_elements,
    linsys,
    ivr,
    ivw,
    body_totnodes,
    wake_panel_sheet_be_map,
) where {TS<:ExternalSolverOptions}

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
        solver_options,
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

"""
function initialize_velocities!(
    solver_options::SolverOptionsType,
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

In-place version of `initialize_velocities`.
"""
function initialize_velocities!(
    solver_options::TS,
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
) where {TS<:ExternalSolverOptions}

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
    rhs = copy(-linsys.b_bf)
    gamb = ImplicitAD.implicit_linear(linsys.A_bb, rhs; lsolve=ldiv!, Af=linsys.A_bb_LU)
    # gamb = ImplicitAD.implicit_linear(
    #     linsys.A_bb, copy(linsys.b_bf); lsolve=ldiv!, Af=linsys.A_bb_LU
    # )
    # gamb = zeros(size(ivr.v_rb, 2) + 2)

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
            is_stator=blade_elements.is_stator[irotor],
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

    return vz_rotor, vtheta_rotor, Cm_wake
end

"""
    initialize_strengths!(
        solver_options::SolverOptionsType,
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
        wake_node_ids_along_center_body_wake_interface,
    )

Initialize strength state variables.

# Arguments
- `solver_options::SolverOptionsType` : solver options type for dispatch
- `Gamr::Vector{Float}` : Rotor circulation state variables (modified in place)
- `sigr::Vector{Float}` : Rotor panel strength state variables (modified in place)
- `gamw::Vector{Float}` : Wake panel strength state variables (modified in place)
- `operating_point::OperatingPoint` : an OperatingPoint object
- `blade_elements::NamedTuple` : A named tuple containing the blade element geometry and airfoil information.
- `linsys::NamedTuple` : A named tuple containing the panel method linear system information.
- `ivr::NamedTuple` : A named tuple containing the unit induced velocities on the rotors
- `ivw::NamedTuple` : A named tuple containing the unit induced velocities on the wake
- `wakeK::Vector{Float}` : geometric constants of wake nodes used in calculating wake strengths
- `body_totnodes::Int` : the total number of panel nodes comprising the duct and center_body geometry
- `wake_nodemap::Matrix{Int}` : an index map of wake panel to the associated node indices
- `wake_endnodeidxs::Matrix{Int}` : the node indices of the start and end of the wake sheets.
- `wake_panel_sheet_be_map::Matrix{Int}` : An index map from the wake panels to the nearest ahead rotor blade element along the wake sheets
- `wake_node_sheet_be_map::Matrix{Int}` : An index map from the wake nodes to the nearest ahead rotor blade element along the wake sheets
- `wake_node_ids_along_casing_wake_interface::type` : An index map indicating which wake nodes interface with the duct wall
- `wake_node_ids_along_center_body_wake_interface::type` : An index map indicating which wake nodes interface with the center_body wall

# Returns
Updates `Gamr`, `sigr`, `gamw`.
"""
function initialize_strengths!(
    solver_options,
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
    wake_node_ids_along_center_body_wake_interface,
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
            is_stator=blade_elements.is_stator[irotor],
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

        if !iszero(blade_elements.is_stator[irotor])
            Gamr[:, irotor] .*= -1.0
        end

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
        wake_node_ids_along_center_body_wake_interface;
    )

    # Gamr struggles to converge if it's not initially positive...
    for g in eachindex(Gamr)
        if Gamr[g] < 0
            Gamr[g] = 0.05
        end
    end

    return Gamr, sigr, gamw
end

"""
    function initialize_strengths!(
        solver_options::CSORSolverOptions,
        Gamr,
        sigr,
        gamw,
        solve_containers,
        operating_point,
        blade_elements,
        wakeK,
        wake_nodemap,
        wake_endnodeidxs,
        wake_panel_sheet_be_map,
        wake_node_sheet_be_map,
        wake_node_ids_along_casing_wake_interface,
        wake_node_ids_along_center_body_wake_interface;
        niter=10,
        rlx=0.5,
    )

Refactored from DFDC SUBROUTINE ROTINITBLD

From the subroutine notes:
Sets reasonable initial circulation using current
rotor blade geometry (chord, beta).
Initial circulations are set w/o induced effects
An iteration is done using the self-induced velocity
from momentum theory to converge an approximate
induced axial velocity

# Arguments
- `solver_options::SolverOptionsType` : solver options type for dispatch
- `Gamr::Vector{Float}` : Rotor circulation state variables (modified in place)
- `sigr::Vector{Float}` : Rotor panel strength state variables (modified in place)
- `gamw::Vector{Float}` : Wake panel strength state variables (modified in place)
- `operating_point::OperatingPoint` : an OperatingPoint object
- `blade_elements::NamedTuple` : A named tuple containing the blade element geometry and airfoil information.
- `wakeK::Vector{Float}` : geometric constants of wake nodes used in calculating wake strengths
- `wake_nodemap::Matrix{Int}` : an index map of wake panel to the associated node indices
- `wake_endnodeidxs::Matrix{Int}` : the node indices of the start and end of the wake sheets.
- `wake_panel_sheet_be_map::Matrix{Int}` : An index map from the wake panels to the nearest ahead rotor blade element along the wake sheets
- `wake_node_sheet_be_map::Matrix{Int}` : An index map from the wake nodes to the nearest ahead rotor blade element along the wake sheets
- `wake_node_ids_along_casing_wake_interface::type` : An index map indicating which wake nodes interface with the duct wall
- `wake_node_ids_along_center_body_wake_interface::type` : An index map indicating which wake nodes interface with the center_body wall

# Keyword Arguments
- `rlx::Float=0.5` : factor for under-relaxation to reduce transients in CL

# Returns
Updated `Gamr`, `sigr`, and `gamw`.
"""
function initialize_strengths!(
    solver_options::CSORSolverOptions,
    Gamr,
    sigr,
    gamw,
    solve_containers,
    operating_point,
    blade_elements,
    wakeK,
    wake_nodemap,
    wake_endnodeidxs,
    wake_panel_sheet_be_map,
    wake_node_sheet_be_map,
    wake_node_ids_along_casing_wake_interface,
    wake_node_ids_along_center_body_wake_interface;
    niter=10,
    rlx=0.5,
    tip_gap=zeros(10), # just make sure it's long enough to not break stuff for now.
)

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
        eltype(blade_elements.twists),
        eltype(blade_elements.chords),
        eltype(blade_elements.Rtip),
        eltype(blade_elements.Rhub),
    )

    # - Rename for Convenience - #
    nbe, nrotor = size(Gamr)
    (; Vinf, rhoinf, muinf, Omega) = operating_point
    (; rotor_panel_centers, B, chords, Rtip, Rhub) = blade_elements

    reset_containers!(solve_containers)

    vz_rotor = -Vinf[] + 1.0

    # loop through rotors
    for irotor in 1:nrotor
        segment_length = (Rtip[irotor] - Rhub[irotor]) / nbe

        # do niter iterations
        for iter in 1:niter
            tsum = 0.0

            # Use upstream circulation to calculate inflow
            if irotor > 1
                @views solve_containers.vtheta_rotor[:, irotor] .=
                    sum(B[1:(irotor - 1)] * Gamr[:, 1:(irotor - 1)]) ./
                    (2.0 * pi * rotor_panel_centers[:, irotor])
            end

            @views @. solve_containers.Cz_rotor[:, irotor] = Vinf[] + vz_rotor

            @views @. solve_containers.Ctheta_rotor[:, irotor] .=
                solve_containers.vtheta_rotor[:, irotor] .-
                rotor_panel_centers[:, irotor] .* Omega[irotor]

            @. solve_containers.Cmag_rotor = sqrt(
                solve_containers.Ctheta_rotor^2 + solve_containers.Cz_rotor^2
            )

            # get cl
            calculate_blade_element_coefficients!(
                solve_containers.cl,
                solve_containers.cd,
                solve_containers.beta1,
                solve_containers.alpha,
                solve_containers.reynolds,
                solve_containers.mach,
                blade_elements,
                solve_containers.Cz_rotor,
                solve_containers.Ctheta_rotor,
                solve_containers.Cmag_rotor,
                operating_point;
            )

            @views BGamr_est =
                0.5 * solve_containers.cl[:, irotor] .*
                solve_containers.Cmag_rotor[:, irotor] .* chords[:, irotor] * B[irotor]

            #TODO: figure out how to handle stators in this case
            # if !iszero(blade_elements.is_stator[irotor])
            #     BGamr_est .*= -1.0
            # end

            delta_BGamr = BGamr_est .- B[irotor] * Gamr[:, irotor]

            Gamr[:, irotor] .+= (rlx * delta_BGamr) / B[irotor]

            @views tsum -= sum(
                B[irotor] .* Gamr[:, irotor] .* rhoinf[] .*
                solve_containers.Ctheta_rotor[:, irotor] .* segment_length[irotor],
            )

            # if tip gap, set dummy circulation to zero
            if tip_gap[irotor] > 0.0
                Gamr[end, irotor] = 0.0
            end

            # use momentum theory estimate of duct induced axial velocity to set VA
            VHSQ = tsum / (rhoinf[] * pi * (Rtip[irotor]^2 - Rhub[irotor]^2))
            if irotor == 1
                vz_rotor = -0.5 * Vinf[] + sqrt((0.5 * Vinf[])^2 + VHSQ)
            end
        end  # for iter

        # - Initialize Wake Strengths - #

        # Set average velocity in duct
        for (wid, wmap) in enumerate(eachrow(wake_panel_sheet_be_map))
            if wmap[2] >= irotor && wmap[2] < irotor + 1
                solve_containers.Cm_wake[wid] = vz_rotor + Vinf[]
            end
        end
    end # for irotor

    # - initialize wake strengths - #
    # TODO: these should be in solve_containers, but need to figure out how to organize that as an input in this case
    Gamma_tilde = zeros(TF, nbe, nrotor)
    H_tilde = zeros(TF, nbe, nrotor)
    deltaGamma2 = zeros(TF, nbe + 1, nrotor)
    deltaH = zeros(TF, nbe + 1, nrotor)
    Cm_avg = zeros(TF, size(gamw)) .= 0

    average_wake_velocities!(
        Cm_avg, solve_containers.Cm_wake, wake_nodemap, wake_endnodeidxs
    )

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
        wake_node_ids_along_center_body_wake_interface;
    )

    return Gamr, gamw, sigr
end
