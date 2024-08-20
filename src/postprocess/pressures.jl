"""
    steady_cp(Vs, Vinf, Vref)

Calculate steady pressure coefficients for a given surface velocity.

# Arguments
- `Vs::Vector{Float}` : the surface velocities
- `Vinf::Vector{Float}` : one element vector with freestream mangnitude
- `Vref::Vector{Float}` : one element vector with reference velocity used for non-dimensionalization

# Returns
- `cp::Vector{Float}` : the steady pressure coefficients
"""
function steady_cp(Vs, Vinf, Vref)
    cp = similar(Vs) .= 0
    return steady_cp!(cp, Vs, Vinf, Vref)
end

"""
    steady_cp!(cp, Vs, Vinf, Vref)

In-place verison of `steady_cp`.
"""
function steady_cp!(cp, Vs, Vinf, Vref)
    cp .= (Vinf^2 .- Vs .^ 2) / Vref^2
    return cp
end

"""
    calculate_entropy_jumps(sigr, Cz_rotor)

Calculate jumps in entropy across the disks.

# Arguments
- `sigr::Matrix{Float}` : rotor source panel strengths
- `Cz_rotor::Vector{Float}` : absolute axial velocity on rotor blade elements

# Returns
- `deltaS::Vector{Float}` : entropy jump across rotor disks
"""
function calculate_entropy_jumps(sigr, Cz_rotor)
    # average sigr's
    sigr_avg = similar(Cz_rotor, size(sigr, 1) - 1, size(sigr, 2)) .= 0
    for (i, s) in enumerate(eachcol(sigr))
        sigr_avg[:, i] = (s[2:end] + s[1:(end - 1)]) / 2.0
    end

    # multiply by Cz_rotor's
    return sigr_avg .* Cz_rotor
end

"""
    delta_cp(deltaH, deltaS, Ctheta, Vref)

Calculate change in pressure coefficient aft of rotor, due to rotor

# Arguments
- `deltaH::Vector{Float}` : Enthalpy jumps across disks
- `deltaS::Vector{Float}` : Entropy jumps across disks`
- `Ctheta::Vector{Float}` : tangenetial velocity
- `Vref::Vector{Float}` : reference velocity for non-dimensionalization

# Returns
- `delta_cp::Vector{Float}` : pressure rises due to rotor disks
"""
function delta_cp(deltaH, deltaS, Ctheta, Vref)
    if isapprox(Vref, 0.0)
        return 0.0
    else
        return (2.0 * (deltaH - deltaS) .- Ctheta .^ 2) / Vref^2
    end
end

"""
    calculate_rotor_jumps(Gamr, Omega, B, sigr, Cz_rotor)

Calculate net circulation and enthalpy and entropy disk jumps

# Arguments
- `Gamr::Matrix{Float}` : Blade element circulation strengths
- `Omega::Vector{Float}` : rotor rotation rates
- `B::Vector{Float}` : blade count for each rotor (usually integers but could be a float)
- `sigr::Matrix{Float}` : rotor source panel strengths
- `Cz_rotor::Vector{Float}` : absolute axial velocity on rotor blade elements

# Returns
- `Gamma_tilde::Matrix{Float}` : net upstream circulation
- `Htilde::Matrix{Float}` : net upstream enthalpy jumps
- `Stilde::Matrix{Float}` : net upstream entropy jumps
"""
function calculate_rotor_jumps(Gamr, Omega, B, sigr, Cz_rotor)

    # - Calculate net circulations - #
    Gamma_tilde = calculate_net_circulation(Gamr, B)

    # - Calculate enthalpy disk jump - #
    Htilde = calculate_enthalpy_jumps(Gamr, Omega, B)

    # - Calculate entropy disk jump - #
    Stilde = calculate_entropy_jumps(sigr, Cz_rotor)

    return Gamma_tilde, Htilde, Stilde
end

"""
    calculate_body_delta_cp!(cp, Gamr, sigr, Cz_rotor, Vref, Omega, B, cpr, casing_panel_ids_aft_of_rotors, centerbody_panel_ids_aft_of_rotors)

Augment surface pressure by change in pressure coefficient due to rotors specifically on the body panels aft of the rotors.

# Arguments
- `cp::Vector{Float}` : steady pressure coeffients, modified in-place to include rotor effects.
- `Gamr::Matrix{Float}` : Blade element circulation strengths
- `sigr::Matrix{Float}` : rotor source panel strengths
- `Cz_rotor::Vector{Float}` : absolute axial velocity on rotor blade elements
- `Vref::Vector{Float}` : one element vector with reference velocity used for non-dimensionalization
- `Omega::Vector{Float}` : rotor rotation rates
- `B::Vector{Float}` : blade count for each rotor (usually integers but could be a float)
- `cpr::Vector{Float}` : control point radial positions of body panels
- `casing_panel_ids_aft_of_rotors::Vector{Int}` : duct indices of control point radial positions aft of rotors
- `centerbody_panel_ids_aft_of_rotors::Vector{Int}` : centerbody indices of control point radial positions aft of rotors
"""
function calculate_body_delta_cp!(
    cp,
    Gamr,
    sigr,
    Cz_rotor,
    Vref,
    Omega,
    B,
    cpr,
    casing_panel_ids_aft_of_rotors,
    centerbody_panel_ids_aft_of_rotors,
)

    ## -- Calculate change in pressure coefficient -- ##

    Gamma_tilde, Htilde, Stilde = calculate_rotor_jumps(Gamr, Omega, B, sigr, Cz_rotor)

    nrotor = size(Gamr, 2)

    for irotor in 1:nrotor

        # - Get the tangential velocities on the bodies - #
        v_theta_duct = calculate_vtheta(
            Gamma_tilde[end, irotor], @view(cpr[1:casing_panel_ids_aft_of_rotors[irotor]])
        )
        v_theta_hub = calculate_vtheta(
            Gamma_tilde[1, irotor],
            @view(cpr[centerbody_panel_ids_aft_of_rotors[irotor]:end])
        )

        # assemble change in cp due to enthalpy and entropy behind rotor(s)
        cp[1:casing_panel_ids_aft_of_rotors[irotor]] += delta_cp(
            Htilde[end, irotor], Stilde[end, irotor], v_theta_duct, Vref
        )
        cp[centerbody_panel_ids_aft_of_rotors[irotor]:end] += delta_cp(
            Htilde[1, irotor], Stilde[1, irotor], v_theta_hub, Vref
        )
    end

    return nothing
end

"""
    calculate_bodywake_delta_cp(Gamr, sigr, Cz_rotor, Vref, Omega, B, cpr; body="duct")

Calculate change in pressure coefficient due to rotors specifically on the body wakes

# Arguments
- `Gamr::Matrix{Float}` : Blade element circulation strengths
- `sigr::Matrix{Float}` : rotor source panel strengths
- `Cz_rotor::Vector{Float}` : absolute axial velocity on rotor blade elements
- `Vref::Vector{Float}` : one element vector with reference velocity used for non-dimensionalization
- `Omega::Vector{Float}` : rotor rotation rates
- `B::Vector{Float}` : blade count for each rotor (usually integers but could be a float)
- `cpr::Vector{Float}` : control point radial positions of body wake "panels"

# Keyword Arguments
- `body::String="duct"` : flag as to whether the body in question is a duct or centerbody.
"""
function calculate_bodywake_delta_cp(Gamr, sigr, Cz_rotor, Vref, Omega, B, cpr; body="duct")

    ## -- Calculate change in pressure coefficient -- ##

    Gamma_tilde, Htilde, Stilde = calculate_rotor_jumps(Gamr, Omega, B, sigr, Cz_rotor)

    # - Get the tangential velocities on the bodies - #
    if body == "duct"
        gt = Gamma_tilde[end]
        ht = Htilde[end]
        st = Stilde[end]
    else
        gt = Gamma_tilde[1, end]
        ht = Htilde[1, end]
        st = Stilde[1, end]
    end

    v_theta_wake = calculate_vtheta(gt, cpr)

    # assemble change in cp due to enthalpy and entropy behind rotor(s)
    deltacp = delta_cp(ht, st, v_theta_wake, Vref)

    return deltacp
end

"""
get_body_cps(
    Vtan_in,
    Vtan_out,
    Gamr,
    sigr,
    Cz_rotor,
    Vinf,
    Vref,
    B,
    Omega,
    casing_panel_ids_aft_of_rotors,
    centerbody_panel_ids_aft_of_rotors,
    controlpoints,
    endpanelidxs,
    zpts,
)

Description

# Arguments
- `Vtan_in::Vector{Float}` : Tangential velocity on the inside of the body panels
- `Vtan_out::Vector{Float}` : Tangential velocity on the outside of the body panels
- `Gamr::Matrix{Float}` : Blade element circulation strengths
- `sigr::Matrix{Float}` : rotor source panel strengths
- `Cz_rotor::Vector{Float}` : absolute axial velocity on rotor blade elements
- `Vinf::Vector{Float}` : one element vector with freestream mangnitude
- `Vref::Vector{Float}` : one element vector with reference velocity used for non-dimensionalization
- `B::Vector{Float}` : blade count for each rotor (usually integers but could be a float)
- `Omega::Vector{Float}` : rotor rotation rates
- `casing_panel_ids_aft_of_rotors::Vector{Int}` : duct indices of control point radial positions aft of rotors
- `centerbody_panel_ids_aft_of_rotors::Vector{Int}` : centerbody indices of control point radial positions aft of rotors
- `controlpoints::Matrix{Float}` : control point locations for each panel
- `endpanelidxs::Matrix{Int}` : the indices of the first and last panels for each body
- `zpts::NamedTuple` : a named tuple containing the z-coordinates of the control points of the duct casing, duct nacelle, and centerbody.

# Returns
- `cp_tuple::NamedTuple` : body surface velocities and various useful breakdowns thereof.
"""
function get_body_cps(
    Vtan_in,
    Vtan_out,
    Gamr,
    sigr,
    Cz_rotor,
    Vinf,
    Vref,
    B,
    Omega,
    casing_panel_ids_aft_of_rotors,
    centerbody_panel_ids_aft_of_rotors,
    controlpoints,
    endpanelidxs,
    zpts;
    isolated_body=false
)
    cp_in = similar(Vtan_in) .= 0
    cp_out = similar(Vtan_in) .= 0
    cp_casing_in = similar(Vtan_in, size(zpts.casing_zpts)) .= 0
    cp_casing_out = similar(Vtan_in, size(zpts.casing_zpts)) .= 0
    cp_nacelle_in = similar(Vtan_in, size(zpts.nacelle_zpts)) .= 0
    cp_nacelle_out = similar(Vtan_in, size(zpts.nacelle_zpts)) .= 0
    cp_centerbody_in = similar(Vtan_in, size(zpts.centerbody_zpts)) .= 0
    cp_centerbody_out = similar(Vtan_in, size(zpts.centerbody_zpts)) .= 0

    cp_tuple = (;
        cp_in,
        cp_out,
        cp_casing_in,
        cp_casing_out,
        cp_nacelle_in,
        cp_nacelle_out,
        cp_centerbody_in,
        cp_centerbody_out,
    )

    return get_body_cps!(
        cp_tuple,
        Vtan_in,
        Vtan_out,
        Gamr,
        sigr,
        Cz_rotor,
        Vinf,
        Vref,
        B,
        Omega,
        casing_panel_ids_aft_of_rotors,
        centerbody_panel_ids_aft_of_rotors,
        controlpoints,
        endpanelidxs,
        zpts;
        isolated_body=isolated_body
    )
end

"""
    get_body_cps!(
        cp_tuple,
        Vtan_in,
        Vtan_out,
        Gamr,
        sigr,
        Cz_rotor,
        Vinf,
        Vref,
        B,
        Omega,
        duct_panel_ids_aft_of_rotors,
        centerbody_panel_ids_aft_of_rotors,
        controlpoints,
        endpanelidxs,
        zpts,
    )

In-place version of `get_body_cps`.
"""
function get_body_cps!(
    cp_tuple,
    Vtan_in,
    Vtan_out,
    Gamr,
    sigr,
    Cz_rotor,
    Vinf,
    Vref,
    B,
    Omega,
    casing_panel_ids_aft_of_rotors,
    centerbody_panel_ids_aft_of_rotors,
    controlpoints,
    endpanelidxs,
    zpts;
    isolated_body=false
)

    # rename for convenience
    (;
        cp_in,             # surface pressure along inside of bodies
        cp_out,            # surface pressure along outside of bodies
        cp_casing_in,      # surface pressure along inside of casing
        cp_nacelle_in,     # surface pressure along inside of nacell
        cp_centerbody_in,  # surface pressure along inside of centerbody
        cp_casing_out,     # surface pressure along outside of casing
        cp_nacelle_out,    # surface pressure along outside of nacelle
        cp_centerbody_out, # surface pressure along outside of centerbody
    ) = cp_tuple

    (; casing_zpts, nacelle_zpts, centerbody_zpts) = zpts

    # - Calculate standard pressure coefficient expression - #
    steady_cp!(cp_in, Vtan_in, Vinf, Vref)
    steady_cp!(cp_out, Vtan_out, Vinf, Vref)

    if !isolated_body
    # - add the change in Cp on the walls due to enthalpy, entropy, and vtheta - #
    calculate_body_delta_cp!(
        cp_out,
        Gamr,
        sigr,
        Cz_rotor,
        Vref,
        Omega,
        B,
        @view(controlpoints[2, :]),
        casing_panel_ids_aft_of_rotors,
        centerbody_panel_ids_aft_of_rotors,
    )
end

    # - Split body strengths into inner/outer duct and hub - #
    split_bodies!(
        cp_casing_in,
        cp_nacelle_in,
        cp_centerbody_in,
        casing_zpts,
        nacelle_zpts,
        centerbody_zpts,
        cp_in,
        controlpoints,
    )
    split_bodies!(
        cp_casing_out,
        cp_nacelle_out,
        cp_centerbody_out,
        casing_zpts,
        nacelle_zpts,
        centerbody_zpts,
        cp_out,
        controlpoints,
    )

    return cp_tuple
end

"""
    get_bodywake_cps(
        Gamr,
        vz_w,
        vr_w,
        gamw,
        vz_r,
        vr_r,
        sigr,
        vz_b,
        vr_b,
        gamb,
        panels,
        Cz_rotor,
        Omega,
        B,
        Vinf,
        Vref;
        body="duct",
    )

Calculate the pressure coefficient distributions on one of the body wakes

# Arguments
- `Gamr::Matrix{Float}` : Blade element circulation strengths
- `vz_w::Matrix{Float}` : unit axial induced velocity of the wake onto the body wake
- `vr_w::Matrix{Float}` : unit radial induced velocity of the wake onto the body wake
- `gamw::Vector{Float}` : wake panel strengths
- `vz_r::Matrix{Float}` : unit axial induced velocity of the rotor onto the body wake
- `vr_r::Matrix{Float}` : unit radial induced velocity of the rotor onto the body wake
- `sigr::Vector{Float}` : rotor panel strengths
- `vz_b::Matrix{Float}` : unit axial induced velocity of the bodies onto the body wake
- `vr_b::Matrix{Float}` : unit radial induced velocity of the bodies onto the body wake
- `gamb::Vector{Float}` : body panel strengths
- `panels::NamedTuple` : A named tuple containing bodywake "panel" geometries
- `Cz_rotor::Vector{Float}` : absolute axial velocity on rotor blade elements
- `Omega::Vector{Float}` : rotor rotation rates
- `B::Vector{Float}` : blade count for each rotor (usually integers but could be a float)
- `Vinf::Vector{Float}` : one element vector containing the velocity magnitude
- `Vref::Vector{Float}` : one element vector with reference velocity used for non-dimensionalization

# Keyword Arguments
- `body::String="duct"` : flag as to whether the body in question is a duct or centerbody.
"""
function get_bodywake_cps(
    Gamr,
    vz_w,
    vr_w,
    gamw,
    vz_r,
    vr_r,
    sigr,
    vz_b,
    vr_b,
    gamb,
    panels,
    Cz_rotor,
    Omega,
    B,
    Vinf,
    Vref;
    body="duct",
)

    # - Get "surface" velocities - #

    # get induced velocities
    vz_bodywake, vr_bodywake = calculate_induced_velocities_on_bodywake(
        vz_w, vr_w, gamw, vz_r, vr_r, sigr, vz_b, vr_b, gamb, Vinf
    )

    # get "surface" velocities
    Vmat = [vz_bodywake vr_bodywake]
    vtan = [dot(v, t) for (v, t) in zip(eachrow(Vmat), panels.tangent)]

    # - Get steady pressure coefficients - #
    cp_steady = steady_cp(vtan, Vinf, Vref)

    # - Get delta cp - #
    deltacp = calculate_bodywake_delta_cp(
        Gamr, sigr, Cz_rotor, Vref, Omega, B, panels.controlpoint[2, :]; body=body
    )

    return cp_steady .+ deltacp, vtan
end

"""
    forces_from_pressure(cp_in, cp_out, panels; rhoinf=1.225, Vref=1.0)

Calculate dimensional and non-dimensional axial force on a single body

# Arguments
- `cp_in::Vector{Float}` : pressure coefficient on inside of body surfaces
- `cp_out::Vector{Float}` : pressure coefficients on outside of body surfaces
- `panels::NamedTuple` : A named tuple containing panel geometry information

# Keyword Arguments
- `rhoinf::Float=1.225` : reference density for non-dimensionalization
- `Vref::Float=1.0` : reference velocity for non-dimensionalization

# Returns
- `thrust::Vector{Float}` : dimensional axial force
- `force_coeff::Vector{Float}` : non-dimensional axial force
"""
function forces_from_pressure(cp_in, cp_out, panels; rhoinf=1.225, Vref=1.0)
    # - initialize - #
    cfx = zeros(eltype(cp_out), Int(panels.nbodies[])) # axial force coefficient (all others are zero for axisymmetric case)
    thrust = similar(cfx) .= 0

    return forces_from_pressure!(thrust, cfx, cp_in, cp_out, panels; rhoinf=rhoinf, Vref=Vref)
end

"""
    forces_from_pressure!(CFx, cfx, cp_in, cp_out, panels; rhoinf=1.225, Vref=1.0)

In-place version of `forces_from_pressure`.
"""
function forces_from_pressure!(thrust, cfx, cp_in, cp_out, panels; rhoinf=1.225, Vref=1.0)

    # - rename for convenience - #
    #just want x-component of normals since it's axisymmetric
    ns = @view(panels.normal[1, :])
    #radial positions
    rs = @view(panels.controlpoint[2, :])
    #panel lengths
    ds = panels.influence_length

    # for each body
    for ib in 1:(Int(panels.nbodies[]))
        # - rectangular integration due to constant panel pressures. - #
        for ip in Int.(panels.endpanelidxs[1, ib]:panels.endpanelidxs[2, ib])
            cfx[ib] += (cp_out[ip] - cp_in[ip]) * ns[ip] * ds[ip] * 2.0 * pi * rs[ip]
        end
    end

    #dimensionalize
    thrust .= cfx .* 0.5 * rhoinf * Vref^2

    #note, thrust is in negative x-direction
    return thrust, cfx
end

"""
    forces_from_TEpanels!(
        thrust, force_coeff, cp_in, cp_out, panels; rhoinf=1.225, Vref=1.0
    )

Add force induced by trailing edge gap panels to total forces.

# Arguments
- `thrust::Vector{Float}` : dimensional axial force
- `force_coeff::Vector{Float}` : non-dimensional axial force
- `cp_in::Vector{Float}` : pressure coefficient on inside of body surfaces
- `cp_out::Vector{Float}` : pressure coefficients on outside of body surfaces
- `panels::NamedTuple` : A named tuple containing panel geometry information

# Keyword Arguments
- `rhoinf::Float=1.225` : reference density for non-dimensionalization
- `Vref::Float=1.0` : reference velocity for non-dimensionalization
"""
function forces_from_TEpanels!(
    thrust, force_coeff, cp_in, cp_out, panels; rhoinf=1.225, Vref=1.0
)

    #dimensionalize
    q = 0.5 * rhoinf * Vref^2

    for i in 1:(Int(panels.nbodies[]))
        if panels.tenode[i, 2, 2] <= eps()
            # if it's the hub, don't average the first and last, but rather just the last
            cpi = cp_in[Int(panels.endpanelidxs[i, 2])]
            cpo = cp_out[Int(panels.endpanelidxs[i, 2])]
        else
            # if it's the duct, then average the first and last panel
            cpi =
                0.5 * (
                    cp_in[Int(panels.endpanelidxs[1, i])] +
                    cp_in[Int(panels.endpanelidxs[2, i])]
                )
            cpo =
                0.5 * (
                    cp_out[Int(panels.endpanelidxs[1, i])] +
                    cp_out[Int(panels.endpanelidxs[2, i])]
                )
        end

        r = 0.5 * sum(panels.tenode[i, :, 2])

        force_coeff[i] +=
            (cpo - cpi) *
            panels.tenormal[1, i] *
            panels.teinfluence_length[i] *
            2.0 *
            pi *
            r

        thrust[i] +=
            q *
            (cpo - cpi) *
            panels.tenormal[1, i] *
            panels.teinfluence_length[i] *
            2.0 *
            pi *
            r
    end

    return thrust, force_coeff
end
