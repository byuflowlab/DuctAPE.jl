"""
Calculate steady pressure coefficient
"""
function steady_cp(Vs, Vinf, Vref)
    cp = similar(Vs) .= 0
    return steady_cp!(cp, Vs, Vinf, Vref)
end

function steady_cp!(cp, Vs, Vinf, Vref)
    cp .= (Vinf^2 .- Vs .^ 2) / Vref^2
    return cp
end
"""
only used in post-process for cp.
expression not in dfdc theory, comes from source code.
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
Calculate change in pressure coefficient aft of rotor, due to rotor
"""
function delta_cp(deltaH, deltaS, Vtheta, Vref)
    if isapprox(Vref, 0.0)
        return 0.0
    else
        return (2.0 * (deltaH - deltaS) .- Vtheta .^ 2) / Vref^2
    end
end

"""
Calculate net circulation and enthalpy and entropy disk jumps
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
Calculate change in pressure coefficient due to rotors specifically on the body panels aft of the rotors
"""
function calculate_body_delta_cp!(cp, Gamr, sigr, Cz_rotor, Vref, Omega, B, cpr, didr, hidr)

    ## -- Calculate change in pressure coefficient -- ##

    Gamma_tilde, Htilde, Stilde = calculate_rotor_jumps(Gamr, Omega, B, sigr, Cz_rotor)

    nrotor = size(Gamr, 2)

    for irotor in 1:nrotor

        # - Get the tangential velocities on the bodies - #
        v_theta_duct = calculate_vtheta(
            Gamma_tilde[end, irotor], @view(cpr[1:didr[irotor]])
        )
        v_theta_hub = calculate_vtheta(Gamma_tilde[1, irotor], @view(cpr[hidr[irotor]:end]))

        # assemble change in cp due to enthalpy and entropy behind rotor(s)
        cp[1:didr[irotor]] += delta_cp(
            Htilde[end, irotor], Stilde[end, irotor], v_theta_duct, Vref
        )
        cp[hidr[irotor]:end] += delta_cp(
            Htilde[1, irotor], Stilde[1, irotor], v_theta_hub, Vref
        )
    end

    return nothing
end

"""
Calculate change in pressure coefficient due to rotors specifically on the body wakes
"""
function calculate_bodywake_delta_cp(Gamr, sigr, Cz_rotor, Vref, Omega, B, r; body="duct")

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

    v_theta_wake = calculate_vtheta(gt, r)

    # assemble change in cp due to enthalpy and entropy behind rotor(s)
    deltacp = delta_cp(ht, st, v_theta_wake, Vref)

    return deltacp
end

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
    didr,
    hidr,
    controlpoints,
    endpanelidxs,
    zpts,
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
        didr,
        hidr,
        controlpoints,
        endpanelidxs,
        zpts,
    )
end

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
    didr,
    hidr,
    controlpoints,
    endpanelidxs,
    zpts,
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

    # - add the change in Cp on the walls due to enthalpy, entropy, and vtheta - #
    calculate_body_delta_cp!(
        cp_out, Gamr, sigr, Cz_rotor, Vref, Omega, B, @view(controlpoints[2, :]), didr, hidr
    )

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
        endpanelidxs,
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
        endpanelidxs,
    )

    return cp_tuple
end

"""
Calculate the pressure coefficient distributions on one of the body wakes
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
Calculate dimensional and non-dimensional axial force on a single body
"""
function forces_from_pressure(cp_in, cp_out, panels; rhoinf=1.225, Vref=1.0)
    # - initialize - #
    cfx = zeros(eltype(cp_out), Int(panels.nbodies[])) # axial force coefficient (all others are zero for axisymmetric case)
    CFx = similar(cfx) .= 0

    return forces_from_pressure!(CFx, cfx, cp_in, cp_out, panels; rhoinf=rhoinf, Vref=Vref)
end

function forces_from_pressure!(CFx, cfx, cp_in, cp_out, panels; rhoinf=1.225, Vref=1.0)

    # - rename for convenience - #
    #just want x-component of normals since it's axisymmetric
    ns = @view(panels.normal[1, :])
    #radial positions
    rs = @view(panels.controlpoint[2, :])
    #panel lengths
    ds = panels.influence_length

    # for each body
    for ib in 1:(Int(panels.nbodies[]))
        # - rectangular integration due to constant panel strengths. - #
        for ip in Int.(panels.endpanelidxs[1, ib]:panels.endpanelidxs[2, ib])
            cfx[ib] += (cp_out[ip] - cp_in[ip]) * ns[ip] * ds[ip] * 2.0 * pi * rs[ip]
        end
    end

    #dimensionalize
    CFx .= cfx .* 0.5 * rhoinf * Vref^2

    #note, thrust is in negative x-direction
    return CFx, cfx
end

"""
Calculate dimensional and non-dimensional axial force on a single body
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
