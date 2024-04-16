function inviscid_rotor_thrust(Ctheta_rotor, Gamma_tilde, rotor_panel_length, rhoinf)
    # initialize
    dTi = similar(Gamma_tilde) .= 0.0
    Tinv = zeros(eltype(Gamma_tilde), size(Gamma_tilde, 2))

    return inviscid_rotor_thrust!(
        Tinv, dTi, Ctheta_rotor, Gamma_tilde, rotor_panel_length, rhoinf
    )
end

function inviscid_rotor_thrust!(
    Tinv, dTi, Ctheta_rotor, Gamma_tilde, rotor_panel_length, rhoinf
)

    # problem dimensions
    nr, nrotor = size(dTi)

    for irotor in 1:nrotor
        for ir in 1:nr
            # section thrust
            dTi[ir, irotor] =
                -rhoinf *
                Gamma_tilde[ir, irotor] *
                Ctheta_rotor[ir, irotor] *
                rotor_panel_length[ir, irotor]
        end
    end

    #sum the section thrust
    Tinv .= sum(dTi; dims=1)

    return Tinv, dTi
end

function viscous_rotor_thrust(
    Cz_rotor, Cmag_rotor, B, chord, rotor_panel_length, cd, rhoinf
)
    #initialize
    dTv = similar(Cz_rotor) .= 0.0
    Tvisc = zeros(eltype(Cz_rotor), size(Cz_rotor, 2))

    return viscous_rotor_thrust!(
        Tvisc, dTv, Cz_rotor, Cmag_rotor, B, chord, rotor_panel_length, cd, rhoinf
    )
end

function viscous_rotor_thrust!(
    Tvisc, dTv, Cz_rotor, Cmag_rotor, B, chord, rotor_panel_length, cd, rhoinf
)

    # get dimensions
    nr, nrotor = size(dTv)

    for irotor in 1:nrotor
        for ir in 1:nr
            # hrwc = 0.5 * rhoinf * Cmag_rotor[ir, irotor] * chord[ir, irotor]
            # bdr = B[irotor] * rotor_panel_length[ir, irotor]
            dTv[ir, irotor] =
                -(0.5 * rhoinf * Cmag_rotor[ir, irotor] * chord[ir, irotor]) *
                cd[ir, irotor] *
                Cz_rotor[ir, irotor] *
                (B[irotor] * rotor_panel_length[ir, irotor])
        end
    end

    Tvisc .= sum(dTv; dims=1)

    return Tvisc, dTv
end

function inviscid_rotor_torque(
    Cz_rotor, rotor_panel_center, rotor_panel_length, Gamma_tilde, rhoinf
)
    # initialize
    dQi = similar(Gamma_tilde) .= 0.0
    Qinv = zeros(eltype(Gamma_tilde), size(Gamma_tilde, 2))

    return inviscid_rotor_torque!(
        Qinv, dQi, Cz_rotor, rotor_panel_center, rotor_panel_length, Gamma_tilde, rhoinf
    )
end

function inviscid_rotor_torque!(
    Qinv, dQi, Cz_rotor, rotor_panel_center, rotor_panel_length, Gamma_tilde, rhoinf
)

    # dimensions
    nr, nrotor = size(dQi)

    for irotor in 1:nrotor
        for ir in 1:nr
            # rdr = rotor_panel_center[ir, irotor] * rotor_panel_length[ir, irotor]
            dQi[ir, irotor] =
                rhoinf *
                Gamma_tilde[ir, irotor] *
                Cz_rotor[ir, irotor] *
                (rotor_panel_center[ir, irotor] * rotor_panel_length[ir, irotor])
        end
    end

    Qinv .= sum(dQi; dims=1)

    return Qinv, dQi
end

function viscous_rotor_torque(
    Wtheta_rotor,
    Cmag_rotor,
    B,
    chord,
    rotor_panel_center,
    rotor_panel_length,
    cd,
    rhoinf;
    TF=nothing,
)
    if isnothing(TF)
        TF = promote_type(
            eltype(Wtheta_rotor),
            eltype(Cmag_rotor),
            eltype(chord),
            eltype(rotor_panel_center),
            eltype(cd),
        )
    end

    dQv = zeros(TF, size(Cmag_rotor))
    Qvisc = zeros(TF, size(B))

    return viscous_rotor_torque!(
        Qvisc,
        dQv,
        Wtheta_rotor,
        Cmag_rotor,
        B,
        chord,
        rotor_panel_center,
        rotor_panel_length,
        cd,
        rhoinf;
        TF=nothing,
    )
end

function viscous_rotor_torque!(
    Qvisc,
    dQv,
    Wtheta_rotor,
    Cmag_rotor,
    B,
    chord,
    rotor_panel_center,
    rotor_panel_length,
    cd,
    rhoinf;
    TF=nothing,
)

    # dimensions
    nr, nrotor = size(dQv)

    # initialize

    for irotor in 1:nrotor
        for ir in 1:nr
            # hrwc = 0.5 * rhoinf * Cmag_rotor[ir, irotor] * chord[ir, irotor]
            # brdr =
            # B[irotor] * rotor_panel_center[ir, irotor] * rotor_panel_length[ir, irotor]
            dQv[ir, irotor] =
                -(0.5 * rhoinf * Cmag_rotor[ir, irotor] * chord[ir, irotor]) *
                cd[ir, irotor] *
                Wtheta_rotor[ir, irotor] *
                (
                    B[irotor] *
                    rotor_panel_center[ir, irotor] *
                    rotor_panel_length[ir, irotor]
                )
        end
    end

    Qvisc .= sum(dQv; dims=1)

    return Qvisc, dQv
end

function rotor_power(Q, dQ, Omega)
    dP = similar(dQ) .= 0.0
    P = similar(Q) .= 0.0

    return rotor_power!(P, dP, Q, dQ, Omega)
end

function rotor_power!(P, dP, Q, dQ, Omega)
    nr, nrotor = size(dP)

    for irotor in 1:nrotor
        for ir in 1:nr
            dP[ir, irotor] = dQ[ir, irotor] * Omega[irotor]
        end

        P[irotor] = Q[irotor] .* Omega[irotor]
    end

    return P, dP
end

function get_total_efficiency(total_thrust, total_power, Vinf)
    TF = promote_type(eltype(total_thrust), eltype(total_power), eltype(Vinf))

    eta = zeros(TF, length(total_thrust))

    return get_total_efficiency!(eta, total_thrust, total_power, Vinf)
end

function get_total_efficiency!(eta, total_thrust, total_power, Vinf)
    for i in 1:length(total_thrust)
        if Vinf <= 0.0 || total_power[i] < eps() || total_thrust[i] <= 0.0
            #do nothing, efficiency can't physically be negative or infinite.
        else
            eta[i] = total_thrust[i] * Vinf / total_power[i]
        end
    end

    return eta
end

function get_induced_efficiency(Tinv, Tduct, Pinv, Vinf)
    TF = promote_type(eltype(Tinv), eltype(Pinv), eltpye(Tduct))
    eta_inv = zeros(TF, size(Tinv))
    return get_induced_efficiency!(eta_inv, Tinv, Tduct, Pinv, Vinf)
end

function get_induced_efficiency!(eta_inv, Tinv, Tduct, Pinv, Vinf)
    for (e, ti, p) in zip(eachrow(eta_inv), Tinv, Pinv)
        if Vinf <= 0.0 || p <= 0.0
            e[1] = 0.0
        else
            e[1] = Vinf * (ti + Tduct) / p
        end
    end
    return eta_inv
end

function get_ideal_efficiency(total_thrust, rhoinf, Vinf, Rref)
    if Vinf != 0.0
        # TC = total_thrust / (0.5 * rhoinf * Vinf^2 * pi * Rref^2)
        return 2.0 / (
            1.0 +
            sqrt(max((total_thrust / (0.5 * rhoinf * Vinf^2 * pi * Rref^2)), -1.0) + 1.0)
        )
    else
        return 0.0
    end
end

function tqpcoeff(thrust, torque, power, rhoinf, Omega, Rref)
    T = promote_type(eltype(thrust), eltype(torque), eltype(Omega))
    CT = zeros(T, length(Omega))
    CQ = zeros(T, length(Omega))
    CP = zeros(T, length(Omega))
    return tqpcoeff!(CT, CQ, CP, thrust, torque, power, rhoinf, Omega, Rref)
end

function tqpcoeff!(CT, CQ, CP, thrust, torque, power, rhoinf, Omega, Rref)
    for (i, o) in enumerate(Omega)
        if isapprox(o, 0.0)
            CT[i] = CQ[i] = CP[i] = 0.0
        else
            # reference diameter
            # D = 2.0 * Rref

            # rototion in rev per second
            # n = o / (2.0 * pi)

            # thrust coefficient
            CT[i] = thrust[i] / (rhoinf * (o / (2.0 * pi))^2 * (2.0 * Rref)^4)

            # torque coefficient
            CQ[i] = torque[i] / (rhoinf * (o / (2.0 * pi))^2 * (2.0 * Rref)^5)

            # power coefficient
            CP[i] = power[i] / (rhoinf * (o / (2.0 * pi))^3 * (2.0 * Rref)^5)
        end
    end

    return CT, CQ, CP
end

function get_blade_loads(Wmag_rotor, beta1, cl, cd, chords, rhoinf)#, Rhub, Rtip, rotor_panel_centers,rotor_panel_lengths ,B, Omega)
    # initialize
    Np = similar(Wmag_rotor) .= 0.0
    Tp = similar(Wmag_rotor) .= 0.0
    return get_blade_loads!(Np, Tp, Wmag_rotor, beta1, cl, cd, chords, rhoinf)#, Rhub, Rtip, rotor_panel_centers,rotor_panel_lengths ,B, Omega)
end

function get_blade_loads!(Np, Tp, Wmag_rotor, beta1, cl, cd, chords, rhoinf, cache)#, Rhub, Rtip, rotor_panel_centers,rotor_panel_lengths ,B, Omega)

    # dimensions
    nr, nrotor = size(Np)

    for irotor in 1:nrotor
        for ir in 1:nr
            # rename for convenience
            cache.cphi[ir, irotor] = cos(beta1[ir, irotor])
            cache.sphi[ir, irotor] = sin(beta1[ir, irotor])

            # resolve lift and drag into normal and tangential coefficients
            cache.cn[ir, irotor] =
                cl[ir, irotor] * cache.cphi[ir, irotor] -
                cd[ir, irotor] * cache.sphi[ir, irotor]
            cache.ct[ir, irotor] =
                cl[ir, irotor] * cache.sphi[ir, irotor] +
                cd[ir, irotor] * cache.cphi[ir, irotor]

            # get the normal and tangential loads per unit length N' and T'
            Np[ir, irotor] =
                cache.cn[ir, irotor] *
                0.5 *
                rhoinf *
                Wmag_rotor[ir, irotor]^2 *
                chords[ir, irotor]
            Tp[ir, irotor] =
                cache.ct[ir, irotor] *
                0.5 *
                rhoinf *
                Wmag_rotor[ir, irotor]^2 *
                chords[ir, irotor]
        end
    end

    # Npfull = [zeros(nrotor)'; Np; zeros(nrotor)']
    # Tpfull = [zeros(nrotor)'; Tp; zeros(nrotor)']

    ## -- Integrate Loads to get Thrust and Torque
    # add hub/tip for complete integration.  loads go to zero at hub/tip.
    # rfull = [Rhub; rotor_panel_centers; Rtip]

    # thrust and torqe distributions
    # thrust = Npfull
    # torque = Tpfull .* rfull

    # integrate Thrust and Torque (trapezoidal)
    # T = B * fm.trapz(rfull, thrust)
    # Q = B * fm.trapz(rfull, torque)
    # - Actually use rectangle rather than trapezoid integration
    # T = B * sum(rotor_panel_lengths.*Np)
    # Q = B * sum(rotor_panel_lengths.* Tp.*rotor_panel_centers)
    # P = Q * Omega

    return Np, Tp
end
