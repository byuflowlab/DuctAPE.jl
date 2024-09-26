"""
    inviscid_rotor_thrust(Ctheta_rotor, Gamma_tilde, rotor_panel_length, rhoinf)

Calculate inviscid rotor thrust.

# Arguments
- `Ctheta_rotor::Vector{Float}` : Absolute tangential velocity on rotor blade elements
- `Gamma_tilde::Matrix{Float}` : net upstream rotor circulation
- `rotor_panel_length::Vector{Float}` : dimensional lengths on which blade element values apply
- `rhoinf::Float` : freestream density

# Returns
- `Tinv::Vector{Float}` : inviscid dimensional thrust
- `dTi::Vector{Float}` : inviscid dimensional thrust distribution
"""
function inviscid_rotor_thrust(Ctheta_rotor, Gamma_tilde, rotor_panel_length, rhoinf)
    # initialize
    dTi = similar(Gamma_tilde) .= 0.0
    Tinv = zeros(eltype(Gamma_tilde), size(Gamma_tilde, 2))

    return inviscid_rotor_thrust!(
        Tinv, dTi, Ctheta_rotor, Gamma_tilde, rotor_panel_length, rhoinf
    )
end

"""
    inviscid_rotor_thrust!(
        Tinv, dTi, Ctheta_rotor, Gamma_tilde, rotor_panel_length, rhoinf
    )

In-place version of `inviscid_rotor_thrust`.
"""
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

"""
    viscous_rotor_thrust(
        Cz_rotor, Cmag_rotor, B, chord, rotor_panel_length, cd, rhoinf
    )

Calculate visous rotor "thrust."

# Arguments
- `Cz_rotor::Vector{Float}` : Absolute axial velocity on rotor blade elements
- `Cmag_rotor::Vector{Float}` : Absolute inflow velocity magnitude on rotor blade elements
- `B::Vector{Float}` : blade count for each rotor (usually integers but could be a float)
- `chord::Vector{Float}` : blade element chord lengths
- `rotor_panel_length::Vector{Float}` : dimensional lengths on which blade element values apply
- `cd::Vector{Float}` : drag coefficient for each blade element
- `rhoinf::Float` : freestream density

# Returns
- `Tvisc::Vector{Float}` : viscous dimensional thrust
- `dTv::Vector{Float}` : viscous dimensional thrust distribution
"""
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

"""
    viscous_rotor_thrust!(
        Tvisc, dTv, Cz_rotor, Cmag_rotor, B, chord, rotor_panel_length, cd, rhoinf
    )

In-place version of `viscous_rotor_thrust`.
"""
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

"""
    inviscid_rotor_torque(
        Cz_rotor, rotor_panel_center, rotor_panel_length, Gamma_tilde, rhoinf
    )

Calculate inviscid rotor torque.

# Arguments
- `Cz_rotor::Vector{Float}` : Absolute axial velocity on rotor blade elements
- `rotor_panel_center::Vector{Float}` : radial location of rotor blade elements
- `rotor_panel_length::Vector{Float}` : dimensional lengths on which blade element values apply
- `Gamma_tilde::Matrix{Float}` : net upstream rotor circulation
- `rhoinf::Float` : freestream density

# Returns
- `Qinv::Vector{Float}` : inviscid dimensional thrust
- `dQi::Vector{Float}` : inviscid dimensional thrust distribution
"""
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

"""
    inviscid_rotor_torque!(
        Qinv, dQi, Cz_rotor, rotor_panel_center, rotor_panel_length, Gamma_tilde, rhoinf
    )

In-place version of `inviscid_rotor_torque`.
"""
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

"""
    viscous_rotor_torque(
        Ctheta_rotor, Cmag_rotor, B, chord, rotor_panel_center, rotor_panel_length, cd, rhoinf
    )

Calculate viscous rotor torque.

# Arguments
- `Ctheta_rotor::Vector{Float}` : Absolute tangential velocity on rotor blade elements
- `Cmag_rotor::Vector{Float}` : Absolute inflow velocity magnitude on rotor blade elements
- `B::Vector{Float}` : blade count for each rotor (usually integers but could be a float)
- `chord::Vector{Float}` : blade element chord lengths
- `rotor_panel_center::Vector{Float}` : radial location of rotor blade elements
- `rotor_panel_length::Vector{Float}` : dimensional lengths on which blade element values apply
- `cd::Vector{Float}` : drag coefficient for each blade element
- `rhoinf::Float` : freestream density

# Returns
- `Qvisc::Vector{Float}` : viscous dimensional thrust
- `dQv::Vector{Float}` : viscous dimensional thrust distribution
"""
function viscous_rotor_torque(
    Ctheta_rotor, Cmag_rotor, B, chord, rotor_panel_center, rotor_panel_length, cd, rhoinf
)
    TF = promote_type(
        eltype(Ctheta_rotor),
        eltype(Cmag_rotor),
        eltype(chord),
        eltype(rotor_panel_center),
        eltype(cd),
    )

    dQv = zeros(TF, size(Cmag_rotor))
    Qvisc = zeros(TF, size(B))

    return viscous_rotor_torque!(
        Qvisc,
        dQv,
        Ctheta_rotor,
        Cmag_rotor,
        B,
        chord,
        rotor_panel_center,
        rotor_panel_length,
        cd,
        rhoinf,
    )
end

"""
    viscous_rotor_torque!(
        Qvisc,
        dQv,
        Ctheta_rotor,
        Cmag_rotor,
        B,
        chord,
        rotor_panel_center,
        rotor_panel_length,
        cd,
        rhoinf
    )

In-place version of `viscous_rotor_torque`.
"""
function viscous_rotor_torque!(
    Qvisc,
    dQv,
    Ctheta_rotor,
    Cmag_rotor,
    B,
    chord,
    rotor_panel_center,
    rotor_panel_length,
    cd,
    rhoinf,
)

    # dimensions
    nr, nrotor = size(dQv)

    for irotor in 1:nrotor
        for ir in 1:nr
            # hrwc = 0.5 * rhoinf * Cmag_rotor[ir, irotor] * chord[ir, irotor]
            # brdr =
            # B[irotor] * rotor_panel_center[ir, irotor] * rotor_panel_length[ir, irotor]
            dQv[ir, irotor] =
                -(0.5 * rhoinf * Cmag_rotor[ir, irotor] * chord[ir, irotor]) *
                cd[ir, irotor] *
                Ctheta_rotor[ir, irotor] *
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

"""
    rotor_power(Q, dQ, Omega)

Calculate power from torque and rotation rate.

# Arguments
- `Q::Vector{Float}` : dimensional thrust
- `dQ::Vector{Float}` : dimensional thrust distribution
- `Omega::Vector{Float}` : rotor rotation rates

# Returns
- `P::Vector{Float}` : dimensional power
- `dP::Vector{Float}` : dimensional thrust distribution
"""
function rotor_power(Q, dQ, Omega)
    dP = similar(dQ) .= 0.0
    P = similar(Q) .= 0.0

    return rotor_power!(P, dP, Q, dQ, Omega)
end

"""
    rotor_power!(P, dP, Q, dQ, Omega)

In-place version of `rotor_power`.
"""
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

"""
    get_total_efficiency(total_thrust, total_power, Vinf)

Get total efficiency.

# Arguments
- `total_thrust::Vector{Float}` : total thrust
- `total_power::Vector{Float}` : total power
- `Vinf::Vector{Float}` : one element vector freestream velocity magnitude

# Returns
- `total_efficiency::Vector{Float} : total efficiency
"""
function get_total_efficiency(total_thrust, total_power, Vinf)
    TF = promote_type(eltype(total_thrust), eltype(total_power), eltype(Vinf))

    eta = zeros(TF, length(total_thrust))

    return get_total_efficiency!(eta, total_thrust, total_power, Vinf)
end

"""
    get_total_efficiency!(eta, total_thrust, total_power, Vinf)

In-place version of `get_total_efficiency`.
"""
function get_total_efficiency!(eta, total_thrust, total_power, Vinf)
    for i in 1:length(total_thrust)
        # if Vinf <= 0.0 || total_power[i] < eps() || total_thrust[i] <= 0.0
        #do nothing, efficiency can't physically be negative or infinite.
        if abs(total_power[i]) < eps()
            # allow negative efficiency
            eta[i] = 0.0
        else
            eta[i] = total_thrust[i] * Vinf / total_power[i]
        end
    end

    return eta
end

"""
    get_induced_efficiency(Tinv, Tduct, Pinv, Vinf)

Get rotor efficiency induced by presence of the duct.

# Arguments
- `Tinv::Vector{Float}` : inviscid dimensional thrust
- `Tduct::Vector{Float}` : duct thrust
- `Pinv::Vector{Float}` : inviscid dimensional power
- `Vinf::Vector{Float}` : one element vector freestream velocity magnitude

# Returns
- `induced_efficiency::Vector{Float}` : rotor efficiency induced by duct
"""
function get_induced_efficiency(Tinv, Tduct, Pinv, Vinf)
    TF = promote_type(eltype(Tinv), eltype(Pinv), eltpye(Tduct))
    eta_inv = zeros(TF, size(Tinv))
    return get_induced_efficiency!(eta_inv, Tinv, Tduct, Pinv, Vinf)
end

"""
    get_induced_efficiency!(eta_inv, Tinv, Tduct, Pinv, Vinf)

In-place version of `get_induced_efficiency`.
"""
function get_induced_efficiency!(eta_inv, Tinv, Tduct, Pinv, Vinf)
    for (e, ti, p) in zip(eachrow(eta_inv), Tinv, Pinv)
        # if Vinf <= 0.0 || p <= 0.0
        if abs(p) <= eps()
            e[1] = 0.0
        else
            e[1] = Vinf * (ti + Tduct) / p
        end
    end
    return eta_inv
end

"""
    get_ideal_efficiency(total_thrust, rhoinf, Vinf, Rref)

Compute ducted fan ideal efficiency

# Arguments
- `total_thrust::Vector{Float}` : total thrust from rotors and duct
- `rhoinf::Float` : freestream density
- `Vinf::Vector{Float}` : one element vector freestream velocity magnitude
- `Rref::Vector{Float}` : one element vector reference rotor tip radius

# Returns
- `ideal_efficiency::Vector{Float}` : ideal ducted fan efficiency
"""
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

"""
    tqpcoeff(thrust, torque, power, rhoinf, Omega, Rref)

Calculate non-dimensional thrust, torque, and power coefficients

# Arguments
- `thrust::Vector{Float}` : dimensional thrust
- `torque::Vector{Float}` : dimensional torque
- `power::Vector{Float}` : dimensional power
- `rhoinf::Float` : freestream density
- `Omega::Vector{Float}` : rotor rotation rates
- `Rref::Vector{Float}` : one element vector reference rotor tip radius

# Returns
- `CT::Vector{Float}` : thrust coefficient
- `CQ::Vector{Float}` : torque coefficient
- `CP::Vector{Float}` : power coefficient
"""
function tqpcoeff(thrust, torque, power, rhoinf, Omega, Rref)
    T = promote_type(eltype(thrust), eltype(torque), eltype(Omega))
    CT = zeros(T, length(Omega))
    CQ = zeros(T, length(Omega))
    CP = zeros(T, length(Omega))
    return tqpcoeff!(CT, CQ, CP, thrust, torque, power, rhoinf, Omega, Rref)
end

"""
    tqpcoeff!(CT, CQ, CP, thrust, torque, power, rhoinf, Omega, Rref)

In-place version of `tqpcoeff`.
"""
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

"""
    get_blade_loads(Cmag_rotor, beta1, cl, cd, chords, rhoinf)

Get loading along blades.

# Arguments
- `Cmag_rotor::Vector{Float}` : blade element inflow magnitudes
- `beta1::Vector{Float}` : blade element inflow angles
- `cl::Vector{Float}` : blade element lift coefficients
- `cd::Vector{Float}` : blade element drag coefficients
- `chords::Vector{Float}` : blade element chord lengths
- `rhoinf::Vector{Float}` : one element freestream density

# Returns
- `Np::Vector{Float}` : blade loading per unit length in the normal direction: N'
- `Tp::Vector{Float}` : blade loading per unit length in the tangential direction: T'
"""
function get_blade_loads(Cmag_rotor, beta1, cl, cd, chords, rhoinf)#, Rhub, Rtip, rotor_panel_centers,rotor_panel_lengths ,B, Omega)
    # initialize
    Np = similar(Cmag_rotor) .= 0.0
    Tp = similar(Cmag_rotor) .= 0.0
    return get_blade_loads!(Np, Tp, Cmag_rotor, beta1, cl, cd, chords, rhoinf)#, Rhub, Rtip, rotor_panel_centers,rotor_panel_lengths ,B, Omega)
end

"""
    get_blade_loads!(Np, Tp, Cmag_rotor, beta1, cl, cd, chords, rhoinf, cache)

In-place version of `get_blade_loads`.
"""
function get_blade_loads!(Np, Tp, Cmag_rotor, beta1, cl, cd, chords, rhoinf, cache)#, Rhub, Rtip, rotor_panel_centers,rotor_panel_lengths ,B, Omega)

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
                Cmag_rotor[ir, irotor]^2 *
                chords[ir, irotor]
            Tp[ir, irotor] =
                cache.ct[ir, irotor] *
                0.5 *
                rhoinf *
                Cmag_rotor[ir, irotor]^2 *
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
