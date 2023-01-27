#=

Functions regarding wake aerodynamics

=#

"""
"""
function calculate_enthalpy_jumps(Gammas, Omegas, Bs)

    # - Rename for Convenience - #
    nr = length(Gammas[1, :])
    nbe = length(Gammas[:, 1])

    # - Initialize output - #
    H_tilde = similar(Gammas)

    # - Loop through rotors - #
    for i in 1:nr
        H_tilde[:, i] .= Omegas[i] * Bs[i] * Gammas[i] / (2.0 * pi)
    end

    # - Return cumulative sum of Enthalpy Jumps  - #
    return cumsum(H_tilde; dims=2)
end

"""
"""
function calculate_net_circulation(Gammas, Bs)

    # - Rename for Convenience - #
    nr = length(Gammas[1, :])
    nbe = length(Gammas[:, 1])

    # - Initialize output - #
    Gamma_tilde = similar(Gammas)

    # - Loop through rotors - #
    for i in 1:nr
        Gamma_tilde[:, i] .= Bs[i] * Gammas[:, i]
    end

    # - Return cumulative sum of net circulations - #
    return cumsum(Gamma_tilde; dims=2)
end

"""
"""
function calculate_meridional_velocities(V_edge, Gamma_tilde, H_tilde, radial_position)

    # - Rename for Convenience - #
    nr = length(Gammas[1, :])
    nbe = length(Gammas[:, 1])
    TF = eltype(Gamma_tilde)

    # - Initialize the output vector - #
    vm = V_edge * ones(TF, nbe)

    # - Loop through the radial stations - #
    for i in (nbe - 1):-1:1
        vm[i] = sqrt(
            vm[i + 1]^2 +
            (1.0 / (2.0 * pi * radial_position[i]))^2 *
            (Gamma_tilde[i]^2 - Gamma_tilde[i + 1]^2) +
            2.0 * (H_tilde[i]^2 - H_tilde[i + 1]^2),
        )
    end

    return vm
end

"""
"""
function calculate_wake_vorticity() end

"""
"""
function solve_linear_system() end
