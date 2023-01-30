#=

Functions regarding wake aerodynamics

=#

"""
DONE. clean up and check
"""
function calculate_enthalpy_jumps(Gammas, blade_elements)

    # - Rename for Convenience - #
    nr = length(Gammas[1, :])
    nbe = length(Gammas[:, 1])

    # - Initialize output - #
    H_tilde = similar(Gammas)

    # - Loop through rotors - #
    for i in 1:nr
        H_tilde[:, i] .=
            blade_elements[i].omega * blade_elements[i].B * Gammas[i] / (2.0 * pi)
    end

    # - Return cumulative sum of Enthalpy Jumps  - #
    return cumsum(H_tilde; dims=2)
end

"""
DONE. clean up and check
"""
function calculate_net_circulation(Gammas, blade_elements)

    # - Rename for Convenience - #
    nr = length(Gammas[1, :])
    nbe = length(Gammas[:, 1])

    # - Initialize output - #
    Gamma_tilde = similar(Gammas)

    # - Loop through rotors - #
    for i in 1:nr
        BGamma[:, i] .= blade_elements[i].B * Gammas[:, i]
    end

    # - Return cumulative sum of net circulations - #
    return BGamma, cumsum(BGamma; dims=2)
end

"""
DONE. clean up and check
"""
function calculate_meridional_velocities(
    body_vortex_strengths,
    duct_panels,
    wake_grid,
    rotoridxs,
    Gamma_tilde,
    H_tilde,
    blade_elements,
)

    # - Rename for Convenience - #
    nr = length(rotoridxs)
    nbe = length(wake_grid[:, 1])
    nx = length(wake_grid[1, :])
    TF = eltype(Gamma_tilde)

    # - Initialize the output vector - #
    vm = zeros(TF, nbe, nx)

    ## -- Set up Edge Velocities -- ##
    x_edge = getindex.(wake_grid[1, :], 1)
    v_surf = reverse(body_vortex_strengths[1, length(duct_panels)])
    vm[1, :] .= fm.akima(getindex.(duct_panels.panel_center, 1), v_surf, x_edge)

    for i in (nbe - 1):-1:1
        for j in 1:nx
            r = findlast(x -> x <= x_edge[j], x_edge[rotoridxs])
            vm[i, j] = sqrt(
                vm[i + 1, j]^2 +
                (1.0 / (2.0 * pi * blade_elements[r].radial_positions[i]))^2 *
                (Gamma_tilde[i + 1, r]^2 - Gamma_tilde[i, r]^2) +
                2.0 * (H_tilde[i + 1, r]^2 - H_tilde[i + 1, r]^2),
            )
        end
    end

    return vm
end

"""
DONE. clean up and check
"""
function calculate_wake_vorticity(vm)
    wake_gammas = similar(vm)
    wake_gammas[2:end, :] = vm[2:end, :] .- vm[1:(end - 1), :]
    wake_gammas[1, :] = vm[1, :]
    return wake_gammas
end
