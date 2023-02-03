#=

Functions regarding wake aerodynamics

=#

"""
DONE. HAS TEST. clean up
"""
function calculate_enthalpy_jumps(Gammas, blade_elements)

    # - Rename for Convenience - #
    nr = length(Gammas[1, :])
    nbe = length(Gammas[:, 1])

    # - Initialize output - #
    H_tilde = similar(Gammas)

    # - Loop through rotors - #
    for i in 1:nr
        @. H_tilde[:, i] =
            blade_elements[i].omega * blade_elements[i].num_blades * Gammas[:, i] /
            (2.0 * pi)
    end

    # - Return cumulative sum of Enthalpy Jumps  - #
    return cumsum(H_tilde; dims=2)
end

"""
DONE. HAS TEST. clean up
"""
function calculate_net_circulation(Gammas, blade_elements)

    # - Rename for Convenience - #
    nr = length(Gammas[1, :])
    nbe = length(Gammas[:, 1])

    # - Initialize output - #
    Gamma_tilde = similar(Gammas)
    BGamma = similar(Gammas)

    # - Loop through rotors - #
    for i in 1:nr
        @. BGamma[:, i] = blade_elements[i].num_blades * Gammas[:, i]
    end

    # - Return cumulative sum of net circulations - #
    return BGamma, cumsum(BGamma; dims=2)
end

"""
DONE. HAS TEST. clean up
"""
function get_surface_velocity(body_vortex_strengths, duct_panels, wake_panels)
    ## -- Set up Edge Velocities -- ##
    x_edge = wake_panels[1].panel_center[:, 1]
    _, duct_inner_idx = findmin(duct_panels.panel_center[:, 1])
    v_surf = -reverse(body_vortex_strengths[1:duct_inner_idx])
    return fm.akima(reverse(duct_panels.panel_center[1:duct_inner_idx, 1]), v_surf, x_edge)
end

"""
DONE. HAS TEST. clean up
"""
function calculate_wake_velocities(
    x_edge, edge_velocity, rotoridxs, Gamma_tilde, H_tilde, blade_elements
)

    # - Rename for Convenience - #
    nbe = length(Gamma_tilde[:, 1])
    nx = length(x_edge)
    TF = eltype(Gamma_tilde)

    # - Initialize the output vector - #
    Vm = zeros(TF, nbe, nx)
    Vm[nbe, :] .= edge_velocity

    if TF != Float64
        # println((p -> p.value).(H_tilde))
        # println((p -> p.value).(Gamma_tilde))
        # println((p -> p.value).(edge_velocity))
    else
        # println(H_tilde)
        # println(Gamma_tilde)
        # println(edge_velocity)
    end

    for j in 1:nx
        for i in (nbe - 1):-1:1
            r = findlast(x -> x <= x_edge[j], x_edge[rotoridxs])
            # if TF != Float64
            #     println(
            #         (
            #             p -> p.value
            #         ).(
            #             Vm[i + 1, j]^2 +
            #             (1.0 / (2.0 * pi * blade_elements[r].radial_positions[i]))^2 *
            #             (Gamma_tilde[i + 1, r]^2 - Gamma_tilde[i, r]^2) +
            #             2.0 * (H_tilde[i + 1, r] - H_tilde[i, r]),
            #         ),
            #     )
            # else
            #     println(
            #         Vm[i + 1, j]^2 +
            #         (1.0 / (2.0 * pi * blade_elements[r].radial_positions[i]))^2 *
            #         (Gamma_tilde[i + 1, r]^2 - Gamma_tilde[i, r]^2) +
            #         2.0 * (H_tilde[i + 1, r] - H_tilde[i, r])
            #     )
            # end
            Vm[i, j] = sqrt(
                Vm[i + 1, j]^2 +
                (1.0 / (2.0 * pi * blade_elements[r].radial_positions[i]))^2 *
                (Gamma_tilde[i + 1, r]^2 - Gamma_tilde[i, r]^2) -
                2.0 * (H_tilde[i + 1, r] - H_tilde[i, r]),
            )
        end
    end

    return Vm
end

"""
DONE. HAS TEST. clean up
assumes circulation is zero at hub.
"""
function calculate_wake_vorticity(vm)
    wake_gammas = similar(vm)
    wake_gammas[2:end, :] = vm[2:end, :] .- vm[1:(end - 1), :]
    wake_gammas[1, :] = vm[1, :]
    return wake_gammas
end
