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

# """
# DONE. HAS TEST. clean up
# """
# function get_surface_velocity(body_vortex_strengths, duct_panels, wake_panels)
#     ## -- Set up Edge Velocities -- ##
#     x_edge = wake_panels[1].panel_center[:, 1]
#     _, duct_inner_idx = findmin(duct_panels.panel_center[:, 1])
#     v_surf = -reverse(body_vortex_strengths[1:duct_inner_idx])
#     return fm.akima(reverse(duct_panels.panel_center[1:duct_inner_idx, 1]), v_surf, x_edge)
# end

"""
"""
function calculate_wake_velocities(
    A_body_to_wake, #vector of matrices size 1 x num wakes
    gamma_body,
    # A_rotor_to_wake, matrix of matrices size num rotors x num wakes
    # rotor_source_strengths,
    A_wake_to_wake, #matrix of matrices size num wakes x num wakes
    gamma_wake,
)

    # - Rename for Convenience - #
    nr = length(gamma_wake[:, 1])
    nx = length(gamma_wake[1, :])

    # - Initialize - #
    Vm = similar(gamma_wake)

    for i in 1:nr
        # - Add Body Induced Velocities - #
        Vm[i, :] .= A_body_to_wake[i] * gamma_body

        # - Add Wake Induced Velocities - #
        #
        for w in 1:nr
            Vm[i, :] .+= A_wake_to_wake[w, i] * gamma_wake[w, :]
        end

        # # - Add Rotor Induced Velocities - #
        # for w in 1:length(rotor_source_strengths[:, 1])
        #     Vm[:, i] .+= A_rotor_to_wake[w, i] * rotor_source_strengths[w, :]
        # end

    end

    return Vm
end

"""
"""
function calculate_wake_vortex_strengths!(
    wake_vortex_strengths, wake_panels, Vm, Vinf, Gamma_tilde, H_tilde, rotoridxs
)

    # - Rename for Convenience - #
    nr = length(Vm[:, 1])
    nx = length(Vm[1, :])

    for r in 1:nr
        for x in 1:nx

            # - Get the Differences in Gamma_tilde and H_tilde - #

            # find the index of the closest rotor in front of current x position
            i = findlast(idx -> idx <= x, rotoridxs)

            # Gamma_tilde difference
            Gamma_diff = Gamma_tilde[r + 1, i]^2 - Gamma_tilde[r, i]^2

            # H_tilde difference
            H_diff = H_tilde[r + 1, i] - H_tilde[r, i]

            # - Get the Average Meridional Velocities - #

            if r == 1 || r == nr
                # if we are at the first or last station, average velocity is just the velocity on the panel
                Vm_avg = Vm[r, x]

            else
                #otherwise, it is an average of the panels above and below the current panel
                Vm_avg = 0.5 * (Vm[r - 1, x] + Vm[r + 1, x])
            end

            # - Set the Wake Vortex Strength - #

            # keep things positive to avoid division by zero
            # Vm_avg = max(0.1 * Vinf, Vm_avg) #dfdc use 0.1, could probably reduce this somewhat

            if Vm_avg <= 0.0
                # try handling zero/negative velocities like this
                wake_vortex_strengths[r, x] *= 0.0
            else
                wake_vortex_strengths[r, x] =
                    1.0 / (2 * Vm_avg) * (
                        -1.0 / (2 * pi * wake_panels[r].panel_center[x, 2]) * Gamma_diff +
                        2 * H_diff
                    )
            end
        end
    end

    return nothing
end

"""
"""
function initialize_wake_vortex_strengths(
    rotor_circulation_strengths, wake_panels, params, rotoridxs
)

    # - Calculate Enthalpy Jumps - #
    H_tilde = calculate_enthalpy_jumps(rotor_circulation_strengths, params.blade_elements)

    # - Calculate Net Circulation - #
    BGamma, Gamma_tilde = calculate_net_circulation(
        rotor_circulation_strengths, params.blade_elements
    )

    # - Get Vm values from thrust method in dfdc - #
    Vm = vm_from_thrust(
        params.freestream,
        params.blade_elements,
        BGamma,
        Gamma_tilde,
        wake_panels,
        rotoridxs,
    )

    # - Initialize matrix of wake vortex strengths - #
    wake_vortex_strengths = similar(Vm)

    # - Use the normal function for calculating the vortex panel strengths - #
    calculate_wake_vortex_strengths!(
        wake_vortex_strengths,
        wake_panels,
        Vm,
        params.freestream.Vinf,
        Gamma_tilde,
        H_tilde,
        rotoridxs,
    )

    return wake_vortex_strengths
end

"""
DONE. HAS TEST. CLEAN UP.
"""
function vm_from_thrust(
    freestream, blade_elements, BGamma, Gamma_tilde, wake_panels, rotoridxs
)

    # - Initialize - #

    TF = eltype(Gamma_tilde)
    nx = length(wake_panels[1].panel_center[:, 1])
    nw = length(wake_panels)

    # average meridional velocity approxmiation
    Vm = freestream.Vinf * ones(TF, blade_elements[1].num_radial_stations - 1, nx)

    # accumulated axial induced velocity
    vi = 0.0

    for x in 1:nx
        r = findlast(idx -> idx <= x, rotoridxs)

        ## -- Calculate Total Thrust of the Rotor -- ##
        T = 0.0

        for j in 1:nw

            # - Blade section length - #
            dr =
                blade_elements[r].radial_positions[j + 1] -
                blade_elements[r].radial_positions[j]

            # - Blade induced velocity - #

            gt_avg = (Gamma_tilde[j + 1, r] + Gamma_tilde[j, r]) / 2.0
            v_theta = gt_avg / (2.0 * pi * wake_panels[j].panel_center[x, 2])

            # - Absolute Velocity - #
            w_theta = v_theta - blade_elements[r].omega * wake_panels[j].panel_center[x, 2]

            # - Thrust contribution - #
            bg_avg = (BGamma[j + 1, r] + BGamma[j, r]) / 2.0
            T -= freestream.rho * bg_avg * w_theta * dr
        end

        ## -- Calculate axial induced velocity of rotor -- ##

        #disk area TODO should this be the annular area??
        A = blade_elements[r].radial_positions[end]^2 * pi

        # thrust induced velocity squared over 2
        vt2_2 = T / (freestream.rho * A)

        #axial induced velocity
        #includes contributions from prior rotors
        vi += sqrt(((freestream.Vinf) / 2.0)^2 + vt2_2)
        vx = vi - freestream.Vinf / 2.0

        # average meridional velocity behind rotor
        Vm[:, x] .= vx + freestream.Vinf
    end

    return Vm
end

