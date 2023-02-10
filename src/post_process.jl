#=

Functions Related to Post Processing of Raw Solution

=#

"""
"""
function general_post_processing(converged_states, params)

    # - Extract values of states - #

    # body_vortex_strengths, wake_vortex_strengths, rotor_circulation_strengths, blade_element_source_strengths = extract_state_variables(
    # converged_states, params
    # )

    nbe = length(params.blade_elements[1].radial_positions)
    nr = params.num_rotors
    nx = params.num_wake_x_panels
    nbp = params.num_body_panels
    n_g_bw = nbp + (nbe - 1) * nx #number of gammas for bodies and wake
    ng = n_g_bw + nbe * nr #number of gammas and Gammas

    #  Body gamma Indices
    body_idx = 1:nbp

    # Wake gamma_theta Indices
    wake_gamma_idx = (nbp + 1):(nbp + (nbe - 1) * nx)

    # Rotor Gamma Indices
    rotor_Gamma_idx = (n_g_bw + 1):(ng)

    # - Extract State Variables - #
    body_vortex_strengths = view(converged_states, body_idx)
    wake_vortex_strengths = reshape(view(converged_states, wake_gamma_idx), (nbe - 1, nx))
    rotor_circulation_strengths = reshape(
        view(converged_states, rotor_Gamma_idx), (nbe, nr)
    )

    # - Body Surface Velocity - #
    duct_control_points = params.body_panels[1].panel_center
    hub_control_points = params.body_panels[2].panel_center
    duct_surface_velocity = body_vortex_strengths[1:length(duct_control_points[:, 1])]
    hub_surface_velocity = body_vortex_strengths[(length(duct_control_points[:, 1]) + 1):end]

    return duct_control_points,
    duct_surface_velocity, hub_control_points,
    hub_surface_velocity
end
