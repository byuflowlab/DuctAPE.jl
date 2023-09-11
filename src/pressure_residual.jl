"""
"""
function steady_cps(vs, vinf, vref)
    return (vinf^2 .- vs .^ 2) / vref^2
end

"""
"""
function calculate_delta_cp(deltaH, deltaS, Vtheta, Vref)
    return (2.0 * (deltaH - deltaS) .- Vtheta .^ 2) / Vref^2
end

"""
"""
function calculate_v_theta(Gamma_tilde, r)
    return Gamma_tilde ./ (2.0 * pi * r)
end

"""
"""
function calculate_delta_cp_TE(Gamr, sigr, Vm_rotor, Omega, B, r, Vref)

    ## -- Calculate change in pressure coefficient -- ##

    # - Calculate net circulations - #
    Gamma_tilde = calculate_net_circulation(Gamr, B)

    # - Calculate enthalpy disk jump - #
    Htilde = calculate_enthalpy_jumps(Gamr, Omega, B)

    # - Calculate entropy disk jump - #
    Stilde = calculate_entropy_jumps(sigr, Vm_rotor[:, end])

    # - Get the tangential velocities on the bodies - #
    Vtheta_body = calculate_v_theta(Gamma_tilde[end], r)

    # assemble change in cp due to enthalpy and entropy behind rotor(s)
    dcp_body = calculate_delta_cp(Htilde[end], Stilde[end], Vtheta_body[end], Vref)

    return dcp_body
end

"""
"""
function cp_residual(states, inputs)

    # extract states
    gamb, gamw, Gamr, sigr = extract_state_variables(states, inputs)

    # extract needed inputs
    Vinf = inputs.freestream.Vinf
    Vref = inputs.reference_parameters.Vref
    B = inputs.blade_elements.B
    Omega = inputs.blade_elements.Omega

    # get Vm_rotor
    _, _, _, _, _, Vm_rotor, _ = calculate_rotor_velocities(Gamr, gamw, sigr, gamb, inputs)

    # get number of panels in duct
    ndpan = length(inputs.body_panels[1].panel_center[:, 1])
    # get radial location of inner duct TE panel
    body_r = inputs.body_panels[1].panel_center[1, 2]

    # - get steady pressure coefficient values - #
    cpductouterTE = steady_cps(gamb[ndpan], Vinf, Vref)
    cpductinnerTE = steady_cps(gamb[1], Vinf, Vref)

    # - Calculate the change in Cp on the walls due to enthalpy, entropy, and vtheta - #
    deltacpduct = calculate_delta_cp_TE(Gamr, sigr, Vm_rotor, Omega, B, body_r, Vref)

    # - add raw and adjusted cp values together - #
    cpductinnerTE += deltacpduct

    # Note: Kutta Condition should just be that the pressure coefficients need to be equal at the trailing edge. the signs are lost in the v^2, but it shouldn't matter.
    return cpductouterTE - cpductinnerTE
end
