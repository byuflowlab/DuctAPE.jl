"""
"""
function steady_cp(vs, vinf, vref)
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
function calculate_delta_cp_TE(Gamr, sigr, Vm_rotor, Omega, B, body_r, wake_r, Vref)

    ## -- Calculate change in pressure coefficient -- ##

    # - Calculate net circulations - #
    Gamma_tilde = calculate_net_circulation(Gamr, B)

    # - Calculate enthalpy disk jump - #
    Htilde = calculate_enthalpy_jumps(Gamr, Omega, B)

    # - Calculate entropy disk jump - #
    Stilde = calculate_entropy_jumps(sigr, Vm_rotor[:, end])

    # - Get the tangential velocities on the bodies - #
    Vtheta_body = calculate_v_theta(Gamma_tilde[end], body_r)
    Vtheta_wake = calculate_v_theta(Gamma_tilde[end], wake_r)

    # assemble change in cp due to enthalpy and entropy behind rotor(s)
    dcp_body = calcualte_delta_cp(Htilde[end], Stilde[end], Vtheta_body[end], Vref)
    dcp_wake = calcualte_delta_cp(Htilde[end], Stilde[end], Vtheta_wake[end], Vref)

    return dcp_body, dcp_wake
end

"""
"""
function cp_residual(
    gamb, gamw, Gamr, sigr, Vm_rotor, Vinf, Vref, B, Omega, body_r, wake_r, ndpan
)
    ## -- get "surface" velocity at duct TE wake panel -- ##
    # - get total vx and vr on wake panel - #
    # - get surface velocity based on panel angle - #

    # - get steady pressure coefficient values - #
    cpductouterTE = steady_cp(gamb[ndpan], Vinf, Vref)
    cpductinnerTE = steady_cp(gamb[1], Vinf, Vref)
    cpductwake = steady_cp(vs_wake, Vinf, Vref)

    # - Calculate the change in Cp on the walls due to enthalpy, entropy, and vtheta - #
    deltacpduct, deltacpwake = calculate_delta_cp_TE(
        Gamr, sigr, Vm_rotor, Vinf, Omega, B, body_r, wake_r, Vref
    )

    # - add raw and adjusted cp values together - #
    cpductinnerTE += deltacpduct
    cpductwake += deltacpwake

    #= TODO: how to apply kutta condition? need to  have the trailing edge pressures sum to zero? but why in DFDC are they all the same value?
      idea 1: return sum of values
      idea 2: add 2 equations: cpinner = cpouter and cpouter = cpwake
      idea 3: repace 1st eqn in linear solve with gamb = gamw, then set cpinner = cpouter here
    =#
    return cpductinnerTE, cpductouterTE, cpductwake
end
