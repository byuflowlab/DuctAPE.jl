"""
This is the function being solved
"""
function rotor_only_residual!(F, init_states, inputs)

    # Initialize outputs
    states = deepcopy(init_states)

    # - Calculated Updated Gamr and sigr Values - #
    update_rotor_states!(states, inputs)

    # - Return Difference in State Variable Values - #
    @. F = init_states - states

    return nothing
end

"""
 This function wraps the residual function in order to allow for additional parameters as inputs
 params.converged is updated in place in this function.
 """
function solve_rotor_only(states, inputs)

    # - Define closure that allows for parameters - #
    rwrap(F, states) = rotor_only_residual!(F, states, inputs)

    # - Call NLsolve function using AD for Jacobian - #
    #= res is of type NLsolve.SolverResults.
    The zero field contains the "solution" to the non-linear solve.
    The converged() function tells us if the solver converged.
    =#
res = NLsolve.nlsolve(
        rwrap,
        states;
        autodiff=:forward,
        method=:newton,
        # iterations=25, #keep iterations low for initial testing/plotting
        iterations=100, #keep iterations low for initial testing/plotting
        show_trace=true,
        # linesearch=BackTracking(; maxstep=1e6),
    )

    println("converged? ", converged(res))
    # - Overwrite the convergence flag in the parameters - #
    inputs.converged[1] = converged(res)

    # - Return the values driving the residual to zero
    return res.zero
end

function update_rotor_states!(states, inputs)

    # - Unpack - #
    Gamr, gamw, sigr = extract_rotor_states(states, inputs)

    rpc = inputs.rotor_panel_centers

    TF = eltype(Gamr)

    # - get the induced velocities at the rotor plane - #
    vx_rotor, vr_rotor, vtheta_rotor = calculate_induced_velocities_on_rotors(
        inputs.blade_elements,
        Gamr,
        inputs.vx_rw,
        inputs.vr_rw,
        gamw,
        inputs.vx_rr,
        inputs.vr_rr,
        sigr,
    )

    Wx_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor = reframe_rotor_velocities(
        vx_rotor, vr_rotor, vtheta_rotor, inputs.Vinf, inputs.blade_elements.Omega, rpc
    )

    vx_wake, vr_wake = calculate_induced_velocities_on_wakes(
        inputs.vx_ww, inputs.vr_ww, gamw, inputs.vx_wr, inputs.vr_wr, sigr
    )

    Wm_wake = reframe_wake_velocities(vx_wake, vr_wake, inputs.Vinf)

    # - Update Gamr - #
    calculate_gamma_sigma!(
        Gamr,
        sigr,
        inputs.blade_elements,
        Wm_rotor,
        Wtheta_rotor,
        Wmag_rotor,
        inputs.freestream,
    )

    Gamma_tilde = calculate_net_circulation(Gamr, inputs.blade_elements[1].B)
    H_tilde = calculate_enthalpy_jumps(
        Gamr, inputs.blade_elements.Omega, inputs.blade_elements.B
    )

    # - update wake strengths - #
    # TODO: update inputs to have wakeK's
    calculate_wake_vortex_strengths!(gamw, Gamr, Wm_wake, inputs)

    return nothing
end

function extract_rotor_states(states, inputs)

    # Problem Dimensions
    nrotor = inputs.num_rotors                             # number of rotors
    nr = inputs.nrotor_panels # number of rotor panels (length for Gamr and sigr)
    nw = inputs.nwake_panels # number of wake panels (length for gamw)

    # State Variable Indices
    iGamr = 1:(nr * nrotor) # rotor circulation strength indices
    igamw = (iGamr[end] + 1):(iGamr[end] + nw)
    isigr = (igamw[end] + 1):(igamw[end] + nr * nrotor) # rotor source strenght indices

    # Extract State variables
    Gamr = reshape(view(states, iGamr), (nr, nrotor)) # rotor circulation strengths
    gamw = view(states, igamw) # wake vortex strengths
    sigr = reshape(view(states, isigr), (nr, nrotor)) # rotor source strengths

    return Gamr, gamw, sigr
    # return Gamr, gamw
end
