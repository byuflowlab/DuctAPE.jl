"""
This is the function being solved
"""
function rotor_only_residual!(F, inputs, params)

    # Initialize outputs
    states = deepcopy(inputs)

    # - Calculated Updated Gamma and Sigma Values - #
    update_rotor_states!(states, params)

    # - Return Difference in State Variable Values - #
    @. F = inputs - states

    return nothing
end

"""
 This function wraps the residual function in order to allow for additional parameters as inputs
 params.converged is updated in place in this function.
 """
function solve_rotor_only(inputs, params)

    # - Define closure that allows for parameters - #
    rwrap(F, inputs) = rotor_only_residual!(F, inputs, params)

    # - Call NLsolve function using AD for Jacobian - #
    #= res is of type NLsolve.SolverResults.
    The zero field contains the "solution" to the non-linear solve.
    The converged() function tells us if the solver converged.
    =#
    res = NLsolve.nlsolve(
        rwrap,
        inputs;
        autodiff=:forward,
        method=:newton,
        # iterations=25, #keep iterations low for initial testing/plotting
        iterations=100, #keep iterations low for initial testing/plotting
        show_trace=true,
        # linesearch=BackTracking(; maxstep=1e6),
    )

    println("converged? ", converged(res))
    # - Overwrite the convergence flag in the parameters - #
    params.converged[1] = converged(res)

    # - Return the values driving the residual to zero
    return res.zero
end

function update_rotor_states!(states, params)

    # - Unpack - #
    Gamma, gamma_theta, sigma = extract_rotor_states(states, params)
    # Gamma, gamma_theta = extract_rotor_states(states, params)

    wake_vortex_strengths = repeat(gamma_theta; inner=(1, params.nxwake))

    rpc = params.rotor_panel_centers

    TF = eltype(Gamma)

    # - get the induced velocities at the rotor plane - #
    vx_rotor, vr_rotor, vtheta_rotor = calculate_induced_velocities_on_rotors(
        params.blade_elements,
        Gamma,
        params.vx_rw,
        params.vr_rw,
        wake_vortex_strengths,
        params.vx_rr,
        params.vr_rr,
        sigma,
    )

    # the axial component also includes the freestream velocity ( see eqn 1.87 in dissertation)
    Wx_rotor = vx_rotor .+ params.Vinf
    # the tangential also includes the negative of the rotation rate (see eqn 1.87 in dissertation)
    Wtheta_rotor = vtheta_rotor .- params.blade_elements[1].Omega .* rpc

    Wm_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2)

    # - Get the inflow magnitude at the rotor as the combination of all the components - #
    Wmag_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2 .+ Wtheta_rotor .^ 2)

    # - Update Gamma - #
    calculate_gamma_sigma!(
        Gamma,
        # similar(Gamma) .= 0,
        sigma,
        params.blade_elements,
        Wm_rotor,
        # Wx_rotor,
        Wtheta_rotor,
        Wmag_rotor,
    )

    Gamma_tilde = calculate_net_circulation(Gamma, params.blade_elements[1].B)
    H_tilde = calculate_enthalpy_jumps(
        Gamma, params.blade_elements[1].Omega, params.blade_elements[1].B
    )

    # - update wake strengths - #
    calculate_wake_vortex_strengths!(
        gamma_theta, params.rotor_panel_edges, Wm_rotor, Gamma_tilde, H_tilde
    )

    return nothing
end

function extract_rotor_states(states, params)

    # Problem Dimensions
    nrotor = params.num_rotors                             # number of rotors
    nr = length(params.blade_elements[1].rbe) # number of radial coordinates

    # State Variable Indices
    iGamr = 1:(nr * nrotor) # rotor circulation strength indices
    igamw = (iGamr[end] + 1):(iGamr[end] + nrotor * (nr + 1))     # wake vortex strength indices
    isigr = (igamw[end] + 1):(igamw[end] + nr * nrotor) # rotor source strenght indices

    # Extract State variables
    Gamr = reshape(view(states, iGamr), (nr, nrotor)) # rotor circulation strengths
    gamw = reshape(view(states, igamw), (nr + 1, nrotor)) # wake circulation strengths
    sigr = reshape(view(states, isigr), (nr, nrotor)) # rotor source strengths

    return Gamr, gamw, sigr
    # return Gamr, gamw
end
