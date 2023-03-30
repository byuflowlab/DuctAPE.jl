"""
This is the function being solved
"""
function residual!(F, state_variables, params)

    # println("input states")
    # println(state_variables[1:17])

    # Initialize outputs
    updated_states = deepcopy(state_variables)

    # - Calculated Updated Gamma and Sigma Values - #
    update_gamma_sigma!(updated_states, params)

    # println("updated_states: ")
    # println(updated_states[1:17])

    # - Return Difference in State Variable Values - #
    @. F = state_variables - updated_states

    return nothing
end

"""
 This function wraps the residual function in order to allow for additional parameters as inputs
 params.converged is updated in place in this function.
 """
function solve!(state_variables, params)

    # - Define closure that allows for parameters - #
    rwrap(F, state_variables) = residual!(F, state_variables, params)

    # - Call NLsolve function using AD for Jacobian - #
    #= res is of type NLsolve.SolverResults.
    The zero field contains the "solution" to the non-linear solve.
    The converged() function tells us if the solver converged.
    =#
    res = NLsolve.nlsolve(
        rwrap,
        state_variables;
        autodiff=:forward,
        method=:newton,
        iterations=25, #keep iterations low for initial testing/plotting
        show_trace=true,
        linesearch=BackTracking(; maxstep=1e6),
    )

    # - Overwrite the convergence flag in the parameters - #
    params.converged[1] = converged(res)

    # - Return the values driving the residual to zero
    return res.zero
end

function update_gamma_sigma!(states, params)

    # - Unpack - #
    Gamma, gamma_theta = extract_state_variables(states, params)

    wake_vortex_strengths = repeat(gamma_theta; inner=(1, params.nxwake))

    rpc = params.rotor_panel_centers

    TF = eltype(Gamma)

    # - get the induced velocities at the rotor plane - #
    vx_rotor, vr_rotor, vtheta_rotor = dt.calculate_induced_velocities_on_rotors(
        params.blade_elements, Gamma, params.vx_rw, params.vr_rw, wake_vortex_strengths
    )

    # the axial component also includes the freestream velocity ( see eqn 1.87 in dissertation)
    Wx_rotor = vx_rotor .+ params.Vinf
    # the tangential also includes the negative of the rotation rate (see eqn 1.87 in dissertation)
    Wtheta_rotor = vtheta_rotor .- params.Omega .* rpc

    Wm_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2)

    # - Get the inflow magnitude at the rotor as the combination of all the components - #
    Wmag_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2 .+ Wtheta_rotor .^ 2)

    # - Update Gamma - #
    # TODO: check that this is outputting as expected. make sure gammas and sigmas aren't flipped.
    # TODO: need to add/update unit test
    dt.calculate_gamma_sigma!(
        Gamma,
        similar(states[params.Gammaidx]),
        params.blade_elements,
        Wm_rotor,
        Wtheta_rotor,
    )

    #TODO: need to unit test these functions
    Gamma_tilde = dt.calculate_net_circulation(Gamma, params.num_blades)
    H_tilde = dt.calculate_enthalpy_jumps(
        Gamma, params.Omega, params.num_blades
    )

    # - update wake strengths - #
    # TODO: need to unit test this function
    dt.calculate_wake_vortex_strengths!(
        gamma_theta,
        params.rotor_panel_edges,
        Wm_rotor,
        Gamma_tilde,
        H_tilde,
    )

    return nothing
end

function extract_state_variables(states, params)

    # Problem Dimensions
    nrotor = params.num_rotors                             # number of rotors
    nr = length(params.blade_elements[1].radial_positions) # number of radial coordinates

    # State Variable Indices
    iΓr = 1:(nr * nrotor) # rotor circulation strength indices
    iΓw = (iΓr[end] + 1):(iΓr[end] + nrotor * (nr + 1))     # wake vortex strength indices

    # Extract State variables
    Γr = reshape(view(states, iΓr), (nr, nrotor)) # rotor circulation strengths
    Γw = reshape(view(states, iΓw), (nr + 1, nrotor)) # wake circulation strengths

    return Γr, Γw
end
