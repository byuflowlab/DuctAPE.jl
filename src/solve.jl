"""
    analyze_propulsor(x, fx=x->x; tol=1e-8, maxiter=100)

Version of `analyze_propeller` designed for sensitivity analysis.  `x` is an input vector
and `fx` is a function which returns the standard input arguments to `analyze_propeller`.
"""
function analyze_propulsor(x, fx=x -> x; tol=1e-8, maxiter=100)

    # convergence flag
    converged = [false]

    # define parameters
    p = (; fx, tol, maxiter, converged)

    # compute state variables (updates convergence flag internally)
    states = ImplicitAD.implicit(solve, residual!, x, p)

    # TODO: post-processing using the converged state variables

    # return solution
    return states, converged[1]
end

"""
    analyze_propulsor(duct_coordinates, hub_coordinates, rotor_parameters, freestream;
        tol=1e-8, maxiter=100)

Finds a converged set of circulation and source and vortex strengths for a ducted propeller.
"""
function analyze_propulsor(
    duct_coordinates,
    hub_coordinates,
    paneling_constants,
    rotor_parameters,
    freestream;
    tol=1e-8,
    maxiter=100,
)

    # use empty input vector
    x = Float64[]

    # use default input function
    fx =
        x -> (;
            duct_coordinates,
            hub_coordinates,
            paneling_constants,
            rotor_parameters,
            freestream,
        )

    # convergence flag
    converged = [false]

    # define parameters
    p = (; fx, tol, maxiter, converged)

    # compute various panel and circulation strenghts (updates convergence flag internally)
    strengths = solve(x, p)

    # TODO: post-processing using the converged strengths

    # return solution
    return strengths, p.converged[1]
end

"""
    solve(x, p)

"""
function solve(x, p)

    # unpack parameters
    (; fx, tol, maxiter, converged) = p

    # unpack inputs
    (; duct_coordinates, hub_coordinates, paneling_constants, rotor_parameters, freestream) = fx(
        x
    )

    # initialize various inputs used in analysis
    inputs = precomputed_inputs(
        duct_coordinates, hub_coordinates, paneling_constants, rotor_parameters, freestream
    )

    # calculate initial guess for state variables
    states = initialize_states(inputs)

    # initialize residual vector
    resid = copy(states)

    # - Define closure that allows for parameters - #
    rwrap(r, states) = residual!(r, states, inputs, p)

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
        iterations=maxiter, #25, #keep iterations low for initial testing/plotting
        show_trace=true,
        # linesearch=BackTracking(; maxstep=1e6),
        linesearch=BackTracking(; maxstep=tol),
    )

    # save convergence information
    converged[1] = converged(res)

    # return solution
    return res.zero
end

"""
    residual!(res, states, inputs, p)

Updates the state variables.
"""
function residual!(res, states, inputs, p)

    # - Extract commonly used items from precomputed inputs - #
    blade_elements = inputs.blade_elements
    rpc = inputs.rotor_panel_centers
    Vinf = inputs.Vinf

    # - Extract States - #
    gamb, gamw, Gamr, sigr = extract_state_variables(states, inputs)

    # - Fill out wake strengths - #
    wake_vortex_strengths = fill_out_wake_strengths(
        gamw, inputs.rotor_indices, inputs.num_wake_x_panels
    )

    # - Get the induced velocities at the rotor plane - #
    vx_rotor, vr_rotor, vtheta_rotor = calculate_induced_velocities_on_rotors(
        blade_elements,
        Gamr,
        inputs.vx_rw,
        inputs.vr_rw,
        wake_vortex_strengths,
        inputs.vx_rr,
        inputs.vr_rr,
        sigr,
        inputs.vx_rb,
        inputs.vr_rb,
        gamb,
    )

    # the axial component also includes the freestream velocity ( see eqn 1.87 in dissertation)
    Wx_rotor = vx_rotor .+ inputs.Vinf

    # the tangential also includes the negative of the rotation rate (see eqn 1.87 in dissertation)
    Wtheta_rotor = vtheta_rotor .- inputs.blade_elements[1].Omega .* rpc

    # meridional component
    Wm_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2)

    # Get the inflow magnitude at the rotor as the combination of all the components
    Wmag_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2 .+ Wtheta_rotor .^ 2)

    # - Update rotor circulation and source panel strengths - #
    calculate_gamma_sigma!(Gamr, sigr, blade_elements, Wm_rotor, Wtheta_rotor, Wmag_rotor)

    # - Calculate net circulation and enthalpy jumps - #
    # TODO: check that your get property override works here for inputting an array of number of blades and rotation rates
    Gamma_tilde = calculate_net_circulation(Gamr, blade_elements.B)
    H_tilde = calculate_enthalpy_jumps(Gamr, blade_elements.Omega, blade_elements.B)

    # - update wake strengths - #
    calculate_wake_vortex_strengths!(
        gamw, inputs.rotor_panel_edges, Wm_rotor, Gamma_tilde, H_tilde
    )

    # - Calculate body vortex strengths - #
    calculate_body_vortex_strengths!(
        gamb,
        inputs.A_bb,
        inputs.b_bf,
        inputs.kutta_idxs,
        inputs.A_bw,
        wake_vortex_strengths,
        inputs.A_br,
        sigr,
    )

    # Update Residual
    res .= states .- res

    return nothing
end

"""
    extract_state_variables(states, inputs)

Extract circulation and source strengths from the state vector
"""
function extract_state_variables(states, inputs)

    # Problem Dimensions
    nb = inputs.num_body_panels                     # number of body panels
    nr, nrotor = size(inputs.rotor_panel_centers)   # number of blade elements and rotors
    nw = nr + 1                                     # number of wake sheets

    # State Variable Indices
    iΓb = 1:nb                                      # body vortex strength indices
    iΓw = (iΓb[end] + 1):(iΓb[end] + nw * nrotor)   # wake vortex strength indices
    iΓr = (iΓw[end] + 1):(iΓw[end] + nr * nrotor)   # rotor circulation strength indices
    iΣr = (iΓr[end] + 1):(iΓr[end] + nr * nrotor)   # rotor source strength indices

    # Extract State variables
    gamb = view(states, iΓb)                        # body vortex strengths
    gamw = reshape(view(states, iΓw), (nw, nrotor)) # wake vortex strengths
    Gamr = reshape(view(states, iΓr), (nr, nrotor)) # rotor circulation strengths
    sigr = reshape(view(states, iΣr), (nr, nrotor)) # rotor circulation strengths

    return gamb, gamw, Gamr, sigr
end
