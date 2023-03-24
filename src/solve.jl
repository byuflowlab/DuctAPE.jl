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

Finds a converged set of circulation and source strengths for a ducted propeller.
"""
function analyze_propulsor(
    duct_coordinates, hub_coordinates, rotor_parameters, freestream; tol=1e-8, maxiter=100
)

    # use empty input vector
    x = Float64[]

    # use default input function
    fx = x -> (duct_coordinates, hub_coordinates, rotor_parameters, freestream)

    # convergence flag
    converged = [false]

    # define parameters
    p = (; fx, tol, maxiter, converged)

    # compute state variables (updates convergence flag internally)
    states = solve(x, p)

    # TODO: post-processing using the converged state variables

    # return solution
    return states, converged[1]
end

"""
    solve(x, p)

Use fixed point iteration to find a converged set of state variables.
"""
function solve(x, p)

    # unpack parameters
    (; fx, tol, maxiter, converged) = p

    # unpack inputs
    (; duct_coordinates, hub_coordinates, rotor_parameters, freestream) = fx(x)

    # initialize parameters
    params = initialize_parameters(
        duct_coordinates, hub_coordinates, rotor_parameters, freestream
    )

    # calculate initial guess for state variables
    states = calculate_initial_states(params)

    # initialize residual vector
    resid = copy(states)

    # set convergence flag to false
    converged[1] = false

    # perform fixed point iteration
    for iter in 1:maxiter

        # store previous states in the residual vector
        resid .= states

        # update state variables
        update_gamma_sigma!(states, params)

        # calculate difference between updated and original states
        resid .= states .- resid

        # check if all state variables are converged
        if all(x -> abs(x) < tol, resid)
            # set convergence flag to true
            converged[1] = true
            # stop iterating
            break
        end
    end

    # return state variables
    return states
end

"""
    residual!(r, y, x, p)

Calculate the residual function (only used for implicit differentiation).
"""
function residual!(r, y, x, p)

    # unpack parameters
    (; fx, tol, maxiter) = p

    # unpack inputs
    (; duct_coordinates, hub_coordinates, rotor_parameters, freestream) = fx(x)

    # initialize parameters
    params = initialize_parameters(
        duct_coordinates, hub_coordinates, rotor_parameters, freestream
    )

    # copy states to the residual vector
    r .= y

    # update states which are stored in the residual vector
    update_gamma_sigma!(r, params)

    # residual is the difference between the updated and original states
    r .-= y

    # return residual
    return r
end

"""
    calculate_initial_states(params)

Calculate an initial guess for the state variables
"""
function calculate_initial_states(params)

    # initialize body vortex strengths (no-rotor linear problem)
    A = params.A_bb # AIC matrix for body to body problem
    b = params.b_fs # right hand side for body to body problem
    Γb = A \ b # get circulation strengths from solving body to body problem

    # initialize blade circulation and source strengths (assume no body influence)
    Γr, Σr = calculate_gamma_sigma(params.blade_elements, params.freestream.Vinf)

    # initialize wake vortex strengths
    Γw = initialize_wake_vortex_strengths(Γr, Σr, params)

    # combine initial states into one vector
    states = vcat(
        Γb, # body vortex sheet strengths
        reduce(vcat, Γw), # wake vortex sheet strengths
        reduce(vcat, Γr), # rotor vortex strengths
        reduce(vcat, Σr), # rotor source strengths
    )

    return states
end

"""
    update_gamma_sigma!(states, params)

Updates the state variables.
"""
function update_gamma_sigma!(states, params)

    # extract parameters
    blade_elements = params.blade_elements
    wake_panels = params.wake_panels
    b_fb = params.bc_freestream_to_body
    Vinf = params.freestream.Vinf
    rotor_indices = params.rotor_idxs

    # extract vortex and source states
    Γb, Γw, Γr, Σr = extract_state_variables(states, params)

    # calculate net circulation
    Γ_tilde = calculate_net_circulation(Γr, num_blades)

    # calculate enthalpy jumps
    H_tilde = calculate_enthalpy_jumps(Γr, Ωr, num_blades)

    # calculate meridional velocities at wakes
    wake_velocities = calculate_wake_velocities(
        vx_wb, vr_wb, Γb, vx_ww, vr_ww, Γw, vx_wr, vr_wr, Σr
    )

    # calculate induced velocity at rotors
    Vm, Vθ = calculate_induced_velocities_on_rotors(
        blade_elements, Γr, vx_rb, vr_rb, Γb, vx_rw, vr_rw, Γw
    )

    # update body vortex strengths
    calculate_body_vortex_strengths!(Γb, A_bb, b_fb, vx_bw, vr_bw, Γw, vx_br, vr_br, Σr)

    # update wake vortex strengths
    calculate_wake_vortex_strengths!(
        Γw, wake_panels, wake_velocities, Γ_tilde, H_tilde, rotor_indices
    )

    # update circulation and source strengths
    calculate_gamma_sigma!(Γr, Σr, blade_elements, Vinf, Vm, Vθ)

    return nothing
end

"""
    extract_state_variables(states, params)

Extract circulation and source strengths from the state vector
"""
function extract_state_variables(states, params)

    # Problem Dimensions
    nb = params.num_body_panels                            # number of body panels
    nrotor = params.num_rotors                             # number of rotors
    nx = params.num_wake_x_panels                          # number of streamwise coordinates
    nr = length(params.blade_elements[1].radial_positions) # number of radial coordinates

    # State Variable Indices
    iΓb = 1:nb                                # body vortex strength indices
    iΓw = (iΓb[end] + 1):(iΓb[end] + nx * nr)     # wake vortex strength indices
    iΓr = (iΓw[end] + 1):(iΓw[end] + nr * nrotor) # rotor circulation strength indices
    iΣr = (iΓr[end] + 1):(iΓr[end] + nr * nrotor) # rotor source strength indices

    # Extract State variables
    Γb = view(states, iΓb)                        # body vortex strengths
    Γw = reshape(view(states, iΓw), (nx, nr))     # wake vortex strengths
    Γr = reshape(view(states, iΓr), (nr, nrotor)) # rotor circulation strengths
    Σr = reshape(view(states, iΣr), (nr, nrotor)) # rotor circulation strengths

    return Γb, Γw, Γr, Σr
end
