function solver(
    r_fun!,
    states;
    # Solver Options
    convergence_tolerance=1e-10,
    max_iterations=500,
    relaxation_parameters=(; nrf=0.4, bt1=0.2, bt2=0.6, pf1=0.4, pf2=0.5),
    state_ids=[[1:length(states)]],
)

    # initialize convergence flag
    converged = [false]

    # initialize iterator
    iter = [0]

    # initialize residual
    r_current = similar(states) .= 0
    r_previous = similar(states) .= 0

    # iterate until converged or maximum allowed iterations
    while !converged && iter[] < max_iterations

        # Calculate Gamr Residuals
        r_fun.r_Gamr!(r_current[state_ids[1]], current_states[state_ids[1]])
        # Update Gamr
        relax_Gamr!(
            current_states[state_ids[1]],
            r_current[state_ids[1]],
            r_previous[state_ids[1]],
            relaxation_parameters,
        )

        # Calculate gamw Residuals
        r_fun.r_gamw!(r_current[state_ids[2]], current_states[state_ids[2]])
        # Update gamw
        relax_gamw!(
            current_states[state_ids[2]],
            r_current[state_ids[2]],
            r_previous[state_ids[2]],
            relaxation_parameters,
        )

        # Calculate sigr Residuals
        r_fun.r_sigr!(r_current[state_ids[3]], current_states[state_ids[3]])
        # Update sigr
        relax_sigr!(
            current_states[state_ids[3]],
            r_current[state_ids[3]],
            r_previous[state_ids[3]],
            relaxation_parameters,
        )

        # Check if residuals are converged
        converged[] = maximum(r) <= convergence_tolerance

        # increment iterator
        iter[] += 1

        # Copy over residual
        r_previous .= r_current
    end

    return (;
        y=states, converged=converged[1], total_iterations=iter[1], residual=r_current
    )
end

function relax_Gamr!(Gamr, r_sub_current, r_sub_previous, relaxation_parameters)
    # update states in place
    return Gamr
end

function relax_gamw!(gamw, r_sub_current, r_sub_previous, relaxation_parameters)
    # update states in place
    return gamw
end

function relax_sigr!(sigr, r_sub_current, r_sub_previous, relaxation_parameters)
    # update states in place
    return sigr
end

function residual!(r, current_states, inputs, constants)

    # Split States
    Gamr, gamw, sigr = split_states(current_states, constants.state_ids)

    # Compute all residuals
    residual_Gamr!(r[constants.state_ids[1]], Gamr, inputs, constants)
    residual_gamw!(r[constants.state_ids[2]], gamw, inputs, constants)
    residual_sigr!(r[constants.state_ids[3]], sigr, inputs, constants)

    return r
end

function residual_Gamr(r, current_Gamr, inputs, constants)
    estimated_Gamr = estimate_Gamr(current_Gamr, inputs, constants)
    @. r = estimated_Gamr - current_Gamr
    return r
end

function residual_gamw(r, current_gamw, inputs, constants)
    estimated_gamw = estimate_gamw(current_gamw, inputs, constants)
    @. r = estimated_gamw - current_gamw
    return r
end

function residual_sigr(r, current_sigr, inputs, constants)
    estimated_sigr = estimate_sigr(current_sigr, inputs, constants)
    @. r = estimated_sigr - current_sigr
    return r
end

function estimate_Gamr(current_Gamr, inputs, constants)
    return estimated_Gamr
end

function estimate_gamw(current_gamw, inputs, constants)
    return estimated_gamw
end

function estimate_sigr(current_sigr, inputs, constants)
    return estimated_sigr
end

function solve(inputs, parameters)

    # Grab parts of the complete residual
    residual_wrapper_Gamr(r, states) = residual_Gamr!(r, states, inputs, parameters)
    residual_wrapper_gamw(r, states) = residual_Gamr!(r, states, inputs, parameters)
    residual_wrapper_sigr(r, states) = residual_Gamr!(r, states, inputs, parameters)

    converged_states = solver(
        (;
            r_Gamr=residual_wrapper_Gamr,
            r_gamw=residual_wrapper_gamw,
            r_sigr=residual_wrapper_sigr,
        ),
        states_initial_guess;
        convergence_tolerance=parameters.solver_options.convergence_tolerance,
        max_iterations=parameters.solver_options.max_iterations,
        relaxation_parameters=parameters.solver_options.relaxation_parameters,
        state_ids=parameters.solver_options.state_ids,
    )

    return converged_states
end

function program(inputs, constants)
    converged_states = implicit(solve, residual!, inputs, constants)
    return converged_states
end
