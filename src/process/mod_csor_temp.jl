"""
    solver(
        r_fun!,
        states;
        convergence_tolerance=1e-10,
        max_iterations=500,
        relaxation_parameters=(; nrf=0.4, bt1=0.2, bt2=0.6, pf1=0.4, pf2=0.5),
        state_ids=[[1:length(states)]],
    )

Iterative solver for coupled state residuals with relaxation applied separately to subsets of the state vector.

# Arguments
- `r_fun!::NamedTuple` : A container of residual functions, expected to have callable fields:
    - `r_Gamr!`(residual_subset, state_subset)
    - `r_gamw!`(residual_subset, state_subset)
    - `r_sigr!`(residual_subset, state_subset)
- `states::AbstractVector{<:Number}` : The combined vector of state variables, updated in-place.

# Keyword Arguments
- `convergence_tolerance::Float64=1e-10` : Absolute tolerance for convergence based on maximum residual magnitude.
- `max_iterations::Int=500` : Maximum number of iterations before solver stops.
- `relaxation_parameters::NamedTuple` : Parameters controlling relaxation steps, keys expected include:
    - `nrf` (nominal relaxation factor)
    - `bt1`, `bt2` (backtrack factors)
    - `pf1`, `pf2` (press forward factors)
- `state_ids::Vector{Vector{Int}}` : Vector of index vectors, each selecting the subset of `states` corresponding to:
    - Gamr indices (e.g. `state_ids[1]`)
    - gamw indices  (e.g. `state_ids[2]`)
    - sigr indices  (e.g. `state_ids[3]`)

# Returns
- A NamedTuple with the following fields:
    - `y` : The updated state vector after convergence or termination.
    - `converged::Bool` : Whether the solver converged within the iteration limit.
    - `total_iterations::Int` : Number of iterations performed.
    - `residual::Vector` : Residual vector at the final iteration.
"""
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

"""
    relax_Gamr!(Gamr, r_sub_current, r_sub_previous, relaxation_parameters)

Relax the circulation state vector `Gamr` by applying a relaxation update based on the current and previous residuals.

# Arguments
- `Gamr::AbstractArray` : Array representing the current circulation strengths, updated in-place.
- `r_sub_current::AbstractArray` : Residual vector for the current iteration.
- `r_sub_previous::AbstractArray` : Residual vector from the previous iteration, used to inform relaxation.
- `relaxation_parameters::NamedTuple` : Parameters controlling the relaxation scheme (e.g., relaxation factors).

# Returns
- `Gamr` (updated in-place).
"""
function relax_Gamr!(Gamr, r_sub_current, r_sub_previous, relaxation_parameters)
    # update states in place
    return Gamr
end

"""
    relax_gamw!(gamw, r_sub_current, r_sub_previous, relaxation_parameters)

Relax the wake vortex strengths state vector `gamw` by applying a relaxation update based on residuals.

# Arguments
- `gamw::AbstractArray` : Array representing current wake vortex strengths, updated in-place.
- `r_sub_current::AbstractArray` : Residual vector for the current iteration.
- `r_sub_previous::AbstractArray` : Residual vector from the previous iteration.
- `relaxation_parameters::NamedTuple` : Parameters controlling relaxation.

# Returns
- `gamw` (updated in-place).
"""
function relax_gamw!(gamw, r_sub_current, r_sub_previous, relaxation_parameters)
    # update states in place
    return gamw
end

"""
    relax_sigr!(sigr, r_sub_current, r_sub_previous, relaxation_parameters)

Relax the rotor source strengths state vector `sigr` by applying a relaxation update based on residuals.

# Arguments
- `sigr::AbstractArray` : Array representing rotor source strengths, updated in-place.
- `r_sub_current::AbstractArray` : Residual vector for current iteration.
- `r_sub_previous::AbstractArray` : Residual vector from previous iteration.
- `relaxation_parameters::NamedTuple` : Parameters controlling relaxation.

# Returns
- `sigr` (updated in-place).
"""
function relax_sigr!(sigr, r_sub_current, r_sub_previous, relaxation_parameters)
    # update states in place
    return sigr
end

"""
    residual!(r, current_states, inputs, constants)

Compute the combined residual vector `r` for all state variables by splitting the full state vector and computing residuals for each subset.

# Arguments
- `r::AbstractVector` : Residual vector to be filled in-place.
- `current_states::AbstractVector` : Current full state vector.
- `inputs` : Inputs to the residual calculations (e.g., flow conditions).
- `constants` : Constants and metadata (e.g., index partitions) used in computations.

# Returns
- The updated residual vector `r`.
"""
function residual!(r, current_states, inputs, constants)

    # Split States
    Gamr, gamw, sigr = split_states(current_states, constants.state_ids)

    # Compute all residuals
    residual_Gamr!(r[constants.state_ids[1]], Gamr, inputs, constants)
    residual_gamw!(r[constants.state_ids[2]], gamw, inputs, constants)
    residual_sigr!(r[constants.state_ids[3]], sigr, inputs, constants)

    return r
end

"""
    residual_Gamr(r, current_Gamr, inputs, constants)

Compute the residual for the circulation strengths `Gamr` as the difference between estimated and current values.

# Arguments
- `r::AbstractVector` : Residual vector for Gamr, updated in-place.
- `current_Gamr::AbstractVector` : Current Gamr state vector.
- `inputs` : Inputs required for estimation.
- `constants` : Constants required for estimation.

# Returns
- The updated residual vector `r`.
"""
function residual_Gamr(r, current_Gamr, inputs, constants)
    estimated_Gamr = estimate_Gamr(current_Gamr, inputs, constants)
    @. r = estimated_Gamr - current_Gamr
    return r
end

"""
    residual_gamw(r, current_gamw, inputs, constants)

Compute the residual for the wake vortex strengths `gamw`.

# Arguments
- `r::AbstractVector` : Residual vector for gamw, updated in-place.
- `current_gamw::AbstractVector` : Current gamw state vector.
- `inputs` : Inputs required for estimation.
- `constants` : Constants required for estimation.

# Returns
- The updated residual vector `r`.
"""
function residual_gamw(r, current_gamw, inputs, constants)
    estimated_gamw = estimate_gamw(current_gamw, inputs, constants)
    @. r = estimated_gamw - current_gamw
    return r
end

"""
    residual_sigr(r, current_sigr, inputs, constants)

Compute the residual for rotor source strengths `sigr`.

# Arguments
- `r::AbstractVector` : Residual vector for sigr, updated in-place.
- `current_sigr::AbstractVector` : Current sigr state vector.
- `inputs` : Inputs required for estimation.
- `constants` : Constants required for estimation.

# Returns
- The updated residual vector `r`.
"""
function residual_sigr(r, current_sigr, inputs, constants)
    estimated_sigr = estimate_sigr(current_sigr, inputs, constants)
    @. r = estimated_sigr - current_sigr
    return r
end

"""
    estimate_Gamr(current_Gamr, inputs, constants)

Estimate updated circulation strengths `Gamr` based on current state and inputs.

# Arguments
- `current_Gamr::AbstractVector` : Current Gamr state.
- `inputs` : Inputs needed for estimation.
- `constants` : Constants needed for estimation.

# Returns
- Estimated Gamr state vector.
"""
function estimate_Gamr(current_Gamr, inputs, constants)
    return estimated_Gamr
end

"""
    estimate_gamw(current_gamw, inputs, constants)

Estimate updated wake vortex strengths `gamw`.

# Arguments
- `current_gamw::AbstractVector` : Current gamw state.
- `inputs` : Inputs needed for estimation.
- `constants` : Constants needed for estimation.

# Returns
- Estimated gamw state vector.
"""
function estimate_gamw(current_gamw, inputs, constants)
    return estimated_gamw
end

"""
    estimate_sigr(current_sigr, inputs, constants)

Estimate updated rotor source strengths `sigr`.

# Arguments
- `current_sigr::AbstractVector` : Current sigr state.
- `inputs` : Inputs needed for estimation.
- `constants` : Constants needed for estimation.

# Returns
- Estimated sigr state vector.
"""
function estimate_sigr(current_sigr, inputs, constants)
    return estimated_sigr
end

"""
    solve(inputs, parameters)

Solve for state variables by iteratively calling the solver with residual functions.

# Arguments
- `inputs` : Inputs required for residual and estimate calculations.
- `parameters` : Parameters including initial guesses, solver options, and state indexing.

# Returns
- The converged states and solver metadata as returned by the solver.
"""
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

"""
    program(inputs, constants)

Wraps the full implicit solve process, calling `solve` and `residual!` with the provided inputs and constants.

# Arguments
- `inputs` : Inputs to the solve routine.
- `constants` : Constants required for the solver.

# Returns
- The final converged state output.
"""
function program(inputs, constants)
    converged_states = implicit(solve, residual!, inputs, constants)
    return converged_states
end
