"""
    mod_COR_solver(
        r_fun!,
        states,
        B,
        state_dims;
        convergence_tolerance=1e-10,
        iteration_limit=500,
        relaxation_parameters=(; nrf=0.4, bt1=0.2, bt2=0.6, pf1=0.4, pf2=0.5, btw=0.6, pfw=1.2),
        verbose=false,
    )

Modified DFDC-style Concurrent Successive Over-Relaxation (CSOR) solver that updates all state variables before relaxing circulation states (`Gamr` and `gamw`).

# Arguments
- `r_fun!::Function` : Residual function to compute the current residual vector given the current state.
- `states::AbstractVector` : Initial guess and current state vector of the solver.
- `B::AbstractVector` : Number of blades on each rotor, used for circulation relaxation scaling.
- `state_dims::NamedTuple` : Dimensions and indexing information for the state variables, used in update and relaxation routines.

# Keyword Arguments
- `convergence_tolerance::Real=1e-10` : Absolute tolerance for convergence of the residual norm.
- `iteration_limit::Int=500` : Maximum number of iterations before stopping.
- `relaxation_parameters::NamedTuple` : Named tuple containing relaxation parameters:
    - `nrf` : nominal relaxation factor
    - `bt1`, `bt2` : backtrack factors for relaxation
    - `pf1`, `pf2` : press forward factors for relaxation
    - `btw`, `pfw` : relaxation factors specific to wake circulation updates
- `verbose::Bool=false` : If true, prints iteration progress and residual norms.

# Returns
A named tuple with fields:
- `y` : The final converged state vector (modified in place).
- `converged` : Boolean indicating whether convergence was achieved.
- `total_iterations` : Number of iterations performed.
- `residual` : Final maximum absolute residual norm.
"""
function mod_COR_solver(
    r_fun!,
    states,
    B,
    state_dims;
    convergence_tolerance=1e-10,
    iteration_limit=500,
    relaxation_parameters=(; nrf=0.4, bt1=0.2, bt2=0.6, pf1=0.4, pf2=0.5, btw=0.6, pfw=1.2),
    verbose=false,
)

    # initialize convergence flag
    converged = [false]

    # initialize iterator
    iter = [0]

    # initialize residual
    r_current = similar(states) .= 0
    r_previous = similar(states) .= 0

    # iterate until converged or maximum allowed iterations
    while !converged[] && iter[] < iteration_limit

        # Calculate Residuals
        # Note: you can pass in intermediate computation caches to the residual via the wrapper that takes in inputs and constants
        r_fun!(r_current, states)

        # Update States
        # Note: you can't pass intermediate computation caches here, unless you allow the solver to take in its own cache
        update_states!(states, r_current, r_previous, B, relaxation_parameters, state_dims)

        # Check if residuals are converged
        converged[] = maximum(smooth_abs(r_current)) <= convergence_tolerance
        if verbose
            println("Iteration: $(iter[])")
            println("max r: ", maximum(smooth_abs(r_current)))
            println("Converged? $(converged[])")
        end

        # increment iterator
        iter[] += 1

        # Copy over residual
        r_previous .= r_current
    end

    return (;
            y=states, converged=converged[1], total_iterations=iter[1], residual=maximum(smooth_abs(r_current))
    )
end

"""
    relax_Gamr_mod!(
        Gamr,
        r_Gamr_current,
        r_Gamr_previous,
        B;
        nrf=0.4,
        bt1=0.2,
        bt2=0.6,
        pf1=0.4,
        pf2=0.5,
    )

Apply relaxed update step to rotor circulation strengths `Gamr` in place.

# Arguments
- `Gamr::AbstractMatrix{<:Real}` : Rotor circulation strengths array (rows = blade elements, columns = rotors), updated in place.
- `r_Gamr_current::AbstractMatrix{<:Real}` : Current iteration circulation differences (Gamr estimate - Gamr current).
- `r_Gamr_previous::AbstractMatrix{<:Real}` : Previous iteration circulation differences, updated in place.
- `B::AbstractVector{<:Real}` : Number of blades on each rotor.

# Keyword Arguments
- `nrf::Real=0.4` : Nominal relaxation factor.
- `bt1::Real=0.2` : Backtrack factor for sign changes on circulation.
- `bt2::Real=0.6` : Secondary backtrack factor applied per blade element.
- `pf1::Real=0.4` : Press forward factor for circulation updates.
- `pf2::Real=0.5` : Secondary press forward factor applied per blade element.

# Returns
- `Gamr` : Updated `Gamr` array (in place).
"""
function relax_Gamr_mod!(
    Gamr, B, r_Gamr_current, r_Gamr_previous; nrf=0.4, bt1=0.2, bt2=0.6, pf1=0.4, pf2=0.5
)

    # initilize
    TF = eltype(Gamr)
    bladeomega = nrf .* ones(TF, size(Gamr, 2))
    omega = nrf .* ones(TF, size(Gamr, 1))
    deltahat = zeros(TF, size(Gamr, 1))
    minBGamr = zeros(TF, size(Gamr, 2))
    maxBGamr = zeros(TF, size(Gamr, 2))

    for (i, (BG, b, delta_prev, delta)) in
        enumerate(zip(eachcol(Gamr), B, eachcol(r_Gamr_previous), eachcol(r_Gamr_current)))

        # multiply by B to get BGamr values
        BG .*= b
        delta .*= b
        delta_prev .*= b

        # - Set the normalization value based on the maximum magnitude value of B*Gamr

        # find max magnitude
        maxBGamr[i], maxi = findmax(BG)
        minBGamr[i], mini = findmin(BG)

        # make sure we don't have any weird jumps
        meang = sum(BG) / length(BG)

        if meang > 0.0 # if mean is positive, make sure maxBGamr[i] is at least 0.1
            maxBGamr[i] = max(maxBGamr[i], 0.1)
        elseif meang < 0.0 # if mean is negative, make sure maxBGamr[i] is at most -0.1
            maxBGamr[i] = min(minBGamr[i], -0.1)
        else # if the average is zero, then set maxBGamr[i] to zero
            maxBGamr[i] = 0.0
        end

        # note: delta = Gamr_estimate .- Gamr_current
        # note: deltahat here is actually 1/deltahat which is the version needed later
        for (j, d) in enumerate(eachrow(deltahat))
            if smooth_abs(delta[j]) < eps()
                d[1] = sign(delta[j]) * sign(maxBGamr[i]) #avoid division by zero
            else
                d[1] = maxBGamr[i] ./ delta[j]
            end
        end

        # get initial relaxation factor
        bladeomega[i], oi = findmin(smooth_abs(deltahat))

        # scale relaxation factor based on if the change and old values are the same sign (back tracking or pressing forward)
        if (nrf / deltahat[oi]) < -bt1
            bladeomega[i] *= bt1
        elseif (nrf / deltahat[oi]) > pf1
            bladeomega[i] *= pf1
        else
            bladeomega[i] = nrf
        end

        # scale blade element relaxation factor
        for (o, d, dp) in zip(eachrow(omega), delta, eachrow(delta_prev))
            # if differences changed sign, use backtrack factor, if sign is same, use press forward factor
            o[1] = bladeomega[i] * (sign(d) != sign(dp[1]) ? bt2 : pf2)

            # save current delta into old one for next iteration
            dp[1] = d
        end

        # relax Gamr for this blade
        BG .+= omega .* delta
        # remove the *b
        BG ./= b
        delta ./= b
        delta_prev ./= b
    end

    return Gamr
end

"""
    relax_gamw_mod!(gamw, r_gamw_current, r_gamw_previous; nrf=0.4, btw=0.6, pfw=1.2)

Apply relaxed step to gamw.

# Arguments
- `gamw::Array{Float}` : Array of rotor circulations (columns = rotors, rows = blade elements), updated in place
- `r_gamw_current::Array{Float}` : Array of current iteration's differences in circulation values
- `r_gamw_previous::Array{Float}` : Array of previous iteration's differences in circulation values, updated in place

# Keyword Arguments
- `nrf::Float=0.4` : nominal relaxation factor
- `bt1::Float=0.2` : backtrack factor 1
- `bt2::Float=0.6` : backtrack factor 2
- `pf1::Float=0.4` : press forward factor 1
- `pf2::Float=0.5` : press forward factor 2

# Returns
- `gamw` : Updated `gamw` array (in place).
"""
function relax_gamw_mod!(gamw, r_gamw_current, r_gamw_previous; nrf=0.4, btw=0.6, pfw=1.2)

    # initilize
    TF = eltype(gamw)
    omega = MVector{1,TF}(nrf)

    for ig in eachindex(gamw)

        # choose relaxation factor based on whether old and new changes are in different or same direction
        omega[] =
            sign(r_gamw_previous[ig] * r_gamw_current[ig]) < 0.0 ? btw * nrf : pfw * nrf

        # update gamw value
        gamw[ig] += r_gamw_current[ig] * omega[]
    end

    return gamw
end

"""
    update_states!(states, r_current, r_previous, B, relaxation_parameters, state_dims)

Update states using DFDC-like relaxation methods.

# Arguments
- `states::Vector{Float}` : current iteration states to update
- `r_current::Vector{Float}` : current iteration residual values
- `r_previous::Vector{Float}` : previous iteration residual values
- `B::Vector{Float}` : number of blades for each rotor
- `relaxation_parameters::NamedTuple` : relaxation parameters
- `state_dims::NamedTuple` : dimensions of the state variables

# Returns
- `states::Vector{Float}` : updates states vector in place.
"""
function update_states!(states, r_current, r_previous, B, relaxation_parameters, state_dims)

    # - Separate out the states and residuals - #
    Gamr, sigr, gamw = extract_state_variables(ModCSORSolverOptions(), states, state_dims)
    r_Gamr_current, r_sigr_current, r_gamw_current = extract_state_variables(
        ModCSORSolverOptions(), r_current, state_dims
    )
    r_Gamr_previous, r_sigr_previous, r_gamw_previous = extract_state_variables(
        ModCSORSolverOptions(), r_previous, state_dims
    )

    # println("before: ", Gamr)
    # - relax Gamr values - #
    relax_Gamr_mod!(
        Gamr,
        B,
        r_Gamr_current,
        r_Gamr_previous;
        nrf=relaxation_parameters.nrf,
        bt1=relaxation_parameters.bt1,
        bt2=relaxation_parameters.bt2,
        pf1=relaxation_parameters.pf1,
        pf2=relaxation_parameters.pf2,
    )
    # println("after: ", Gamr)

    # relax gamw values
    relax_gamw_mod!(
        gamw,
        r_gamw_current,
        r_gamw_previous;
        nrf=relaxation_parameters.nrf,
        btw=relaxation_parameters.btw,
        pfw=relaxation_parameters.pfw,
    )

    # - relax sigr values using same method as Gamr - #
    relax_Gamr_mod!(
        sigr,
        B,
        r_sigr_current,
        r_sigr_previous;
        nrf=relaxation_parameters.nrf,
        bt1=relaxation_parameters.bt1,
        bt2=relaxation_parameters.bt2,
        pf1=relaxation_parameters.pf1,
        pf2=relaxation_parameters.pf2,
    )

    return states
end
