# """
#     analyze_propulsor(x, fx=x->x; maximum_linesearch_step_size=1e-8, iteration_limit=100)

# Version of `analyze_propeller` designed for sensitivity analysis.  `x` is an input vector
# and `fx` is a function which returns the standard input arguments to `analyze_propeller`.
# """
# function analyze_propulsor(x, fx=x -> x; maximum_linesearch_step_size=1e-8, iteration_limit=15)

#     # convergence flag
#     converged = [false]

#     # define parameters
#     p = (; fx, maximum_linesearch_step_size, iteration_limit, converged)

#     # compute state variables (updates convergence flag internally)
#     states = ImplicitAD.implicit(solve, residual!, x, p)

#     # TODO: post-processing using the converged state variables

#     # return solution
#     return states, converged[1]
# end

function analyze_propulsor(
    inputs; debug=false, maximum_linesearch_step_size=1e6, iteration_limit=100, ftol=1e-8
)
    initial_states = initialize_states(inputs)
    initials = copy(initial_states)

    # - Define closure that allows for parameters - #
    rwrap(r, states) = residual!(r, states, inputs, [])

    # - Call NLsolve function using AD for Jacobian - #
    #= res is of type NLsolve.SolverResults.
    The zero field contains the "solution" to the non-linear solve.
    The converged() function tells us if the solver converged.
    =#
    res = NLsolve.nlsolve(
        rwrap,
        initial_states;
        autodiff=:forward,
        method=:newton,
        iterations=iteration_limit,
        ftol=ftol,
        show_trace=true,
        linesearch=BackTracking(; maxstep=maximum_linesearch_step_size),
    )

    # converged[1] = res.f_converged

    # return solution
    return res.zero, initials, res.f_converged
end

"""
    analyze_propulsor(duct_coordinates, hub_coordinates, rotor_parameters, freestream;
        maximum_linesearch_step_size=1e-8, iteration_limit=100)

"""
function analyze_propulsor(
    duct_coordinates,
    hub_coordinates,
    paneling_constants,
    rotor_parameters,
    freestream,
    reference_parameters;
    debug=false,
    verbose=false,
    maximum_linesearch_step_size=1e6,
    iteration_limit=100,
    ftol=1e-8,
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
            reference_parameters,
        )

    # convergence flag
    converged = [false]

    # define parameters
    p = (;
        fx, maximum_linesearch_step_size, iteration_limit, ftol, converged, debug, verbose
    )

    # compute various panel and circulation strenghts (updates convergence flag internally)
    strengths, inputs, initials = solve(x, p)
    if debug
        println("NLSolve Complete")
    end

    # post-processing using the converged strengths
    out = post_process(strengths, inputs)

    # return solution
    return out, strengths, inputs, initials, p.converged[1]
end

"""
    solve(x, p)

"""
function solve(x, p)

    # unpack parameters
    (; fx, maximum_linesearch_step_size, iteration_limit, ftol, converged, debug, verbose) =
        p

    # unpack inputs
    (; duct_coordinates, hub_coordinates, paneling_constants, rotor_parameters, freestream, reference_parameters) = fx(
        x
    )

    # initialize various inputs used in analysis
    # repanels bodies and rotors, generates wake "grid", precomputes influence matrices, etc.
    inputs = precomputed_inputs(
        duct_coordinates,
        hub_coordinates,
        paneling_constants,
        rotor_parameters,
        freestream,
        reference_parameters;
        debug=debug,
    )

    # calculate initial guess for state variables
    initial_states = initialize_states(inputs)
    initials = copy(initial_states)

    # - Define closure that allows for parameters - #
    rwrap(r, states) = residual!(r, states, inputs, p)

    # - Call NLsolve function using AD for Jacobian - #
    #= res is of type NLsolve.SolverResults.
    The zero field contains the "solution" to the non-linear solve.
    The converged() function tells us if the solver converged.
    =#
    res = NLsolve.nlsolve(
        rwrap,
        initial_states;
        autodiff=:forward,
        method=:newton,
        iterations=iteration_limit,
        show_trace=verbose,
        linesearch=BackTracking(; maxstep=maximum_linesearch_step_size),
        ftol=ftol,
    )

    # save convergence information
    # converged[1] = NLsolve.converged(res)
    converged[1] = res.f_converged

    # return solution
    return res.zero, inputs, initials
end

"""
    residual!(res, states, inputs, p)

Updates the state variables.
"""
function residual!(res, states, inputs, p)
    updated_states = copy(states)

    update_strengths!(updated_states, inputs, p)

    # Update Residual
    @. res = updated_states - states
    # TODO: need to add the pressure residual here,
    # TODO: this adds one more equation than there are states.
    # TODO; remove the first state residual (associated with the inner duct TE panel strength) and replace with the pressure coefficient residual.
    # @. res[2:end] = updated_states[2:end] - states[2:end]
    # res[1] = cp_residual(states, inputs)

    return nothing
end

function residual(states, inputs, p)
    res = Inf .* ones(eltype(states), length(states) + 1)

    residual!(res, states, inputs, p)

    return res
end

function update_strengths!(states, inputs, p)

    # - Extract states - #
    mub, gamw, Gamr, sigr = extract_state_variables(states, inputs)

    ### --- Get Velocities Before Updating States --- ###
    # - Get velocities at rotor planes - #
    _, _, _, _, Wtheta_rotor, Wm_rotor, Wmag_rotor = calculate_rotor_velocities(
        Gamr, gamw, sigr, mub, inputs
    )

    # - Get velocities on wake panels - #
    Wm_wake = calculate_wake_velocities(gamw, sigr, mub, inputs)

    # - Generate raw RHS, viz. velocities on body, (before updating state dependencies) - #
    RHS = update_RHS(inputs.b_bf, inputs.A_bw, gamw, inputs.A_br, sigr)

    # - Calculate body vortex strengths (before updating state dependencies) - #
    # solve_body_strengths!(
    #     # mub, inputs.A_bb, RHS, inputs.prescribedpanels, inputs.body_doublet_panels.nbodies
    #     mub, inputs.mured, inputs.A_bb, RHS,
    #     inputs.LHSlsq, inputs.LHSlsqlu, inputs.RHSlsq, inputs.tLHSred, inputs.b_bf0,
    #     inputs.prescribedpanels
    # )

    strengths = solve_body_strengths(
        inputs.A_bb,
        RHS,
        inputs.LHSlsq,
        inputs.LHSlsqlu,
        inputs.tLHSred,
        inputs.b_bf0,
        inputs.prescribedpanels,
    )
    mub .= strengths

    # - Calculate wake vortex strengths (before updating state dependencies) - #
    calculate_wake_vortex_strengths!(gamw, Gamr, Wm_wake, inputs)

    # - Update rotor circulation and source panel strengths - #
    calculate_gamma_sigma!(
        Gamr,
        sigr,
        inputs.blade_elements,
        Wm_rotor,
        Wtheta_rotor,
        Wmag_rotor,
        inputs.freestream;
    )

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
    nx = inputs.num_wake_x_panels                   # number of wake panels per sheet

    # State Variable Indices
    iμb = 1:nb                                      # body vortex strength indices
    iΓw = (iμb[end] + 1):(iμb[end] + nw * nx)       # wake vortex strength indices
    iΓr = (iΓw[end] + 1):(iΓw[end] + nr * nrotor)   # rotor circulation strength indices
    iΣr = (iΓr[end] + 1):(iΓr[end] + nr * nrotor)   # rotor source strength indices

    # Extract State variables
    mub = view(states, iμb)                        # body vortex strengths
    gamw = view(states, iΓw)     # wake vortex strengths
    Gamr = reshape(view(states, iΓr), (nr, nrotor)) # rotor circulation strengths
    sigr = reshape(view(states, iΣr), (nr, nrotor)) # rotor circulation strengths

    return mub, gamw, Gamr, sigr
end

################################################################################
#TODO: move to initialize.jl
# PREPROCESSING OF LINEAR SYSTEM
################################################################################
"""
    calc_Alu!(Apivot::AbstractMatrix, A::AbstractMatrix) -> Alu

Returns the LU decomposition of `A` using `Apivot` as storage memory to pivot
leaving `A` unchanged.
"""
function calc_Alu!(Apivot, A::AbstractMatrix{T}) where {T}

    # Prepare pivot array
    calc_Avalue!(Apivot, A)

    # LU decomposition
    Alu = lu!(Apivot)

    return Alu
end

"""
    calc_Alu!(A::AbstractMatrix) -> Alu

Returns the LU decomposition of `A`. If `A` does not carry Dual nor TrackedReal
numbers, computation will be done in-place using `A`; hence `A` should not be
reused for multiple solves or for implicit differentiation (use `calc_Alu(A)`
and `calc_Alu!(Apivot, A)` instead).
"""
function calc_Alu!(A::AbstractMatrix{T}) where {T}

    # Allocate memory for pivot
    if T <: ImplicitAD.ForwardDiff.Dual || T <: ImplicitAD.ReverseDiff.TrackedReal  # Automatic differentiation case
        Tprimal = T.parameters[T <: ImplicitAD.ForwardDiff.Dual ? 2 : 1]
        Apivot = zeros(Tprimal, size(A))

        # LU decomposition
        Alu = calc_Alu!(Apivot, A)

    else
        # LU decomposition
        Alu = LA.lu!(A)
    end

    return Alu
end

"""
    calc_Alu(A::AbstractMatrix) -> Alu

Returns the LU decomposition of `A`.
"""
function calc_Alu(A::AbstractMatrix{T}) where {T}

    # Allocate memory for pivot
    if T <: ImplicitAD.ForwardDiff.Dual || T <: ImplicitAD.ReverseDiff.TrackedReal  # Automatic differentiation case
        Tprimal = T.parameters[T <: ImplicitAD.ForwardDiff.Dual ? 2 : 1]
        Apivot = zeros(Tprimal, size(A))

    else
        Apivot = zeros(T, size(A))
    end

    # LU decomposition
    return calc_Alu!(Apivot, A)
end

"""
    calc_Avalue!(Avalue::AbstractMatrix, A::AbstractMatrix)

Copy the primal values of `A` into `Avalue`.
"""
function calc_Avalue!(Avalue, A::AbstractMatrix{T}) where {T}
    if T <: ImplicitAD.ForwardDiff.Dual || T <: ImplicitAD.ReverseDiff.TrackedReal  # Automatic differentiation case

        # Extract primal values of A
        value = if T <: ImplicitAD.ForwardDiff.Dual
            ImplicitAD.ForwardDiff.value
        else
            ImplicitAD.ForwardDiff.value
        end
        map!(value, Avalue, A)

    else                                # Normal case
        # Deep copy A
        Avalue .= A
    end

    return Avalue
end

"""
    calc_Avalue(A::AbstractMatrix)

Return the primal values of matrix `A`, which is simply `A` if the elements
of `A` are not Dual nor TrackedReal numbers.
"""
function calc_Avalue(A::AbstractMatrix{T}) where {T}
    if T <: ImplicitAD.ForwardDiff.Dual || T <: ImplicitAD.ReverseDiff.TrackedReal  # Automatic differentiation case
        Tprimal = T.parameters[T <: ImplicitAD.ForwardDiff.Dual ? 2 : 1]
        Avalue = zeros(Tprimal, size(A))
        calc_Avalue!(Avalue, A)

        return Avalue
    else                                # Normal case
        return A
    end
end
#### END OF LINEAR-SOLVER PREPROCESSING ########################################