# """
# rotor only solve, calls nominal solve function but sets gamb to nothing
# """
# function solve!(
#     inputs, Gamr, Gamr_est, gamw, gamw_est; nosource=true, maxiter=1e2, verbose=false
# )
#     return solve!(
#         inputs,
#         Gamr,
#         Gamr_est,
#         gamw,
#         gamw_est,
#         nothing;
#         nosource=nosource,
#         maxiter=maxiter,
#         verbose=verbose,
#     )
# end

"""

"""
function solve!(
    inputs,
    Gamr,
    Gamr_est,
    sigr,
    sigr_est,
    gamw,
    gamw_est,
    gamb=nothing;
    nosource=false,
    maxiter=1e2,
    verbose=false,
    nrf=0.5,
    bt1=0.2,
    bt2=0.6,
    pf1=0.4,
    pf2=0.5,
    btw=0.6,
    pfw=1.2,
    f_circ=1e-3, #DFDC values
    # f_sig=1e-3, #tighter tolerances don't seem to help much
    f_dgamw=2e-4, #DFDC Values
    # f_circ=1e-10, #tighter tolerances don't seem to help much
    # f_dgamw=2e-12, # tighter tolderances don't seem to help much
)
    freestream = inputs.freestream

    # initialize convergence criteria
    TF = eltype(Gamr)
    maxBGamr = MVector{1,TF}(0.0)
    maxdeltaBGamr = MVector{1,TF}(0.0)
    maxsigr = MVector{1,TF}(0.0)
    maxdeltasigr = MVector{1,TF}(0.0)
    maxdeltagamw = MVector{1,TF}(0.0)
    conv = MVector{1,Bool}(false)
    iter = 0

    # initialize differences
    deltaG_prev = Gamr_est .- Gamr
    deltaG = similar(deltaG_prev) .= 0.0
    deltaS_prev = sigr_est .- sigr
    deltaS = similar(deltaS_prev) .= 0.0
    deltag_prev = gamw_est .- gamw
    deltag = similar(deltag_prev) .= 0.0

    # zero out sigr if not used
    if nosource
        sigr .= 0.0
    end

    # loop until converged or max iterations are reached
    while !conv[] && iter <= maxiter
        # update iteration number
        iter += 1
        if verbose
            println("Iteration $(iter):")
        end

        # Solve linear system if including
        if !isnothing(gamb)
            # in place solve for gamb
            # TODO: update this function to match new methods
            calculate_body_vortex_strengths!(
                gamb,
                inputs.A_bb,
                inputs.b_bf,
                inputs.kidx,
                inputs.A_bw,
                gamw,
                inputs.A_br,
                sigr,
                inputs.ductwakeidx,
            )

            # Update rotor blade element velocities
            _, _, _, _, Wtheta_rotor, Wm_rotor, Wmag_rotor = calculate_rotor_velocities(
                Gamr, gamw, sigr, gamb, inputs
            )
        else

            # Update rotor blade element velocities
            _, _, _, _, Wtheta_rotor, Wm_rotor, Wmag_rotor = calculate_rotor_velocities(
                Gamr, gamw, sigr, inputs
            )
        end

        # in-place solve for Gamr, updating Gamr_est
        calculate_gamma_sigma!(
            Gamr_est,
            sigr,
            inputs.blade_elements,
            Wm_rotor,
            Wtheta_rotor,
            Wmag_rotor,
            freestream;
            debug=false,
            verbose=false,
        )

        # get difference between estimated Gamr and old Gamr
        @. deltaG = Gamr_est - Gamr

        # relax Gamr values
        relax_Gamr!(
            Gamr,
            deltaG_prev,
            deltaG,
            maxBGamr,
            maxdeltaBGamr,
            inputs.blade_elements.B;
            nrf=nrf,
            bt1=bt1,
            bt2=bt2,
            pf1=pf1,
            pf2=pf2,
        )

        if nosource
            sigr .= 0.0
        else
            # get difference between estimated Gamr and old Gamr
            @. deltaS = sigr_est - sigr

            # relax Gamr values
            relax_sigr!(
                sigr,
                deltaS_prev,
                deltaS,
                maxsigr,
                maxdeltasigr;
                nrf=nrf,
                bt1=bt1,
                bt2=bt2,
                pf1=pf1,
                pf2=pf2,
            )
        end

        # Update Vm_avg in wake using new Gamr and sigma
        if !isnothing(gamb)
            Wm_wake = calculate_wake_velocities(gamw, sigr, gamb, inputs)
        else
            Wm_wake = calculate_wake_velocities(gamw, sigr, inputs)
        end

        # in-place solve for gamw, update gamw_est
        calculate_wake_vortex_strengths!(gamw_est, Gamr, Wm_wake, inputs)

        # get difference beetween estimated gamw and old gamw
        @. deltag = gamw_est - gamw

        # relax gamw values
        relax_gamw!(gamw, deltag_prev, deltag, maxdeltagamw; nrf=nrf, btw=btw, pfw=pfw)

        # converged?
        if nosource
        check_convergence!(
            conv,
            maxBGamr,
            maxdeltaBGamr,
            maxdeltagamw,
            # inputs.reference_parameters.Vref;
            inputs.Vconv[1],
            f_circ=f_circ,
            f_dgamw=f_dgamw,
            verbose=verbose,
        )
        else
        check_convergence!(
            conv,
            maxBGamr,
            maxdeltaBGamr,
            maxsigr,
            maxdeltasigr,
            maxdeltagamw,
            # inputs.reference_parameters.Vref;
            inputs.Vconv[1];
            f_circ=f_circ,
            f_dgamw=f_dgamw,
            verbose=verbose,
        )
    end

        # print iteration information if verbose is true
    end

    inputs.converged[1] = conv[1]

    return Gamr, sigr, gamw
end

"""
- `Gamr::Array{Float}` : Array of rotor circulations (columns = rotors, rows = blade elements), updated in place
- `delta_prev_mat::Array{Float}` : Array of previous iteration's differences in circulation values, updated in place
- `delta_mat::Array{Float}` : Array of current iteration's differences in circulation values
- `B::Vector{Float}` : number of blades on each rotor
- `nrf::Float=0.5` : nominal relaxation factor
- `bt1::Float=0.2` : backtrack factor 1
- `bt2::Float=0.6` : backtrack factor 2
- `pf1::Float=0.4` : press forward factor 1
- `pf2::Float=0.5` : press forward factor 2
"""
function relax_Gamr!(
    Gamr,
    delta_prev_mat,
    delta_mat,
    maxBGamr,
    maxdeltaBGamr,
    B;
    nrf=0.5,
    bt1=0.2,
    bt2=0.6,
    pf1=0.4,
    pf2=0.5,
    test=false,
)

    # initilize
    TF = eltype(Gamr)
    bladeomega = nrf .* ones(TF, size(Gamr, 2))
    omega = nrf .* ones(TF, size(Gamr, 1))
    deltahat = zeros(TF, size(Gamr, 1))

    for (i, (BG, b, delta_prev, delta)) in
        enumerate(zip(eachcol(Gamr), B, eachcol(delta_prev_mat), eachcol(delta_mat)))

        # multiply by B to get BGamr values
        BG.*=b
        delta .*= b

        # - Set the normalization value based on the maximum magnitude value of B*Gamr

        # find max magnitude
        maxBGamr[i], mi = findmax(abs.(BG))

        # maintain sign of original value
        maxBGamr[i] *= sign(BG[mi])

        # make sure we don't have any weird jumps
        meang = sum(BG) / length(BG)
        if meang > 0.0 # if mean is positive, make sure maxBGamr[i] is at least 0.1
            maxBGamr[i] = max(maxBGamr[i], 0.1)
        elseif meang < 0.0 # if mean is negative, make sure maxBGamr[i] is at most -0.1
            maxBGamr[i] = min(maxBGamr[i], -0.1)
        else # if the average is zero, then set maxBGamr[i] to zero
            maxBGamr[i] = 0.0
        end

        # note: delta = Gamr_new .- Gamr
        # note: deltahat here is actually 1/deltahat which is the version needed later
        for (j, d) in enumerate(eachrow(deltahat))
            if abs(delta[j]) < eps()
                d[1] = sign(delta[j])*sign(maxdBGamr[i]) #avoid division by zero
            else
                d[1] =  maxBGamr[i]./delta[j]
            end
        end

        # get initial relaxation factor
        bladeomega[i], oi = findmin(abs.(deltahat))

        # scale relaxation factor based on if the change and old values are the same sign (back tracking or pressing forward)
        if -0.2 > nrf/(deltahat[oi]*maxBGamr[i]) > 0.4
            bladeomega[i] *= sign(deltahat[oi]) < 0.0 ? bt1 : pf1
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

        # save max difference for convergence criteria
        # maxdeltaBGamr[i] = maximum(delta)
        maxdeltaBGamr[i], mdi = findmax(abs.(delta))
        # maintain sign of original value
        maxdeltaBGamr[i] *= sign(delta[mdi])

        # relax Gamr for this blade
        BG .+= omega .* delta
        # remove the *b
        BG ./= b
    end

    if test
        return Gamr, bladeomega, omega
    else
    return Gamr
end
end

"""
- `sigr::Array{Float}` : Array of rotor circulations (columns = rotors, rows = blade elements), updated in place
- `delta_prev_mat::Array{Float}` : Array of previous iteration's differences in circulation values, updated in place
- `delta_mat::Array{Float}` : Array of current iteration's differences in circulation values
- `B::Vector{Float}` : number of blades on each rotor
- `nrf::Float=0.5` : nominal relaxation factor
- `bt1::Float=0.2` : backtrack factor 1
- `bt2::Float=0.6` : backtrack factor 2
- `pf1::Float=0.4` : press forward factor 1
- `pf2::Float=0.5` : press forward factor 2
"""
function relax_sigr!(
    sigr,
    delta_prev_mat,
    delta_mat,
    maxsigr,
    maxdeltasigr;
    nrf=0.5,
    bt1=0.2,
    bt2=0.6,
    pf1=0.4,
    pf2=0.5,
    test=false,
)

    # initilize
    TF = eltype(sigr)
    omega = nrf .* ones(TF, size(sigr, 1))
    deltahat = zeros(TF, size(sigr, 1))
    bladeomega = MVector{1,TF}(0.5)

    for (i, (S, delta_prev, delta)) in
        enumerate(zip(eachcol(sigr), eachcol(delta_prev_mat), eachcol(delta_mat)))
        # - Set the normalization value based on the maximum magnitude value of B*sigr

        # find max magnitude
        maxsigr[i], mi = findmax(abs.(S))

        # maintain sign of original value
        maxsigr[i] *= sign(S[mi])

        # make sure we don't have any weird jumps
        meang = sum(S) / length(S)
        if meang > 0.0 # if mean is positive, make sure maxsigr[i] is at least 0.1
            maxsigr[i] = max(maxsigr[i], 0.1)
        elseif meang < 0.0 # if mean is negative, make sure maxsigr[i] is at most -0.1
            maxsigr[i] = min(maxsigr[i], -0.1)
        else # if the average is zero, then set maxsigr[i] to zero
            maxsigr[i] = 0.0
        end

        # note: delta = sigr_new .- sigr
        for (j, d) in enumerate(eachrow(deltahat))
            if delta[j] < eps()
                d[1] = 0.0
            else
                d[1] = maxsigr[i] ./ delta[j]
            end
        end

        # get initial relaxation factor
        bladeomega[i], oi = findmin(abs.(deltahat))

        # scale relaxation factor based on if the change and old values are the same sign (back tracking or pressing forward)
        bladeomega[i] *= sign(deltahat[oi]) < 0.0 ? bt1 : pf1

        # scale blade element relaxation factor
        for (o, d, dp) in zip(eachrow(omega), delta, eachrow(delta_prev))
            # if differences changed sign, use backtrack factor, if sign is same, use press forward factor
            o[1] = bladeomega[i] * (sign(d) != sign(dp[1]) ? bt2 : pf2)

            # save current delta into old one for next iteration
            dp[1] = d
        end

        # save max relaxation factor for convergence criteria
        maxdeltasigr[i] = maximum(omega)

        # relax sigr for this blade
        S .+= omega .* delta
    end

    if test
    else
    return sigr
end
end

"""
- `gamw::Array{Float}` : Array of rotor circulations (columns = rotors, rows = blade elements), updated in place
- `delta_prev_mat::Array{Float}` : Array of previous iteration's differences in circulation values, updated in place
- `delta_mat::Array{Float}` : Array of current iteration's differences in circulation values
- `B::Vector{Float}` : number of blades on each rotor
- `nrf::Float=0.5` : nominal relaxation factor
- `bt1::Float=0.2` : backtrack factor 1
- `bt2::Float=0.6` : backtrack factor 2
- `pf1::Float=0.4` : press forward factor 1
- `pf2::Float=0.5` : press forward factor 2
"""
function relax_gamw!(gamw, delta_prev, delta, maxdeltagamw; nrf=0.5, btw=0.6, pfw=1.2, test=false)

    # initilize
    TF = eltype(gamw)
    omega = MVector{1,TF}(0.5)

    # update delta gamw max for convergence criteria
    maxdeltagamw[], mi = findmax(abs.(delta))
    maxdeltagamw[] *= sign(delta[mi])

    # use delta_prev as a place holder for relaxation factor criteria
    delta_prev .*= delta

    for ig in eachindex(gamw)

        # choose relaxation factor based on whether old and new changes are in different or same direction
        # note that delta_prev at this point = delta_prev .* delta
        omega[] = sign(delta_prev[ig]) < 0.0 ? btw * nrf : pfw * nrf

        # save delta_prev for next iteration
        delta_prev[ig] = delta[ig] * omega[]

        # update gamw value
        gamw[ig] += delta_prev[ig]
    end

    if test
        return gamw, omega
    else
    return gamw
end
end

"""
"""
function check_convergence!(
    conv,
    maxBGamr,
    maxdeltaBGamr,
    maxsigr,
    maxdeltasigr,
    maxdeltagamw,
    Vref;
    f_circ=1e-3,
    f_dgamw=2e-4,
    verbose=false,
)

    # find max ratio among blades and use that for convergence
    _, idG = findmax(maxdeltaBGamr ./ maxBGamr)
    _, idS = findmax(maxdeltasigr ./ maxsigr)

    # set convergence flag
    conv[] =
        abs(maxdeltaBGamr[idG]) < f_circ * abs(maxBGamr[idG]) &&
        abs(maxdeltasigr[idS]) < f_circ * abs(maxsigr[idS]) &&
        maxdeltagamw[] < f_dgamw * Vref # abs already taken care of
        # abs(maxdeltagamw[]) < f_dgamw * Vref

    # # set convergence flag, note: this is how dfdc does it, without regard to which rotor for the Gamr values
    # conv[] =
    # max(abs.(maxdeltaBGamr)...) < f_circ * max(abs.(maxBGamr)...) &&
    #     maxdeltagamw[] < f_dgamw * Vref

    if verbose
        # ff = Printf.Format("%1.6f    %1.6f")
        # fs = Printf.Format("%-16s    %-16s")
        @printf "\t%-16s      %-16s" "maxdBGamr"  "f*maxBGamr"
        println()
        @printf "\t%1.16f    %1.16f" abs(maxdeltaBGamr[idG])  f_circ * abs(maxBGamr[idG])
        println()
        @printf "\t%-16s      %-16s" "maxdsigr"  "f*maxsigr"
        println()
        @printf "\t%1.16f    %1.16f" abs(maxdeltasigr[idS])  f_circ * abs(maxsigr[idS])
        println()
        @printf "\t%-16s      %-16s" "maxdgamw"  "f*Vconv"
        println()
        @printf "\t%1.16f    %1.16f" abs(maxdeltagamw[])  f_dgamw * Vref
        println()
        println()
    end

    return conv
end

function check_convergence!(
    conv,
    maxBGamr,
    maxdeltaBGamr,
    maxdeltagamw,
    Vref;
    f_circ=1e-3,
    f_dgamw=2e-4,
    verbose=false,
)

    # find max ratio among blades and use that for convergence
    _, idG = findmax(maxdeltaBGamr ./ maxBGamr)

    # set convergence flag
    conv[] =
        abs(maxdeltaBGamr[idG]) < f_circ * abs(maxBGamr[idG]) &&
        maxdeltagamw[] < f_dgamw * Vref # abs already taken care of
        # abs(maxdeltagamw[]) < f_dgamw * Vref

    # # set convergence flag, note: this is how dfdc does it, without regard to which rotor for the Gamr values
    # conv[] =
    # max(abs.(maxdeltaBGamr)...) < f_circ * max(abs.(maxBGamr)...) &&
    #     maxdeltagamw[] < f_dgamw * Vref


    return conv
end

# OLD STUFF
#
#
#
#
#
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
    calculate_Gamr_sigma!(
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
