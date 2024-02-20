######################################################################
#                                                                    #
#                        FUNCTION SET HEADER                         #
#                                                                    #
######################################################################

# """
# notional objective function for an optimization
# """
# function objective!(constraints, design_vars, optimization_parameters; cache=make_a_cache)
#     parameters = taylors_fancy_package_function(optimization_parameters, design_vars)
#     fan_params = f(parameters)
#     aero_outs = analyze_propulsor(fan_params; cache=cache)
#     objective = fun(aero_outs)
#     func!(constraints, aero_outs)
#     return objective
# end

"""
"""
function define_constants(;
    # - General Parameters - #
    verbose=false,
    # - NLSolve Parameters - #
    # nlsolve parameters
    nlsolve_method=:newton,
    nlsolve_autodiff=:forward,
    nlsolve_ftol=1e-8, #1e-8 is nlsolve default
    nlsolve_iteration_limit=1e3, #1000 is nlsolve default
    show_trace=false,
    # line search parameters
    dolinesearch=true,
    maximum_linesearch_step_size=1e6,
    # - Post-processing Options - #
    write_outputs=false,
    outfile="outputs.jl",
    checkoutfileexists=false,
    tuple_name="outs",
)
    return (;
        # - General Parameters - #
        verbose,
        # - NLSolve Parameters - #
        # nlsolve parameters
        nlsolve_method,
        nlsolve_autodiff,
        nlsolve_ftol,
        nlsolve_iteration_limit,
        show_trace,
        # line search parameters
        nlsolve_linesearch=if dolinesearch
            LineSearches.BackTracking(; maxstep=maximum_linesearch_step_size)
        else
            LineSearches.Static()
        end,
    )
end

"""
"""
function analyze_propulsor(fan_params, constants; cache=generate_cache(fan_params))

    # - combine cache and constants - #
    const_cache = (; constants..., cache, fx=x -> (fan_params))

    # - Vectorize fan_params - #
    inputs = vectorize!(fan_params, const_cache)

    # - Initialize Aero and Converge Gamr and gamw - #
    Gamr_gamw_sigr = implicit(solve_iad!, iad_res_fun!, inputs, const_cache)

    # - Post-Process - #
    # re-populate cache in this scope so that the derivatives between the fan_params and the post-processing are visible.
    # also for debugging to make sure you're not magically passing things into the cache inside the solve and then using them outside without tracking it.
    # note: will be computing gamb inside the post-processing function anyway so don't need to call the whole residual probably, but may still need to depending on what you want available (e.g. if we put rotor velocities in the cache, then it would make sense to run the whole residual again).
    populate_cache!(cache, fan_params)
    # iad_res_fun!(similar(Gamr_gamw_sigr) .= 0, Gamr_gamw_sigr, inputs, const_cache)

    # de-vectorize Gamr and gamw
    Gamr, gamw, sigr = extract_state_vars(Gamr_gamw_sigr, const_cache.cache)

    # do the rest of the post-processing
    outs = post_process_iad(Gamr, gamw, sigr, const_cache)

    return outs
end

#---------------------------------#
#              SOLVE              #
#---------------------------------#

"""
"""
function solve_iad!(r, differentiated_variables, const_cache)

    ##### ----- Do everything in here ----- #####
    # Note: put all the precomputation in here, but not in the nlsolve residual.
    # Note: it will have to be done twice, but oh well.
    # - Extract constants - #
    (;
        fx,
        nlsolve_method,
        nlsolve_autodiff,
        nlsolve_linesearch,
        nlsolve_ftol,
        nlsolve_iteration_limit,
        verbose,
        cache,
    ) = const_cache

    # - Get fan_params back out using fx(x) in constants - #
    (; fan_params) = fx(differentiated_variables)

    # - Do precomputations - #
    populate_cache!(cache, fan_params)

    # - Initialize Aero - #
    # TODO: put initial guess generation here.
    # Note: you can solve the isolated body system and use that for the initial guess here because it's just an initial guess, it doesn't matter what it is relative to the differentiated variables.
    Gamr, gamw, sigr = initialize_aero_iad(cache)

    # - Wrap residual - #
    # TODO: the differentiated_variables don't need to get passed in as anything do they? I can just pass in an empty vector here.
    function rwrap!(r, state_variables)
        # return nls_res_fun!(r, state_variables, differentiated_variables, constants)
        return nls_res_fun!(r, state_variables, cache)
    end

    res = NLsolve.nlsolve(
        rwrap!,
        [Gamr; gamw; sigr]; # initial states guess
        method=nlsolve_method, #:newton
        autodiff=nlsolve_autodiff, #:forward
        linesearch=nlsolve_linesearch, #BackTracking(; maxstep=maximum_linesearch_step_size), #inside constants
        linsolve=(x, A, b) -> x .= ImplicitAD.implicit_linear(A, b), #used in newton method, unused otherwise so fine to keep it here.
        ftol=nlsolve_ftol,
        iterations=nlsolve_iteration_limit, #inside constants
        show_trace=verbose, #inside constants
    )

    r .= res.zero

    return res.zero
end

#---------------------------------#
#             RESIDUALS           #
#---------------------------------#

"""
This is the residual that gets passed into ImplicitAD, and includes all the precomputation stuff.
basically just copy and paste the solve function but replace the nlsolve with a simple call the the nlsolve residual function.
"""
function iad_res_fun!(r, state_variables, differentiated_variables, const_cache)

    # Do all the initialization and setup stuff
    # call the nlsolve residual function
    # return the residual and various outputs if you want to return them

    # - Extract constants - #
    (;
        fx,
        nlsolve_method,
        nlsolve_autodiff,
        nlsolve_linesearch,
        nlsolve_ftol,
        nlsolve_iteration_limit,
        verbose,
        cache,
    ) = const_cache

    # - Get fan_params back out using fx(x) in constants - #
    (; fan_params) = fx(differentiated_variables)

    # - Do precomputations - #
    populate_cache!(cache, fan_params)

    # - Initialize Aero - #
    Gamr, gamw, sigr = extract_state_vars(state_variables)

    # - Call nlsolve residual using dummy r - #
    return nls_res_fun!(r, [Gamr; gamw; sigr], cache)
end

"""
This is the residual that gets passed into nlsolve, and is just for converging Gamr and gamw.
"""
function nls_res_fun!(r, state_variables, cache)
    # function nls_res_fun!(r, state_variables, differentiated_variables, cache)

    # r are the residual values the solver is trying to converge by changing the state variables
    # state_variables are the values being changed by the solver to converge the residuals
    # differentiated_variables are any values that may be affected by potential optimization design variables and are NOT changed by the solver.
    # cache are CONSTANT parameters (i.e. unaffected by any potential optimization design variables) and are NOT changed by the solver. though in this case cache is a big container to avoid lots of allocations

    # - Separate out the Circulation and Wake Panel Strengths (USE VIEWS) - #
    # note, cache needs to include the dimension information for number of panels, rotors, etc.
    Gamr, gamw, sigr = extract_state_vars(state_variables, cache)

    # - Zero out State Estimates - #
    cache.Gamr_est .= 0.0
    cache.gamw_est .= 0.0
    cache.sigr_est .= 0.0

    # - Estimate New States - #
    estimate_states!(
        cache.Gamr_est, cache.gamw_est, cache.sigr_est, Gamr, gamw, sigr, cache
    )

    # - Get final Residual Values - #
    @. r = [cache.Gamr_est - Gamr; cache.gamw_est - gamw; cache.sigr_est - sigr]

    return r
end

"""
"""
function estimate_states!(Gamr_est, gamw_est, sigr_est, Gamr, gamw, sigr, cache)

    # TODO: NOTHING IN THE CACHE IS ALLOWED TO BE USED BETWEEN ITERATIONS
    # TODO: ANY UPDATES TO THE CACHE MUST BE BASED ON THINGS DONE IN THIS ITERATION
    # TODO: all the strengths need to be state variables then...

    # - Solve Linear System for [dummy] Body Strengths - #
    # Note: body strengths are a direct function of gamw, and sigr, so they don't need to be state variables
    calculate_body_vortex_strengths!(
        cache.gamb,
        cache.A_bb, # note that the AIC matrices here are constant relative to nlsolve
        cache.b_bf,
        gamw,
        cache.A_bw,
        cache.A_pw,
        sigr,
        cache.A_br,
        cache.A_pr,
        cache.RHS;
    )

    # Update rotor blade element velocities with body influence
    # Note: we only use the values in gamb associated with the node strengths
    # note: could put these velocities inside the cache and do this in place
    # TODO: refactor functions that take in cache (inputs) to more explicitly take in the specific inputs they need
    _, _, _, Wz_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor = calculate_rotor_velocities(
        Gamr, gamw, cache.sigr, cache.gamb[1:(cache.body_vortex_panels.totnode)], cache
    ) # note that the blade elements and such inside cache are constant relative to nlsolve

    # - Calculate Blade Element Values - #
    # note: could put cl and cd inside the cache and do this in place
    cl, cd = calculate_blade_element_coefficients(
        cache.blade_elements, Wz_rotor, Wtheta_rotor, Wmag_rotor, cache.freestream;
    )

    # - Estimate Gamr - #
    # in-place solve for Gamr, updating Gamr_est
    calculate_rotor_circulation_strengths!(Gamr_est, Wmag_rotor, cache.blade_elements, cl)

    # - Calculate Velocities on Wake Panels - #
    # note: could put Wm_wake inside the cache and do this in place
    # TODO: refactor functions that take in cache (inputs) to more explicitly take in the specific inputs they need
    Wm_wake = calculate_wake_velocities(gamw, sigr, cache)

    # - Estimate Wake Panel Strengths - #
    # in-place solve for gamw, update gamw_est
    # TODO: refactor functions that take in cache (inputs) to more explicitly take in the specific inputs they need
    calculate_wake_vortex_strengths!(gamw_est, Gamr, Wm_wake, cache)

    # - Estimate Rotor Panel Strengths - #
    # in-place solver for sigr, update sigr_est
    calculate_rotor_source_strengths!(
        sigr_est, Wmag_rotor, cache.blade_elements, cd, cache.freestream.rhoinf
    )

    return Gamr_est, gamw_est, sigr_est
end

##---------------------------------#
##            RELAXATION           #
##---------------------------------#

#function relax_gamw!(gamw, ; nrf=0.4, btw=0.6, pfw=1.2)

#    # initilize
#    TF = eltype(gamw)
#    omega = Array{TF,0}(undef)
#    omega[] = nrf

#    # use delta_prev as a place holder for relaxation factor criteria
#    # TODO: how to access these?
#    delta_prev .*= delta

#    for ig in eachindex(gamw)

#        # choose relaxation factor based on whether old and new changes are in different or same direction
#        # note that delta_prev at this point = delta_prev .* delta
#        omega[] = sign(delta_prev[ig]) < 0.0 ? btw * nrf : pfw * nrf

#        # save delta_prev for next iteration
#        delta_prev[ig] = delta[ig] * omega[]

#        # update gamw value
#        gamw[ig] += delta_prev[ig]
#    end

#    return gamw
#end

#function relax_Gamr!(Gamr, ; nrf=0.4, bt1=0.2, bt2=0.6, pf1=0.4, pf2=0.5)

#    # initilize
#    TF = eltype(Gamr)
#    bladeomega = nrf .* ones(TF, size(Gamr, 2))
#    omega = nrf .* ones(TF, size(Gamr, 1))
#    deltahat = zeros(TF, size(Gamr, 1))

#    #TODO; how to access B, delta_prev_mat, and delta_mat?
#    #TODO: should they get put inside the constants?

#    for (i, (BG, b, delta_prev, delta)) in
#        enumerate(zip(eachcol(Gamr), B, eachcol(delta_prev_mat), eachcol(delta_mat)))

#        # multiply by B to get BGamr values
#        BG .*= b
#        delta .*= b
#        delta_prev .*= b

#        # - Set the normalization value based on the maximum magnitude value of B*Gamr

#        # find max magnitude
#        maxBGamr[i], mi = findmax(abs.(BG))

#        # maintain sign of original value
#        maxBGamr[i] *= sign(BG[mi])

#        # make sure we don't have any weird jumps
#        meang = sum(BG) / length(BG)
#        if meang > 0.0 # if mean is positive, make sure maxBGamr[i] is at least 0.1
#            maxBGamr[i] = max(maxBGamr[i], 0.1)
#        elseif meang < 0.0 # if mean is negative, make sure maxBGamr[i] is at most -0.1
#            maxBGamr[i] = min(maxBGamr[i], -0.1)
#        else # if the average is zero, then set maxBGamr[i] to zero
#            maxBGamr[i] = 0.0
#        end

#        # note: delta = Gamr_new .- Gamr
#        # note: deltahat here is actually 1/deltahat which is the version needed later
#        for (j, d) in enumerate(eachrow(deltahat))
#            if abs(delta[j]) < eps()
#                d[1] = sign(delta[j]) * sign(maxBGamr[i]) #avoid division by zero
#            else
#                d[1] = maxBGamr[i] ./ delta[j]
#            end
#        end

#        # get initial relaxation factor
#        bladeomega[i], oi = findmin(abs.(deltahat))

#        # scale relaxation factor based on if the change and old values are the same sign (back tracking or pressing forward)
#        if (nrf / deltahat[oi]) < -bt1
#            bladeomega[i] *= bt1
#        elseif (nrf / deltahat[oi]) > pf1
#            bladeomega[i] *= pf1
#        else
#            bladeomega[i] = nrf
#        end

#        # scale blade element relaxation factor
#        for (o, d, dp) in zip(eachrow(omega), delta, eachrow(delta_prev))
#            # if differences changed sign, use backtrack factor, if sign is same, use press forward factor
#            o[1] = bladeomega[i] * (sign(d) != sign(dp[1]) ? bt2 : pf2)

#            # save current delta into old one for next iteration
#            dp[1] = d
#        end

#        # relax Gamr for this blade
#        BG .+= omega .* delta
#        # remove the *b
#        BG ./= b
#        delta ./= b
#        delta_prev ./= b
#    end

#    return Gamr
#end

#"""
#"""
#function relax_state_variables!(state_variables)

#    # - Extract States into Gamr and gamw using views - #

#    # - Relax Gamr - #
#    relax_Gamr!(Gamr)

#    # - Relax gamw - #
#    relax_gamw!(gamw)

#    return state_variables
#end

#"""
#"""
#function converge_state_variables(nls_res_fun!, state_variables)

#    # - Initialize - #
#    TF = eltype(state_variables)
#    r = zeros(TF, 2)
#    conv = [false]
#    iter = 0

#    # initialize differences
#    # TODO: hold these in the constants?
#    # TODO; need these, but how to get them inside the residual and the relaxation functions?
#    deltaG_prev = copy(Gamr)
#    deltaG = similar(deltaG_prev) .= 0.0
#    deltag_prev = copy(gamw)
#    deltag = similar(deltag_prev) .= 0.0

#    # loop until converged or max iterations are reached
#    while !conv[] && iter <= maxiter

#        # update iteration number
#        iter += 1
#        if verbose
#            println("Iteration $(iter):")
#        end

#        # Call Residual
#        nls_res_fun!(
#            r, state_variables, differentiated_variables, constants
#        )

#        # Relax Values
#        # TODO: how to get the required additional variables in here, do they go into constants?
#        relax_state_variables!(state_variables)

#        # Check Convergence
#        check_convergence(r)
#    end
#end

