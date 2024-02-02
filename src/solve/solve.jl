"""
"""
function solve!(
    inputs,
    Gamr,
    Gamr_est,
    sigr,
    sigr_est,
    gamw,
    gamw_est;
    nosource=true,
    maxiter=1e2,
    verbose=false,
    nrf=0.4,
    bt1=0.2,
    bt2=0.6,
    pf1=0.4,
    pf2=0.5,
    btw=0.6,
    pfw=1.2,
    f_circ=1e-3, #DFDC values
    f_dgamw=2e-4, #DFDC Values
    # f_circ=1e-10, # tighter tolerances don't appear to help much in accuracy
    # f_dgamw=2e-12,
    f_sigr=1e-2, # DFDC doesn't converge on sigr
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

    #TODO: check how DFDC does this part.
    # initialize differences
    deltaG_prev = Gamr_est .- Gamr
    deltaG = similar(deltaG_prev) .= 0.0
    deltaS_prev = sigr_est .- sigr
    deltaS = similar(deltaS_prev) .= 0.0
    deltag_prev = gamw_est .- gamw
    deltag = similar(deltag_prev) .= 0.0

    # # zero out sigr before first iteration
    # if nosource
    #     sigr .= 0.0
    # end

    # loop until converged or max iterations are reached
    while !conv[] && iter <= maxiter
        # update iteration number
        iter += 1
        if verbose
            println("Iteration $(iter):")
        end

        ##### ----- Solve linear system if including ----- #####
        if !isnothing(inputs.gamb)
            # in place solve for gamb
            calculate_body_vortex_strengths!(
                inputs.gamb,
                inputs.A_bb,
                inputs.b_bf,
                gamw,
                inputs.A_bw,
                inputs.A_pw,
                sigr,
                inputs.A_br,
                inputs.A_pr,
                inputs.RHS;
                post=false,
            )

            # Update rotor blade element velocities with body influence
            # note: only use the values in gamb associated with the node strengths
            _, _, _, Wz_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor = calculate_rotor_velocities(
                Gamr, gamw, sigr, inputs.gamb[1:(inputs.body_vortex_panels.totnode)], inputs
            )
        else

            # Update rotor blade element velocities without body influence
            _, _, _, Wz_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor = calculate_rotor_velocities(
                Gamr, gamw, sigr, inputs
            )
        end

        ##### ----- Calculate Blade Element Values ----- #####
        # calculate lift and drag coefficients along blades
        cl, cd = calculate_blade_element_coefficients(
            inputs.blade_elements, Wz_rotor, Wtheta_rotor, Wmag_rotor, inputs.freestream;
        )

        ##### ----- Estimate and Relax Gamr ----- #####
        # in-place solve for Gamr, updating Gamr_est
        calculate_rotor_circulation_strengths!(
            Gamr_est, Wmag_rotor, inputs.blade_elements, cl
        )

        # get difference between estimated Gamr and old Gamr
        @. deltaG = Gamr_est - Gamr

        # Keep for now, used in test development
        # println("BGX = ", inputs.blade_elements[1].B*Gamr_est[:,1])
        # println("BGAM = ", inputs.blade_elements[1].B*Gamr[:,1])
        # println("BGMAX = ", maximum(inputs.blade_elements[1].B*Gamr[:,1]))
        # println("NRC = ", length(Gamr))
        # println("DBGOLD =", deltaG_prev)

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

        ##### ----- Estimate and Relax gamw ----- #####

        # Update Vm_avg in wake using new Gamr and sigma
        if !isnothing(inputs.gamb)
            Wm_wake = calculate_wake_velocities(
                gamw, sigr, inputs.gamb[1:(inputs.body_vortex_panels.totnode)], inputs
            )
        else
            Wm_wake = calculate_wake_velocities(gamw, sigr, inputs)
        end

        # in-place solve for gamw, update gamw_est
        calculate_wake_vortex_strengths!(gamw_est, Gamr, Wm_wake, inputs)

        # get difference beetween estimated gamw and old gamw
        @. deltag = gamw_est - gamw

        # Keep for now, used in test development
        # println("GAMTH = ", gamw_est)
        # println("GTH = ", gamw)
        # println("NPTOT = ", length(gamw))
        # println()
        # println("DGOLD =", deltag_prev)

        # relax gamw values
        relax_gamw!(gamw, deltag_prev, deltag, maxdeltagamw; nrf=nrf, btw=btw, pfw=pfw)

        ##### ----- Update sigr ----- #####
        if nosource
            # Update rotor blade element velocities without body influence
            _, _, _, _, _, _, Wmag_rotor = calculate_rotor_velocities(
                Gamr, gamw, sigr, inputs
            )

            # update sigr in place
            calculate_rotor_source_strengths!(
                sigr, Wmag_rotor, inputs.blade_elements, cd, inputs.freestream.rhoinf
            )

        else
            # - Relax Sigr if including in "residual" - #
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

        # converged?
        if nosource
            check_convergence!(
                conv,
                maxBGamr,
                maxdeltaBGamr,
                maxdeltagamw,
                inputs.Vconv[1];
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
                inputs.Vconv[1];
                f_circ=f_circ,
                f_dgamw=f_dgamw,
                f_sigr=f_sigr,
                verbose=verbose,
            )
        end
    end

    inputs.converged[1] = conv[1]
    inputs.iterations[1] = iter

    return Gamr, sigr, gamw
end

"""
# Arguments:
- `Gamr::Array{Float}` : Array of rotor circulations (columns = rotors, rows = blade elements), updated in place
- `delta_prev_mat::Array{Float}` : Array of previous iteration's differences in circulation values, updated in place
- `delta_mat::Array{Float}` : Array of current iteration's differences in circulation values
- `B::Vector{Float}` : number of blades on each rotor
- `nrf::Float=0.4` : nominal relaxation factor
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
    nrf=0.4,
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
        BG .*= b
        delta .*= b
        delta_prev .*= b

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
                d[1] = sign(delta[j]) * sign(maxBGamr[i]) #avoid division by zero
            else
                d[1] = maxBGamr[i] ./ delta[j]
            end
        end

        # get initial relaxation factor
        bladeomega[i], oi = findmin(abs.(deltahat))

        # scale relaxation factor based on if the change and old values are the same sign (back tracking or pressing forward)
        if -0.2 > (nrf / deltahat[oi]) || (nrf / deltahat[oi]) > 0.4
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
        delta ./= b
        delta_prev ./= b
    end

    if test
        return Gamr, bladeomega, omega
    else
        return Gamr
    end
end

"""
# Arguments:
- `sigr::Array{Float}` : Array of rotor circulations (columns = rotors, rows = blade elements), updated in place
- `delta_prev_mat::Array{Float}` : Array of previous iteration's differences in circulation values, updated in place
- `delta_mat::Array{Float}` : Array of current iteration's differences in circulation values
- `B::Vector{Float}` : number of blades on each rotor
- `nrf::Float=0.4` : nominal relaxation factor
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
    nrf=0.4,
    bt1=0.2,
    bt2=0.6,
    pf1=0.4,
    pf2=0.5,
    test=false,
)

    # initilize
    TF = eltype(sigr)
    bladeomega = nrf .* ones(TF, size(sigr, 2))
    omega = nrf .* ones(TF, size(sigr, 1))
    deltahat = zeros(TF, size(sigr, 1))

    for (i, (BS, delta_prev, delta)) in
        enumerate(zip(eachcol(sigr), eachcol(delta_prev_mat), eachcol(delta_mat)))

        # - Set the normalization value based on the maximum magnitude value of B*sigr

        # find max magnitude
        maxsigr[i], mi = findmax(abs.(BS))

        # maintain sign of original value
        maxsigr[i] *= sign(BS[mi])

        # make sure we don't have any weird jumps
        meang = sum(BS) / length(BS)
        if meang > 0.0 # if mean is positive, make sure maxsigr[i] is at least 0.1
            maxsigr[i] = max(maxsigr[i], 0.1)
        elseif meang < 0.0 # if mean is negative, make sure maxsigr[i] is at most -0.1
            maxsigr[i] = min(maxsigr[i], -0.1)
        else # if the average is zero, then set maxsigr[i] to zero
            maxsigr[i] = 0.0
        end

        # note: delta = sigr_new .- sigr
        # note: deltahat here is actually 1/deltahat which is the version needed later
        for (j, d) in enumerate(eachrow(deltahat))
            if abs(delta[j]) < eps()
                d[1] = sign(delta[j]) * sign(maxdBsigr[i]) #avoid division by zero
            else
                d[1] = maxsigr[i] ./ delta[j]
            end
        end

        # get initial relaxation factor
        bladeomega[i], oi = findmin(abs.(deltahat))

        # scale relaxation factor based on if the change and old values are the same sign (back tracking or pressing forward)
        if -0.2 > nrf / (deltahat[oi] * maxsigr[i]) > 0.4
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
        # maxdeltasigr[i] = maximum(delta)
        maxdeltasigr[i], mdi = findmax(abs.(delta))
        # maintain sign of original value
        maxdeltasigr[i] *= sign(delta[mdi])

        # relax sigr for this blade
        BS .+= omega .* delta
    end

    if test
        return sigr, bladeomega, omega
    else
        return sigr
    end
end

"""
# Arguments:
- `gamw::Array{Float}` : Array of rotor circulations (columns = rotors, rows = blade elements), updated in place
- `delta_prev_mat::Array{Float}` : Array of previous iteration's differences in circulation values, updated in place
- `delta_mat::Array{Float}` : Array of current iteration's differences in circulation values
- `B::Vector{Float}` : number of blades on each rotor
- `nrf::Float=0.4` : nominal relaxation factor
- `bt1::Float=0.2` : backtrack factor 1
- `bt2::Float=0.6` : backtrack factor 2
- `pf1::Float=0.4` : press forward factor 1
- `pf2::Float=0.5` : press forward factor 2
"""
function relax_gamw!(
    gamw, delta_prev, delta, maxdeltagamw; nrf=0.4, btw=0.6, pfw=1.2, test=false
)

    # initilize
    TF = eltype(gamw)
    omega = MVector{1,TF}(nrf)

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
    f_sigr=1e-3,
    verbose=false,
)

    # find max ratio among blades and use that for convergence
    _, idG = findmax(maxdeltaBGamr ./ maxBGamr)
    _, idS = findmax(maxdeltasigr ./ maxsigr)

    # set convergence flag
    conv[] =
        abs(maxdeltaBGamr[idG]) < f_circ * abs(maxBGamr[idG]) &&
        abs(maxdeltasigr[idS]) < f_sigr * abs(maxsigr[idS]) &&
        maxdeltagamw[] < f_dgamw * Vref # abs already taken care of
    # abs(maxdeltagamw[]) < f_dgamw * Vref

    # # set convergence flag, note: this is how dfdc does it, without regard to which rotor for the Gamr values
    # conv[] =
    # max(abs.(maxdeltaBGamr)...) < f_circ * max(abs.(maxBGamr)...) &&
    #     maxdeltagamw[] < f_dgamw * Vref

    if verbose
        @printf "\t%-16s      %-16s" "maxdBGamr" "fG*maxBGamr"
        println()
        @printf "\t%1.16f    %1.16f" abs(maxdeltaBGamr[idG]) f_circ * abs(maxBGamr[idG])
        println()
        @printf "\t%-16s      %-16s" "maxdsigr" "fs*maxsigr"
        println()
        @printf "\t%1.16f    %1.16f" abs(maxdeltasigr[idS]) f_sigr * abs(maxsigr[idS])
        println()
        @printf "\t%-16s      %-16s" "maxdgamw" "fg*Vconv"
        println()
        @printf "\t%1.16f    %1.16f" abs(maxdeltagamw[]) f_dgamw * Vref
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

    if verbose
        @printf "\t%-16s      %-16s" "maxdBGamr" "fG*maxBGamr"
        println()
        @printf "\t%1.16f    %1.16f" abs(maxdeltaBGamr[idG]) f_circ * abs(maxBGamr[idG])
        println()
        @printf "\t%-16s      %-16s" "maxdgamw" "fg*Vconv"
        println()
        @printf "\t%1.16f    %1.16f" abs(maxdeltagamw[]) f_dgamw * Vref
        println()
        println()
    end

    return conv
end
