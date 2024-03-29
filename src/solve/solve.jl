"""
"""
function solve!(state_variables, const_cache; verbose=false)

    ### --- Unpack contants and caches --- ###
    (;
        # solve options
        solve_options,
        # Constant Parameters
        airfoils,
        A_bb_LU,
        idmaps,
        state_dims,
        # Cache(s)
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        solve_container_cache,
        solve_container_cache_dims,
    ) = const_cache

    # - Separate out the state variables - #
    Gamr, sigr, gamw = extract_state_variables(solve_options, state_variables, state_dims)

    # - separate out solve parameters here as well - #
    solve_parameter_tuple = withdraw_solve_parameter_cache(
        solve_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )

    # unpack parameters
    (;
        operating_point, # freestream, Omega, etc.
        ivr,             # induced velocities on rotor panels
        ivw,             # induced velocities on wake panels
        linsys,          # includes AIC's for linear system
        blade_elements,  # includes blade element geometry
        wakeK,           # geometric "constants" for wake node strength calculation
    ) = solve_parameter_tuple

    # - Extract and Reset Cache - #
    # get cache vector of correct types
    solve_container_cache_vector = @views PreallocationTools.get_tmp(
        solve_container_cache, eltype(state_variables)(1.0)
    )
    # reset cache
    solve_container_cache_vector .= 0
    solve_containers = withdraw_solve_container_cache(
        solve_options, solve_container_cache_vector, solve_container_cache_dims
    )

    # - initialize convergence criteria - #
    @. solve_containers.deltaG_prev = solve_containers.Gamr_est - Gamr
    @. solve_containers.deltaG = 0.0
    @. solve_containers.deltag_prev = solve_containers.gamw_est - gamw
    @. solve_containers.deltag = 0.0

    TF = eltype(state_variables)

    #Note, using MVectors here is faster than putting these in a cache
    maxBGamr = MVector{1,TF}(0.0)
    maxdeltaBGamr = MVector{1,TF}(0.0)
    maxsigr = MVector{1,TF}(0.0)
    maxdeltasigr = MVector{1,TF}(0.0)
    maxdeltagamw = MVector{1,TF}(0.0)
    conv = solve_options.converged
    # conv = MVector{1,Bool}(false)
    iter = 0

    # loop until converged or max iterations are reached
    while !conv[] && iter <= solve_options.maxiter
        # update iteration number
        iter += 1
        if verbose
            println("Iteration $(iter):")
        end

        # - Solve Linear System - #
        # in place solve for gamb
        calculate_body_vortex_strengths!(
            solve_containers.gamb,
            A_bb_LU,
            linsys.b_bf,
            gamw,
            linsys.A_bw,
            linsys.A_pw,
            sigr,
            linsys.A_br,
            linsys.A_pr,
            linsys.A_bb,
            solve_containers.rhs;
            post=false,
        )

        # - Update rotor blade element velocities with body influence - #
        calculate_induced_velocities_on_rotors!(
            solve_containers.vz_rotor,
            solve_containers.vtheta_rotor,
            Gamr,
            gamw,
            sigr,
            @view(solve_containers.gamb[1:(idmaps.body_totnodes)]),
            @view(ivr.v_rw[:, :, 1]),
            @view(ivr.v_rr[:, :, 1]),
            @view(ivr.v_rb[:, :, 1]),
            blade_elements.B,
            blade_elements.rotor_panel_centers,
        )

        # - Get Absolute Rotor Velocities - #
        reframe_rotor_velocities!(
            solve_containers.Cz_rotor,
            solve_containers.Ctheta_rotor,
            solve_containers.Cmag_rotor,
            solve_containers.vz_rotor,
            solve_containers.vtheta_rotor,
            operating_point.Vinf[],
            operating_point.Omega,
            blade_elements.rotor_panel_centers,
        )

        # - Calculate Blade Element Values - #
        calculate_blade_element_coefficients!(
            solve_containers.cl,
            solve_containers.cd,
            solve_containers.beta1,
            solve_containers.alpha,
            solve_containers.reynolds,
            solve_containers.mach,
            (; blade_elements..., airfoils...),
            solve_containers.Cz_rotor,
            solve_containers.Ctheta_rotor,
            solve_containers.Cmag_rotor,
            operating_point;
            verbose=verbose,
        )

        ##### ----- Estimate and Relax Gamr ----- #####
        # - Calculate Blade Element Circulation - #
        calculate_rotor_circulation_strengths!(
            solve_containers.Gamr_est,
            solve_containers.Cmag_rotor,
            blade_elements.chords,
            solve_containers.cl,
        )

        # - get difference between estimated Gamr and old Gamr - #
        @. solve_containers.deltaG = solve_containers.Gamr_est - Gamr

        # Keep for now, used in test development
        # println("BGX = ", inputs.blade_elements[1].B*Gamr_est[:,1])
        # println("BGAM = ", inputs.blade_elements[1].B*Gamr[:,1])
        # println("BGMAX = ", maximum(inputs.blade_elements[1].B*Gamr[:,1]))
        # println("NRC = ", length(Gamr))
        # println("DBGOLD =", deltaG_prev)

        # - relax Gamr values - #
        relax_Gamr!(
            Gamr,
            solve_containers.deltaG_prev,
            solve_containers.deltaG,
            maxBGamr,
            maxdeltaBGamr,
            blade_elements.B;
            nrf=solve_options.nrf,
            bt1=solve_options.bt1,
            bt2=solve_options.bt2,
            pf1=solve_options.pf1,
            pf2=solve_options.pf2,
        )

        ##### ----- Estimate and Relax gamw ----- #####

        # - Calculate Velocities on Wake Panels - #
        calculate_wake_velocities!(
            solve_containers.Cm_wake,
            solve_containers.vz_wake,
            solve_containers.vr_wake,
            solve_containers.gamw_est,
            sigr,
            @view(solve_containers.gamb[1:(idmaps.body_totnodes)]),
            ivw,
            operating_point.Vinf[],
        )

        # - Get Average Wake Velocities - #
        # currently has 5 allocations
        average_wake_velocities!(
            solve_containers.Cm_avg,
            solve_containers.Cm_wake,
            idmaps.wake_nodemap,
            idmaps.wake_endnodeidxs,
        )

        # - Estimate wake strengths - #
        # in-place solve for gamw, update gamw_est
        calculate_wake_vortex_strengths!(
            solve_containers.gamw_est,
            solve_containers.Gamma_tilde,
            solve_containers.H_tilde,
            solve_containers.deltaGamma2,
            solve_containers.deltaH,
            Gamr,
            solve_containers.Cm_avg,
            blade_elements.B,
            operating_point.Omega,
            wakeK,
            idmaps.wake_node_sheet_be_map,
            idmaps.wake_node_ids_along_casing_wake_interface,
            idmaps.wake_node_ids_along_centerbody_wake_interface;
        )

        # - get difference between estimated gamw and old gamw - #
        @. solve_containers.deltag = solve_containers.gamw_est - gamw

        # Keep for now, used in test development
        # println("GAMTH = ", gamw_est)
        # println("GTH = ", gamw)
        # println("NPTOT = ", length(gamw))
        # println()
        # println("DGOLD =", deltag_prev)

        # relax gamw values
        relax_gamw!(
            gamw,
            solve_containers.deltag_prev,
            solve_containers.deltag,
            maxdeltagamw;
            nrf=solve_options.nrf,
            btw=solve_options.btw,
            pfw=solve_options.pfw,
        )

        ##### ----- Update sigr ----- #####
        # Update rotor blade element velocities without body influence
        # - Update rotor blade element velocities with body influence - #
        calculate_induced_velocities_on_rotors!(
            solve_containers.vz_rotor,
            solve_containers.vtheta_rotor,
            Gamr,
            gamw,
            sigr,
            @view(solve_containers.gamb[1:(idmaps.body_totnodes)]),
            @view(ivr.v_rw[:, :, 1]),
            @view(ivr.v_rr[:, :, 1]),
            @view(ivr.v_rb[:, :, 1]),
            blade_elements.B,
            blade_elements.rotor_panel_centers,
        )

        # - Get Absolute Rotor Velocities - #
        reframe_rotor_velocities!(
            solve_containers.Cz_rotor,
            solve_containers.Ctheta_rotor,
            solve_containers.Cmag_rotor,
            solve_containers.vz_rotor,
            solve_containers.vtheta_rotor,
            operating_point.Vinf[],
            operating_point.Omega,
            blade_elements.rotor_panel_centers,
        )

        # - Calculate Rotor Panel Strengths - #
        calculate_rotor_source_strengths!(
            sigr,
            solve_containers.Cmag_rotor,
            blade_elements.chords,
            blade_elements.B,
            solve_containers.cd,
        )

        # converged?
        check_convergence!(
            conv,
            maxBGamr,
            maxdeltaBGamr,
            maxdeltagamw;
            Vconv=solve_options.Vconv[1],
            use_abstol=solve_options.use_abstol,
            f_circ=solve_options.f_circ,
            f_dgamw=solve_options.f_dgamw,
            verbose=verbose,
        )
    end

    return state_variables
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
        # if -bt1 > (nrf / deltahat[oi]) || (nrf / deltahat[oi]) > pf1
        #     bladeomega[i] *= sign(deltahat[oi]) < 0.0 ? bt1 : pf1
        # else
        #     bladeomega[i] = nrf
        # end

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
    # omega = TF[nrf]

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

function check_convergence!(
    conv,
    maxBGamr,
    maxdeltaBGamr,
    maxdeltagamw;
    use_abstol=true,
    Vconv=1.0,
    f_circ=1e-3,
    f_dgamw=2e-4,
    verbose=false,
)

    # find max ratio among blades and use that for convergence
    _, idG = findmax(maxdeltaBGamr ./ maxBGamr)

    # set convergence flag
    if use_abstol
        conv[] = abs(maxdeltaBGamr[idG]) < f_circ && abs(maxdeltagamw[]) < f_dgamw
    else
        conv[] =
            abs(maxdeltaBGamr[idG]) < f_circ * abs(maxBGamr[idG]) &&
            abs(maxdeltagamw[]) < f_dgamw * Vconv
    end

    # # set convergence flag, note: this is how dfdc does it, without regard to which rotor for the Gamr values
    # conv[] =
    # max(abs.(maxdeltaBGamr)...) < f_circ * max(abs.(maxBGamr)...) &&
    #     maxdeltagamw[] < f_dgamw * Vconv

    if verbose
        if use_abstol
            @printf "\t%-16s      %-16s" "maxdBGamr" "fG"
            println()
            @printf "\t%1.16f    %1.16f" abs(maxdeltaBGamr[idG]) f_circ
            println()
            @printf "\t%-16s      %-16s" "maxdgamw" "fg"
            println()
            @printf "\t%1.16f    %1.16f" abs(maxdeltagamw[]) f_dgamw
            println()
            println()
        else
            @printf "\t%-16s      %-16s" "maxdBGamr" "fG*maxBGamr"
            println()
            @printf "\t%1.16f    %1.16f" abs(maxdeltaBGamr[idG]) f_circ * abs(maxBGamr[idG])
            println()
            @printf "\t%-16s      %-16s" "maxdgamw" "fg*Vconv"
            println()
            @printf "\t%1.16f    %1.16f" abs(maxdeltagamw[]) f_dgamw * Vconv
            println()
            println()
        end
    end

    return conv
end
