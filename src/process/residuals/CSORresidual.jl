"""
"""
function CSOR_residual!(resid, state_variables, sensitivity_parameters, constants)

    # - Extract constants - #
    (;
        # solve options for dispatch
        solver_options,
        # comp,
        airfoils,                   # airfoils
        A_bb_LU,                    # linear system left hand side LU decomposition
        idmaps,                     # book keeping items
        solve_parameter_cache_dims, # dimensions for shaping the view of the parameter cache
        solve_container_cache,      # cache for solve_containers used in solve
        solve_container_cache_dims, # dimensions for shaping the view of the solve cache
    ) = constants

    # - Separate out the state variables - #
    Gamr, sigr, gamw = extract_state_variables(
        solver_options, state_variables, solve_parameter_cache_dims.state_dims
    )

    # separate out sensitivity_parameters here as well
    solve_parameter_tuple = withdraw_solve_parameter_cache(
        solver_options, sensitivity_parameters, solve_parameter_cache_dims
    )

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
        solve_container_cache,
        promote_type(eltype(state_variables), eltype(sensitivity_parameters))(1.0),
    )
    # zero out contents of solve_containers to avoid any potential contamination issues
    solve_container_cache_vector .= 0
    solve_containers = withdraw_solve_container_cache(
        solver_options, solve_container_cache_vector, solve_container_cache_dims
    )

    # - Estimate New States - #
    compute_CSOR_residual!(
        resid,
        solver_options,
        solve_containers,
        Gamr,
        sigr,
        gamw,
        operating_point,
        ivr,
        ivw,
        (; linsys..., A_bb_LU),
        (; blade_elements..., airfoils...),
        wakeK,
        idmaps;
        verbose=solver_options.verbose,
    )

    return state_variables
end

"""
"""
function compute_CSOR_residual!(
    resid,
    solver_options,
    solve_containers,
    Gamr,
    sigr,
    gamw,
    operating_point,
    ivr,
    ivw,
    linsys,
    blade_elements,
    wakeK,
    idmaps;
    verbose=false,
)

    # - Adjust Relaxation Factors (decrease relaxation as solution converges) - #
    nrf, bt1, bt2, pf1, pf2, nrfw, btw, pfw = apply_relaxation_schedule(
        resid, solver_options
    )

    # - Rename for Convenience - #
    (; maxBGamr, maxdeltaBGamr, maxdeltagamw) = solve_containers

    # - Solve Linear System - #
    # in place solve for gamb
    calculate_body_vortex_strengths!(
        solve_containers.gamb,
        linsys.A_bb_LU,
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
        blade_elements,
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
        nrf=nrf,
        bt1=bt1,
        bt2=bt2,
        pf1=pf1,
        pf2=pf2,
    )

    ##### ----- Estimate and Relax gamw ----- #####

    # - Calculate Velocities on Wake Panels - #
    calculate_wake_velocities!(
        solve_containers.Cm_wake,
        solve_containers.vz_wake,
        solve_containers.vr_wake,
        gamw,
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
        nrf=nrfw,
        btw=btw,
        pfw=pfw,
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
    # TODO: Cmag_rotor doesn't quite match up here with previous version, why?
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

    update_CSOR_residual_values!(
        solver_options.convergence_type,
        resid,
        maxBGamr,
        maxdeltaBGamr,
        maxdeltagamw,
        solver_options.Vconv,
    )

    return nothing
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

"""
"""
function apply_relaxation_schedule(
    resid::AbstractArray, solver_options::TS
) where {TS<:SolverOptionsType}
    # Apply relaxation schedule to Circulation relaxation factors
    nrf = apply_relaxation_schedule(
        resid[1], solver_options.nrf, solver_options.relaxation_schedule
    )

    bt1 = solver_options.bt1 * nrf / solver_options.nrf
    bt2 = solver_options.bt2 * nrf / solver_options.nrf
    pf1 = solver_options.pf1 * nrf / solver_options.nrf
    pf2 = solver_options.pf2 * nrf / solver_options.nrf

    # Apply relaxation schedule to wake strength relaxation factors
    nrfw = apply_relaxation_schedule(
        resid[2], solver_options.nrf, solver_options.relaxation_schedule
    )

    btw = solver_options.btw * nrf / solver_options.nrf
    pfw = solver_options.pfw * nrf / solver_options.nrf

    return nrf, bt1, bt2, pf1, pf2, nrfw, btw, pfw
end

"""
"""
function apply_relaxation_schedule(resid, nominal, schedule)
    rf = linear_transform(
        (0, 1), (nominal, 1), FLOWMath.linear(schedule[1], schedule[2], resid)
    )

    return rf
end

"""
"""
function update_CSOR_residual_values!(
    convergence_type::Relative, resid, maxBGamr, maxdeltaBGamr, maxdeltagamw, Vconv
)
    resid[1] = maximum(abs, maxdeltaBGamr ./ maxBGamr)
    resid[2] = abs.(maxdeltagamw[] / Vconv[])

    return resid
end

"""
"""
function update_CSOR_residual_values!(
    convergence_type::Absolute, resid, maxBGamr, maxdeltaBGamr, maxdeltagamw, Vconv
)
    resid[1] = maximum(abs, maxdeltaBGamr)
    resid[2] = abs.(maxdeltagamw[])

    return resid
end

function check_CSOR_convergence!(
    conv, resid; f_circ=1e-3, f_dgamw=2e-4, convergence_type=Relative(), verbose=false
)

    # set convergence flag
    conv[] = resid[1] < f_circ && resid[2] < f_dgamw

    if verbose
        if typeof(convergence_type) <: Relative
            @printf "\t%-16s      %-16s" "max(ΔBGamr/BGamr)" "fG"
        else
            @printf "\t%-16s      %-16s" "max(ΔBGamr)" "fG"
        end
        println()
        @printf "\t%1.16f    %1.16f" maximum(resid[1:(end - 1)]) f_circ
        println()
        if convergence_type <: Relative
            @printf "\t%-16s      %-16s" "max(Δgamw/Vconv)" "fg"
        else
            @printf "\t%-16s      %-16s" "max(Δgamw)" "fg"
        end
        println()
        @printf "\t%1.16f    %1.16f" resid[end] f_dgamw
        println()
        println()
    end

    return conv
end

