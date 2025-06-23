"""
    CSOR_residual!(
        resid,
        state_variables,
        sensitivity_parameters,
        constants
    )

Compute the in-place residual vector for the CSOR (Coupled Surface and Rotor) solver method.

This function extracts the current state variables (circulations, source strengths, wake strengths),
updates aerodynamic states and induced velocities, and calculates the nonlinear residual used
for iterative solution.

# Arguments
- `resid::Vector{Float}` : Residual vector to be updated in-place.
- `state_variables::Vector{Float}` : Current solution state variables (e.g., circulation, source strengths, wake strengths).
- `sensitivity_parameters::Vector{Float}` : Parameters with respect to which sensitivities or derivatives are computed.
- `constants::NamedTuple` : Contains solver options, linear system data, geometry, and cache structures needed for evaluation.

# Returns
- `resid::Vector{Float}` : The updated residual vector (modified in place).

# Notes
- Handles type dispatch between forward-mode autodiff variables and standard Float64 state variables.
- Uses cached memory containers to avoid repeated allocations.
- Calls `compute_CSOR_residual!` internally to perform the full nonlinear residual evaluation.
- Supports multipoint or multi-parameter solver setups via `multipoint_index`.
"""
function CSOR_residual!(resid, state_variables, sensitivity_parameters, constants)

    # - Extract constants - #
    (;
        # solve options for dispatch
        solver_options,
        # comp,
        A_bb_LU,                    # linear system left hand side LU decomposition
        idmaps,                     # book keeping items
        solve_parameter_cache_dims, # dimensions for shaping the view of the parameter cache
        solve_container_cache,      # cache for solve_containers used in solve
        solve_container_cache_dims, # dimensions for shaping the view of the solve cache
        multipoint_index,
    ) = constants

    # - Separate out the state variables - #
    if eltype(sensitivity_parameters) != eltype(state_variables) &&
        eltype(state_variables) <: Float64
        for k in propertynames(solve_parameter_cache_dims)
            println(k, ": ", getfield(solve_parameter_cache_dims, k))
        end
        println("state var type: ", eltype(state_variables))
        println("param type: ", eltype(sensitivity_parameters))

        svar = extract_initial_guess(
            solver_options, sensitivity_parameters, solve_parameter_cache_dims.state_dims
        )
        csvar = copy(svar)
        Gamr, sigr, gamw = extract_state_variables(
            solver_options, svar, solve_parameter_cache_dims.state_dims
        )
        Gamrval, sigrval, gamwval = extract_state_variables(
            solver_options, state_variables, solve_parameter_cache_dims.state_dims
        )

        Gamr .= Gamrval
        sigr .= sigrval
        gamw .= gamwval

    else
        Gamr, sigr, gamw = extract_state_variables(
            solver_options, state_variables, solve_parameter_cache_dims.state_dims
        )
    end

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
        blade_elements,
        wakeK,
        idmaps,
        multipoint_index;
        verbose=solver_options.verbose,
    )

    if eltype(sensitivity_parameters) != eltype(state_variables) &&
        eltype(state_variables) <: Float64
        state_variables .= ForwardDiff.value.(svar)
        svar .= csvar
    end

    return resid
end

"""
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
        linsys,
        blade_elements,
        wakeK,
        idmaps,
        multipoint_index;
        verbose=false,
    )

Compute and update the residual vector for the coupled surface and rotor (CSOR) panel method solver.

This function performs a full update cycle including:
- Solving for body vortex strengths.
- Calculating induced velocities on rotors and wakes.
- Computing blade element aerodynamic coefficients.
- Estimating and relaxing rotor circulation (`Gamr`) and wake strengths (`gamw`).
- Updating rotor source strengths (`sigr`).
- Updating the convergence residual vector.

# Arguments
- `resid::Vector{Float}` : Residual vector to be updated, used in solver convergence checks.
- `solver_options::SolverOptionsType` : Contains solver settings and relaxation parameters.
- `solve_containers::NamedTuple` : Cache for intermediate solution variables and temporary arrays.
- `Gamr` : Matrix of blade element circulation strengths.
- `sigr` : Matrix of rotor source panel strengths.
- `gamw` : Vector or matrix of wake vortex strengths.
- `operating_point::NamedTuple` : Operating conditions (e.g., inflow velocity, rotation rate).
- `ivr::NamedTuple` : Unit induced velocity influence coefficients on rotors.
- `ivw::NamedTuple` : Unit induced velocity influence coefficients on wakes.
- `linsys::NamedTuple` : Linear system matrices and vectors for panel method solution.
- `blade_elements::NamedTuple` : Geometry and aerodynamic data of blade elements.
- `wakeK::Vector{Float}` : Geometric constants related to wake panels.
- `idmaps::NamedTuple` : Index maps for body, rotor, and wake panels.
- `multipoint_index` : Indexing for multipoint or multiphase solver setups.

# Keyword Arguments
- `verbose::Bool=false` : Enable verbose debug printing.

# Returns
- `nothing` (updates residual and related containers in place)

# Notes
- The function orchestrates a nonlinear iteration step, including relaxation of circulation and wake strengths.
- Residual updates use the specified convergence criteria from `solver_options`.
- Intended for use inside a nonlinear solver loop for iterative solution of the coupled rotor-wake-body system.
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
    idmaps,
    multipoint_index;
    verbose=false,
)

    # - Adjust Relaxation Factors (decrease relaxation as solution converges) - #
    # nrf, bt1, bt2, pf1, pf2, nrfw, btw, pfw = apply_relaxation_schedule(
    #     resid, solver_options
    # )

    (;nrf, bt1, bt2, pf1, pf2, btw, pfw) = solver_options
    nrfw = nrf

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
        blade_elements
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
        idmaps.wake_node_ids_along_center_body_wake_interface;
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
        solver_options.Vconv[multipoint_index[]],
    )

    return nothing
end

"""
    relax_Gamr!(
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
    ) -> Matrix{Float}

Apply adaptive relaxation to rotor circulation (`Gamr`) using previous and current update directions. The method supports per-rotor scaling and convergence monitoring via `B * Gamr`.

# Arguments
- `Gamr::Matrix{Float}` : Rotor circulation matrix, updated in-place. Dimensions: (blade elements × rotors).
- `delta_prev_mat::Matrix{Float}` : Previous update differences (`ΔΓ`), updated in-place for use in the next iteration.
- `delta_mat::Matrix{Float}` : Current iteration's update to circulation values.
- `maxBGamr::Vector{Float}` : Output vector storing the maximum value of `B * Gamr` for each rotor (used for normalization).
- `maxdeltaBGamr::Vector{Float}` : Output vector storing the maximum change in `B * Gamr` for convergence monitoring.
- `B::Vector{Float}` : Number of blades per rotor.

# Keyword Arguments
- `nrf::Float=0.4` : Nominal relaxation factor (base scaling).
- `bt1::Float=0.2` : Backtrack factor 1 (used to reduce step size when update direction is unstable).
- `bt2::Float=0.6` : Backtrack factor 2 (used to reduce elementwise update if signs differ).
- `pf1::Float=0.4` : Press-forward factor 1 (used when update direction is stable).
- `pf2::Float=0.5` : Press-forward factor 2 (elementwise enhancement if direction remains consistent).
- `test::Bool=false` : If true, returns additional diagnostic values (`bladeomega`, `omega`).

# Returns
- If `test=false` (default): returns updated `Gamr`.
- If `test=true`: returns a tuple `(Gamr, bladeomega, omega)`, where:
  - `bladeomega::Vector{Float}` : Relaxation factor used per rotor.
  - `omega::Matrix{Float}` : Relaxation factor used per blade element.
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
    minBGamr = zeros(TF, size(Gamr, 2))

    for (i, (BG, b, delta_prev, delta)) in
        enumerate(zip(eachcol(Gamr), B, eachcol(delta_prev_mat), eachcol(delta_mat)))

        # multiply by B to get BGamr values
        BG .*= b
        delta .*= b
        delta_prev .*= b

        # - Set the normalization value based on the maximum magnitude value of B*Gamr

        # find max magnitude
        maxBGamr[i], maxi = findmax(BG)
        minBGamr[i], mini = findmin(BG)

        # # maintain sign of original value
        # maxBGamr[i] *= sign(BG[maxi])

        # make sure we don't have any weird jumps
        meang = sum(BG) / length(BG)

        if meang > 0.0 # if mean is positive, make sure maxBGamr[i] is at least 0.1
            maxBGamr[i] = max(maxBGamr[i], 0.1)
        elseif meang < 0.0 # if mean is negative, make sure maxBGamr[i] is at most -0.1
            maxBGamr[i] = min(minBGamr[i], -0.1)
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
    relax_gamw!(
        gamw,
        delta_prev,
        delta,
        maxdeltagamw;
        nrf=0.4,
        btw=0.6,
        pfw=1.2,
        test=false
    ) -> Array{Float64}

Apply a relaxed update to the wake circulation strengths (`gamw`) using an adaptive relaxation strategy.

This function modifies `gamw` in place based on the current and previous update directions to improve convergence. It optionally returns the current relaxation factor for testing or diagnostics.

# Arguments
- `gamw::AbstractArray{<:Real}` : Current wake circulation values, updated in place.
- `delta_prev::AbstractArray{<:Real}` : Storage of previous update steps. Will be overwritten with the new scaled updates.
- `delta::AbstractArray{<:Real}` : Current change in `gamw` from solver or residual.
- `maxdeltagamw::Ref{<:Real}` : A scalar reference holding the maximum change in `gamw`, updated in-place for convergence monitoring.

# Keyword Arguments
- `nrf::Float64=0.4` : Nominal relaxation factor.
- `btw::Float64=0.6` : Backtracking multiplier (applied when update direction changes).
- `pfw::Float64=1.2` : Press-forward multiplier (applied when update direction remains consistent).
- `test::Bool=false` : If true, also returns the final relaxation factor used (`omega`).

# Returns
- `gamw::Array` : Updated `gamw` array (also modified in place).
- If `test=true`, also returns:
  - `omega::SVector{1, Float64}` : Last relaxation factor used (mainly for diagnostics).
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

# """
#     apply_relaxation_schedule(
#         resid::AbstractVector, solver_options::TS
#     ) where {TS<:SolverOptionsType}
#
# Apply custom relaxation schedule to all relaxation factor inputs based on residual values.
#
# # Arguments
# - `resid::AbstractVector{Float}` : current residual values
# - `solver_options::SolverOptionsType` : SolverOptions containing relaxation schedule
#
# # Returns
# - `nrf::Float` : nominal relaxation factor
# - `bt1::Float` : backtrack factor 1
# - `bt2::Float` : backtrack factor 2
# - `pf1::Float` : press forward factor 1
# - `pf2::Float` : press forward factor 2
# """
# function apply_relaxation_schedule(
#     resid::AbstractArray, solver_options::TS
# ) where {TS<:SolverOptionsType}
#     # Apply relaxation schedule to Circulation relaxation factors
#     nrf = apply_relaxation_schedule(
#         resid[1], solver_options.nrf, solver_options.relaxation_schedule
#     )
#
#     bt1 = solver_options.bt1 * nrf / solver_options.nrf
#     bt2 = solver_options.bt2 * nrf / solver_options.nrf
#     pf1 = solver_options.pf1 * nrf / solver_options.nrf
#     pf2 = solver_options.pf2 * nrf / solver_options.nrf
#
#     # Apply relaxation schedule to wake strength relaxation factors
#     nrfw = apply_relaxation_schedule(
#         resid[2], solver_options.nrf, solver_options.relaxation_schedule
#     )
#
#     btw = solver_options.btw * nrf / solver_options.nrf
#     pfw = solver_options.pfw * nrf / solver_options.nrf
#
#     return nrf, bt1, bt2, pf1, pf2, nrfw, btw, pfw
# end
#
# """
#     apply_relaxation_schedule(resid, nominal, schedule)
#
# Apply custom relaxation schedule to a single relaxation factor input.
#
# # Arguments
# - `resid::Float` : residual value
# - `nominal::Float` : nominal relaxation value
# - `schedule::AbstractVector{AbstractVector{Float}}` : values between which to interpolate to scale the nominal relaxation value.
#
# # Returns
# - `rf::Float` : the updated relaxation factor
# """
# function apply_relaxation_schedule(resid, nominal, schedule)
#     rf = linear_transform(
#         (0, 1), (nominal, 1), FLOWMath.linear(schedule[1], schedule[2], resid)
#     )
#
#     return rf
# end

"""
    update_CSOR_residual_values!(
        convergence_type::ConvergenceType, resid, maxBGamr, maxdeltaBGamr, maxdeltagamw, Vconv
    )

Update CSOR residual values in place.

# Arguments
- `convergence_type::ConvergenceType` : used for dispatch of relative or absolute residual values.
- `resid::Vector{Float}` : residual values modified in place
- `maxBGamr::Float` : Maximum value of B*Gamr among all blade elements
- `maxdeltaBGamr::Float` : Maximum change in B*Gamr between iterations among all blade elements
- `maxdeltagamw::Vector{Float}` : Maximum change in gamw among all wake nodes (one element)
- `Vconv::Float` : Reference velocity upon which the relative convergence criteria is based (one element)
"""
function update_CSOR_residual_values!(
    convergence_type::Relative, resid, maxBGamr, maxdeltaBGamr, maxdeltagamw, Vconv
)
    resid[1] = maximum(abs, maxdeltaBGamr ./ maxBGamr)
    resid[2] = abs.(maxdeltagamw[] / Vconv[])

    return resid
end

"""
    update_CSOR_residual_values!(convergence_type::Absolute, resid, maxBGamr, maxdeltaBGamr, maxdeltagamw, Vconv)

In place version of update_CSOR_residual_values.
"""
function update_CSOR_residual_values!(
    convergence_type::Absolute, resid, maxBGamr, maxdeltaBGamr, maxdeltagamw, Vconv
)
    resid[1] = maximum(abs, maxdeltaBGamr)
    resid[2] = abs.(maxdeltagamw[])

    return resid
end

"""
    check_CSOR_convergence!(
        conv, resid; f_circ=1e-3, f_dgamw=2e-4, convergence_type=Relative(), verbose=false
    )

Description

# Arguments
- `conv::Vector{Float}` : container holding convergence flag
- `resid::Vector{Float}` : residual vector

# Keyword Arguments
- `f_circ::Float=1e-3` : convergence criteria for circulation residual
- `f_dgamw::Float=2e-4` : convergence criteria for wake strength residual
- `convergence_type::ConvergenceType=Relative()` : convergence type (absolute or relative) for print statements
- `verbose::Bool=false` : flag for verbose print statements
"""
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