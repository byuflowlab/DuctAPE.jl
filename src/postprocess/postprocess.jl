"""
"""
function run_residual!(
    solver_options::TS,
    state_variables,
    state_dims,
    solve_container_cache,
    solve_container_cache_dims,
    operating_point,
    ivr,
    ivw,
    linsys,
    blade_elements,
    wakeK,
    idmaps,
) where {TS<:ExternalSolverOptions}

    #=
      NOTE: we want to get all the intermediate values available to user if desired.
      The solve_containers cache will contain all the intermediate values after running the estimate states function.
    =#
    # - Separate out the state variables - #
    vz_rotor, vtheta_rotor, Cm_wake = extract_state_variables(
        solver_options, state_variables, state_dims
    )

    # - Extract and Reset Cache - #
    # get cache vector of correct types
    solve_container_cache_vec = @views PreallocationTools.get_tmp(
        solve_container_cache, state_variables
    )
    solve_containers = withdraw_solve_container_cache(
        solver_options, solve_container_cache_vec, solve_container_cache_dims
    )
    reset_containers!(solve_containers) #note: also zeros out state estimates

    # - Estimate New States - #
    # currently has 280 allocations
    estimate_states!(
        solve_containers,
        vz_rotor,
        vtheta_rotor,
        Cm_wake,
        operating_point,
        ivr,
        ivw,
        linsys,
        blade_elements,
        wakeK,
        idmaps,
    )

    return (; vz_rotor, vtheta_rotor, Cm_wake, solve_containers...)
end

"""
"""
function run_residual!(
    solver_options::CSORSolverOptions,
    state_variables,
    state_dims,
    solve_container_cache,
    solve_container_cache_dims,
    operating_point,
    ivr,
    ivw,
    linsys,
    blade_elements,
    wakeK,
    idmaps,
)

    #=
      NOTE: we want to get all the intermediate values available to user if desired.
      The solve_containers cache will contain all the intermediate values after running the insides of the residual function
    =#
    # - Separate out the state variables - #
    Gamr, sigr, gamw = extract_state_variables(solver_options, state_variables, state_dims)

    # - Extract and Reset Cache - #
    # get cache vector of correct types
    solve_container_cache_vec = @views PreallocationTools.get_tmp(
        solve_container_cache, state_variables
    )
    solve_container_cache_vec .= 0
    solve_containers = withdraw_solve_container_cache(
        solver_options, solve_container_cache_vec, solve_container_cache_dims
    )

    # - Run Residual - #
    compute_CSOR_residual!(
        zeros(2),
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

    return (; Gamr, sigr, gamw, solve_containers...)
end

"""
"""
function post_process(
    solver_options,
    solve_container_caching,
    state_variables,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
    operating_point,
    reference_parameters,
    ivb,
    A_bb_LU,
    airfoils,
    panels,
    idmaps,
    problem_dimensions;
    write_outputs=false,
    outfile="outputs.jl",
    checkoutfileexists=false,
    output_tuple_name="outs",
    verbose=false,
)

    ### --- SETUP --- ###

    # - Extract Operating Point - #
    (; Vinf, asound, muinf, rhoinf, Omega) = operating_point

    # - Extract Reference Parameters - #
    (; Vref, Rref) = reference_parameters

    # - Extract Panels - #
    (; body_vortex_panels, rotor_source_panels, wake_vortex_panels) = panels

    # rename rotor panel lengths
    rotor_panel_lengths = rotor_source_panels.influence_length

    # - Extract Solve Parameter Cache - #
    solve_parameter_tuple = withdraw_solve_parameter_cache(
        solver_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )
    (; ivr, ivw, linsys, blade_elements, wakeK) = solve_parameter_tuple

    # put airfoils in blade elements and LU decomp into linsys
    blade_elements = (; blade_elements..., airfoils...)
    linsys = (; linsys..., A_bb_LU)

    # - Extract Solve Container Cache - #
    (; solve_container_cache, solve_container_cache_dims) = solve_container_caching

    # - Run Residual to get intermediate values - #
    res_vals = run_residual!(
        solver_options,
        state_variables,
        solve_parameter_cache_dims.state_dims,
        solve_container_cache,
        solve_container_cache_dims,
        operating_point,
        ivr,
        ivw,
        linsys,
        blade_elements,
        wakeK,
        idmaps,
    )

    (;
        Gamr,
        sigr,
        gamw,
        gamb,
        alpha,
        beta1,
        vz_rotor,
        Cz_rotor,
        vtheta_rotor,
        Ctheta_rotor,
        Cmag_rotor,
        cl,
        cd,
        Gamma_tilde,
        Cm_wake,
    ) = res_vals

    ### ----- ROTOR OUTPUTS ----- ###
    # rename for convenience
    (; nbe, nrotor) = problem_dimensions
    rotor_panel_lengths = reshape(rotor_panel_lengths, (nbe, nrotor))
    rotor_panel_centers = blade_elements.rotor_panel_centers
    B = blade_elements.B
    chords = blade_elements.chords

    # - Rotor Thrust - #
    # inviscid thrust
    rotor_inviscid_thrust, rotor_inviscid_thrust_dist = inviscid_rotor_thrust(
        Ctheta_rotor, Gamma_tilde, rotor_panel_lengths, rhoinf[1]
    )

    # viscous thrust
    rotor_viscous_thrust, rotor_viscous_thrust_dist = viscous_rotor_thrust(
        Cz_rotor, Cmag_rotor, B, chords, rotor_panel_lengths, cd, rhoinf[1]
    )

    # total thrust
    rotor_thrust = rotor_inviscid_thrust .+ rotor_viscous_thrust

    # - Rotor Torque - #
    # inviscid torque
    rotor_inviscid_torque, rotor_inviscid_torque_dist = inviscid_rotor_torque(
        Cz_rotor, rotor_panel_centers, rotor_panel_lengths, Gamma_tilde, rhoinf[1]
    )

    # viscous torque
    rotor_viscous_torque, rotor_viscous_torque_dist = viscous_rotor_torque(
        Ctheta_rotor,
        Cmag_rotor,
        B,
        chords,
        rotor_panel_centers,
        rotor_panel_lengths,
        cd,
        rhoinf[1],
    )

    # total torque
    rotor_torque = rotor_inviscid_torque .+ rotor_viscous_torque

    # - Rotor Power - #
    # inviscid power
    rotor_inviscid_power = inviscid_rotor_power(rotor_inviscid_torque', Omega)

    rotor_inviscid_power_dist = similar(rotor_inviscid_torque_dist) .= 0.0
    for ir in 1:nrotor
        rotor_inviscid_power_dist[:, ir] = inviscid_rotor_power(
            @view(rotor_inviscid_torque_dist[:, ir]), Omega[ir]
        )
    end

    # viscous power
    rotor_viscous_power = viscous_rotor_power(rotor_viscous_torque', Omega)
    rotor_viscous_power_dist = similar(rotor_viscous_torque_dist) .= 0.0
    for ir in 1:nrotor
        rotor_inviscid_power_dist[:, ir] = viscous_rotor_power(
            @view(rotor_inviscid_torque_dist[:, ir]), Omega[ir]
        )
    end

    # total power
    rotor_power = rotor_inviscid_power .+ rotor_viscous_power

    # - Rotor Performance Coefficients - #
    rotor_CT, rotor_CQ, rotor_CP = tqpcoeff(
        rotor_thrust, rotor_torque, rotor_power, rhoinf[1], Omega, Rref[1]
    )

    # - Rotor Efficiency - #
    rotor_efficiency = get_total_efficiency(rotor_thrust, rotor_power, Vinf[1])

    # - Blade Loading - #
    blade_normal_force_per_unit_span, blade_tangential_force_per_unit_span = get_blade_loads(
        Cmag_rotor, beta1, cl, cd, chords, rhoinf[1]
    )

    ### --- BODY OUTPUTS --- ###
    # - Surface Velocity on Bodies - #
    # TODO: update any tests for this function
    vtan_tuple = get_body_tangential_velocities(
        gamb,
        gamw,
        sigr,
        ivb,
        Vinf[1],
        body_vortex_panels.totnode,
        body_vortex_panels.totpanel,
        body_vortex_panels.nnode,
        body_vortex_panels.npanel,
        body_vortex_panels.tangent,
        body_vortex_panels.controlpoint,
        body_vortex_panels.endpanelidxs,
        idmaps.wake_node_ids_along_centerbody_wake_interface,
        idmaps.wake_node_ids_along_casing_wake_interface,
        idmaps.centerbody_panel_ids_along_centerbody_wake_interface,
        idmaps.duct_panel_ids_along_centerbody_wake_interface,
    )

    # extract surface velocity outputs
    (;
        # Totals and Components:
        Vtot_in,             # total velocity inside bodies
        Vtot_out,            # total velocity outside bodies
        Vtan_in,             # tangent velocity inside bodies
        Vtan_out,            # tangent velocity outside bodies
        Vtot_prejump,        # total velocity before boundary jumps are added
        vtot_body,           # total velocity induced by bodies
        vtot_jump,           # total velocity due to jump across boundary
        vtot_wake,           # total velocity induced by wake
        vtot_rotors,         # total velocity induced by rotors
        # Splits:
        vtan_casing_in,      # tangent velocity along casing inside of duct body
        vtan_casing_out,     # tangent velocity along casing outside of duct body
        casing_zpts,         # axial locations of control points along casing
        vtan_nacelle_in,     # tangent velocity along nacelle inside of duct body
        vtan_nacelle_out,    # tangent velocity along nacelle outside of duct body
        nacelle_zpts,        # axial locations of control points along nacelle
        vtan_centerbody_in,  # tangent velocity inside of centerbody
        vtan_centerbody_out, # tangent velocity outside of centerbody
        centerbody_zpts,     # axial locations of control points along centerbody
    ) = vtan_tuple

    # - Surface Pressure on Bodies - #
    cp_tuple = get_body_cps(
        Vtan_in,
        Vtan_out,
        Gamr,
        sigr,
        Cz_rotor,
        Vinf[1],
        Vref[1],
        B,
        Omega,
        idmaps.id_of_first_casing_panel_aft_of_each_rotor,
        idmaps.id_of_first_centerbody_panel_aft_of_each_rotor,
        body_vortex_panels.controlpoint,
        body_vortex_panels.endpanelidxs,
    )

    # extract surface pressure outputs
    (;
        cp_in,             # surface pressure along inside of bodies
        cp_out,            # surface pressure along outside of bodies
        cp_casing_in,      # surface pressure along inside of casing
        cp_nacelle_in,     # surface pressure along inside of nacell
        cp_centerbody_in,  # surface pressure along inside of centerbody
        cp_casing_out,     # surface pressure along outside of casing
        cp_nacelle_out,    # surface pressure along outside of nacelle
        cp_centerbody_out, # surface pressure along outside of centerbody
    ) = cp_tuple

    # - Calculate Thrust from Bodies - #
    body_thrust, body_force_coeff = forces_from_pressure(
        cp_in, cp_out, body_vortex_panels; rhoinf=rhoinf[1], Vref=Vref[1]
    )

    # add thrust from trailing edge panels on bodies
    forces_from_TEpanels!(
        body_thrust,
        body_force_coeff,
        cp_in,
        cp_out,
        body_vortex_panels;
        rhoinf=rhoinf[1],
        Vref=Vref[1],
    )

    ### --- TOTAL OUTPUTS --- ###

    # - Total Thrust - #
    total_thrust = sum([rotor_inviscid_thrust'; rotor_viscous_thrust'; body_thrust])

    # - Total Torque - #
    total_torque = sum([rotor_inviscid_torque; rotor_viscous_torque])

    # - Total Power - #
    total_power = sum([rotor_inviscid_power; rotor_viscous_power])

    # - Total Efficiency - #
    total_efficiency = get_total_efficiency(total_thrust, total_power, Vinf[1])
    ideal_efficiency = get_ideal_efficiency(total_thrust, rhoinf[1], Vinf[1], Rref[1])

    # - Body Induced Efficiency - #
    induced_efficiency = [
        get_induced_efficiency(
            rotor_inviscid_thrust[ir], body_thrust, rotor_inviscid_power[ir], Vinf[1]
        ) for ir in 1:nrotor
    ]

    # - Total Thrust and Torque Coefficients - #
    total_CT, total_CQ, total_CP = tqpcoeff(
        total_thrust, total_torque, total_power, rhoinf[1], Omega, Rref[1]
    )

    outs = (;
        # - States - #
        # TODO: instead of states, call these something else so that the outputs are generalized
        # - Wake Values - #
        wake=(; panel_strengths=gamw),
        # - Body Values - #
        bodies=(;
            # panel strengths
            panel_strengths=gamb[1:(idmaps.body_totnodes)],
            # body thrust
            total_thrust=sum(body_thrust),
            thrust_comp=body_thrust,
            induced_efficiency,
            # surface pressures
            cp_in,
            cp_out,
            cp_casing_in,
            cp_casing_out,
            casing_zpts,
            cp_nacelle_in,
            cp_nacelle_out,
            nacelle_zpts,
            cp_centerbody_in,
            cp_centerbody_out,
            centerbody_zpts,
            #individual body velocity contributions
            Vtot_in,
            Vtot_out,
            Vtot_prejump,
            vtot_body,
            vtot_jump,
            vtot_wake,
            vtot_rotors,
            Vtan_in,
            Vtan_out,
            vtan_casing_in,
            vtan_casing_out,
            vtan_nacelle_in,
            vtan_nacelle_out,
            vtan_centerbody_in,
            vtan_centerbody_out,
        ),
        # - Rotor Values - #
        rotors=(;
            circulation=Gamr,
            panel_strengths=sigr,
            efficiency=rotor_efficiency,
            inviscid_thrust=rotor_inviscid_thrust,
            inviscid_thrust_dist=rotor_inviscid_thrust_dist,
            viscous_thrust=rotor_viscous_thrust,
            viscous_thrust_dist=rotor_viscous_thrust_dist,
            thrust=rotor_thrust,
            CT=rotor_CT,
            # rotor torque
            inviscid_torque=rotor_inviscid_torque,
            inviscid_torque_dist=rotor_inviscid_torque_dist,
            viscous_torque=rotor_viscous_torque,
            viscous_torque_dist=rotor_viscous_torque_dist,
            torque=rotor_torque,
            CQ=rotor_CQ,
            # rotor power
            inviscid_power=rotor_inviscid_power,
            inviscid_power_dist=rotor_inviscid_power_dist,
            viscous_power=rotor_viscous_power,
            viscous_power_dist=rotor_viscous_power_dist,
            power=rotor_power,
            CP=rotor_CP,
            # - Blade Element Values - #
            cl,
            cd,
            alpha,
            beta1,
            blade_normal_force_per_unit_span,
            blade_tangential_force_per_unit_span,
        ),
        # - Total Values - #
        totals=(;
            thrust=total_thrust,
            torque=total_torque,
            power=total_power,
            CT=total_CT[1],
            CQ=total_CQ[1],
            CP=total_CP[1],
            total_efficiency=total_efficiency[1],
            ideal_efficiency,
        ),
        # - Intermediate Values from Residual - #
        intermediate_solve_values=(;
            vz_rotor,
            vtheta_rotor,
            Cm_wake,
            res_vals.reynolds,
            res_vals.mach,
            res_vals.Cz_rotor,
            res_vals.Ctheta_rotor,
            res_vals.Cmag_rotor,
            res_vals.Gamma_tilde,
            res_vals.H_tilde,
            res_vals.deltaGamma2,
            res_vals.deltaH,
            res_vals.vz_wake,
            res_vals.vr_wake,
            res_vals.Cm_avg,
        ),
    )

    if write_outputs
        write_data(
            outs,
            outfile;
            output_tuple_name=output_tuple_name,
            checkoutfileexists=checkoutfileexists,
        )
    end

    return outs
end

######################################################################
#                                                                    #
#                        Velocity Functions                          #
#                                                                    #
######################################################################

"""
"""
function get_body_tangential_velocities(
    gamb,
    gamw,
    sigr,
    ivb,
    Vinf,
    totnode,
    totpanel,
    nnode,
    npanel,
    tangent,
    controlpoints,
    endpanelidxs,
    wake_panel_ids_along_centerbody_wake_interface,
    wake_panel_ids_along_casing_wake_interface,
    centerbody_panel_ids_along_centerbody_wake_interface,
    duct_panel_ids_along_centerbody_wake_interface,
)

    # - setup - #
    nws, nrotor = size(sigr)
    TF = promote_type(eltype(gamb), eltype(gamw), eltype(sigr))
    (; v_bb, v_br, v_bw) = ivb

    # rename for convenience
    hwi = wake_panel_ids_along_centerbody_wake_interface
    dwi = wake_panel_ids_along_casing_wake_interface
    whi = centerbody_panel_ids_along_centerbody_wake_interface
    wdi = duct_panel_ids_along_centerbody_wake_interface

    # TODO also consider including the body wakes here as well.

    # - initialize total velocity - #
    Vtot = zeros(TF, 2, totpanel)

    # - Velocity Contributions from body - #

    vtot_body = similar(Vtot) .= 0
    for (i, vt) in enumerate(eachrow(vtot_body))
        vt .= @views v_bb[:, :, i] * gamb[1:size(v_bb, 2)]
    end
    Vtot .+= vtot_body

    # - Velocity Contributions from wake - #
    vtot_wake = similar(Vtot) .= 0
    for (i, vt) in enumerate(eachrow(vtot_wake))
        vt .= @views v_bw[:, :, i] * gamw
    end
    Vtot .+= vtot_wake # opposite sign from linear solve

    # - Velocity Contributions from rotors - #
    vtot_rotors = similar(Vtot) .= 0.0
    for jrotor in 1:nrotor
        rotorrange = (nws * (jrotor - 1) + 1):(nws * jrotor)
        for (i, vt) in enumerate(eachrow(vtot_rotors))
            vt .+= @views v_br[:, rotorrange, i] * sigr[rotorrange]
        end
    end
    Vtot .+= vtot_rotors # opposite sign from linear solve

    # - Influence from Freestream - #
    Vtot[1, :] .+= Vinf # opposite sign from linear solve
    Vtot_prejump = copy(Vtot)

    # - Add in Jump Term - #
    # duct
    jumpduct = (gamb[1:(npanel[1])] + gamb[2:(nnode[1])]) / 2

    # wake panels interfacing with duct
    jumpduct[wdi] .+= (gamw[dwi[1]:(dwi[end] - 1)] + gamw[(dwi[1] + 1):dwi[end]]) / 2.0

    # center body panels
    jumphub = (gamb[(nnode[1] + 1):(totnode - 1)] + gamb[(nnode[1] + 2):(totnode)]) / 2.0

    # wake panels interfacing with center body
    jumphub[whi] .+= (gamw[hwi[1]:(hwi[end] - 1)] + gamw[(hwi[1] + 1):hwi[end]]) / 2.0

    body_jump_term = [jumpduct; jumphub]

    vtot_jump = similar(Vtot) .= 0.0
    for (vt, tan) in zip(eachrow(vtot_jump), eachrow(tangent))
        vt .+= body_jump_term .* tan ./ 2.0
    end

    # assign velocities to each side of the panel
    Vtot_out = Vtot .- vtot_jump # outer side of boundary
    Vtot_in = Vtot .+ vtot_jump # inner side of boundary

    # Get the magnitude of the sum of the velocities and this is the surface velocity since the body velocities have been solved to eliminate the normal components in the summed velocities
    Vtan_out = sqrt.(Vtot_out[1, :] .^ 2 .+ Vtot_out[2, :] .^ 2)
    Vtan_in = sqrt.(Vtot_in[1, :] .^ 2 .+ Vtot_in[2, :] .^ 2)

    # - Split Velocities associates with inner and outer duct and hub - #
    # total tangential velocities
    vtan_casing_out, vtan_nacelle_out, vtan_centerbody_out, casing_zpts, nacelle_zpts, centerbody_zpts = split_bodies(
        Vtan_out, controlpoints, endpanelidxs
    )
    vtan_casing_in, vtan_nacelle_in, vtan_centerbody_in, casing_zpts, nacelle_zpts, centerbody_zpts = split_bodies(
        Vtan_in, controlpoints, endpanelidxs
    )

    return (;
        # Totals and Components:
        Vtan_in,
        Vtot_in,
        Vtan_out,
        Vtot_out,
        Vtot_prejump,
        vtot_body,
        vtot_jump,
        vtot_wake,
        vtot_rotors,
        # Splits:
        vtan_casing_in,
        vtan_casing_out,
        casing_zpts,
        vtan_nacelle_in,
        vtan_nacelle_out,
        nacelle_zpts,
        vtan_centerbody_in,
        vtan_centerbody_out,
        centerbody_zpts,
    )
end

"""
Calculate tangential velocity for a given net circulation and radial location
"""
function calculate_vtheta(Gamma_tilde, r)
    TF = promote_type(eltype(Gamma_tilde), eltype(r))
    vtheta = zeros(TF, length(r))

    for i in 1:length(r)
        if isapprox(r[i], 0.0)
            vtheta[i] = 0.0
        else
            vtheta[i] = Gamma_tilde ./ (2.0 * pi * r[i])
        end
    end

    return vtheta
end

"""
Calculate the induced velocities on one of the body wakes (unit velocity inputs determine which one)
"""
function calculate_induced_velocities_on_bodywake(
    vz_w, vr_w, gamw, vz_r, vr_r, sigr, vz_b, vr_b, gamb, Vinf
)

    # problem dimensions
    nrotor = size(sigr, 2) # number of rotors
    np = size(vz_b, 1) # number of panels in bodywake

    # initialize outputs
    vz = Vinf * ones(eltype(gamw), np) # axial induced velocity
    vr = zeros(eltype(gamw), np) # radial induced velocity

    # add body induced velocities
    @views vz[:] .+= vz_b * gamb
    @views vr[:] .+= vr_b * gamb

    # add wake induced velocities
    @views vz[:] .+= vz_w * gamw
    @views vr[:] .+= vr_w * gamw

    # add rotor induced velocities
    for jrotor in 1:nrotor
        @views vz[:] .+= vz_r[jrotor] * sigr[:, jrotor]
        @views vr[:] .+= vr_r[jrotor] * sigr[:, jrotor]
    end

    # return raw induced velocities
    return vz, vr
end

######################################################################
#                                                                    #
#                        Pressure Functions                          #
#                                                                    #
######################################################################

"""
Calculate steady pressure coefficient
"""
function steady_cp(Vs, Vinf, Vref)
    return (Vinf^2 .- Vs .^ 2) / Vref^2
end

"""
only used in post-process for cp.
expression not in dfdc theory, comes from source code.
"""
function calculate_entropy_jumps(sigr, Cz_rotor)
    # average sigr's
    sigr_avg = similar(Cz_rotor, size(sigr, 1) - 1, size(sigr, 2)) .= 0
    for (i, s) in enumerate(eachcol(sigr))
        sigr_avg[:, i] = (s[2:end] + s[1:(end - 1)]) / 2.0
    end

    # multiply by Cz_rotor's
    return sigr_avg .* Cz_rotor
end

"""
Calculate change in pressure coefficient aft of rotor, due to rotor
"""
function delta_cp(deltaH, deltaS, Vtheta, Vref)
    if isapprox(Vref, 0.0)
        return 0.0
    else
        return (2.0 * (deltaH - deltaS) .- Vtheta .^ 2) / Vref^2
    end
end

"""
Calculate net circulation and enthalpy and entropy disk jumps
"""
function calculate_rotor_jumps(Gamr, Omega, B, sigr, Cz_rotor)

    # - Calculate net circulations - #
    Gamma_tilde = calculate_net_circulation(Gamr, B)

    # - Calculate enthalpy disk jump - #
    Htilde = calculate_enthalpy_jumps(Gamr, Omega, B)

    # - Calculate entropy disk jump - #
    Stilde = calculate_entropy_jumps(sigr, Cz_rotor)

    return Gamma_tilde, Htilde, Stilde
end

"""
Calculate change in pressure coefficient due to rotors specifically on the body panels aft of the rotors
"""
function calculate_body_delta_cp!(cp, Gamr, sigr, Cz_rotor, Vref, Omega, B, cpr, didr, hidr)

    ## -- Calculate change in pressure coefficient -- ##

    Gamma_tilde, Htilde, Stilde = calculate_rotor_jumps(Gamr, Omega, B, sigr, Cz_rotor)

    nrotor = size(Gamr, 2)

    for irotor in 1:nrotor

        # - Get the tangential velocities on the bodies - #
        v_theta_duct = calculate_vtheta(
            Gamma_tilde[end, irotor], @view(cpr[1:didr[irotor]])
        )
        v_theta_hub = calculate_vtheta(Gamma_tilde[1, irotor], @view(cpr[hidr[irotor]:end]))

        # assemble change in cp due to enthalpy and entropy behind rotor(s)
        cp[1:didr[irotor]] += delta_cp(
            Htilde[end, irotor], Stilde[end, irotor], v_theta_duct, Vref
        )
        cp[hidr[irotor]:end] += delta_cp(
            Htilde[1, irotor], Stilde[1, irotor], v_theta_hub, Vref
        )
    end

    return nothing
end

"""
Calculate change in pressure coefficient due to rotors specifically on the body wakes
"""
function calculate_bodywake_delta_cp(Gamr, sigr, Cz_rotor, Vref, Omega, B, r; body="duct")

    ## -- Calculate change in pressure coefficient -- ##

    Gamma_tilde, Htilde, Stilde = calculate_rotor_jumps(Gamr, Omega, B, sigr, Cz_rotor)

    # - Get the tangential velocities on the bodies - #
    if body == "duct"
        gt = Gamma_tilde[end]
        ht = Htilde[end]
        st = Stilde[end]
    else
        gt = Gamma_tilde[1, end]
        ht = Htilde[1, end]
        st = Stilde[1, end]
    end

    v_theta_wake = calculate_vtheta(gt, r)

    # assemble change in cp due to enthalpy and entropy behind rotor(s)
    deltacp = delta_cp(ht, st, v_theta_wake, Vref)

    return deltacp
end

"""
calculate pressure coefficient distribution on duct/hub walls
formulation taken from DFDC source code. TODO: derive where the expressions came from.

"""
function get_body_cps(
    Vtan_in,
    Vtan_out,
    Gamr,
    sigr,
    Cz_rotor,
    Vinf,
    Vref,
    B,
    Omega,
    didr,
    hidr,
    controlpoints,
    endpanelidxs,
)

    # - Calculate standard pressure coefficient expression - #
    cp_in = steady_cp(Vtan_in, Vinf, Vref)
    cp_out = steady_cp(Vtan_out, Vinf, Vref)

    # - add the change in Cp on the walls due to enthalpy, entropy, and vtheta - #
    calculate_body_delta_cp!(
        cp_out, Gamr, sigr, Cz_rotor, Vref, Omega, B, @view(controlpoints[2, :]), didr, hidr
    )

    # - Split body strengths into inner/outer duct and hub - #
    cp_casing_in, cp_nacelle_in, cp_centerbody_in, casing_zpts, nacelle_zpts, centerbody_zpts = split_bodies(
        cp_in, controlpoints, endpanelidxs
    )
    cp_casing_out, cp_nacelle_out, cp_centerbody_out, _, _, _ = split_bodies(
        cp_out, controlpoints, endpanelidxs
    )

    return (;
        cp_in,
        cp_out,
        cp_casing_in,
        cp_nacelle_in,
        cp_centerbody_in,
        casing_zpts,
        nacelle_zpts,
        centerbody_zpts,
        cp_casing_out,
        cp_nacelle_out,
        cp_centerbody_out,
    )
end

"""
Calculate the pressure coefficient distributions on one of the body wakes
"""
function get_bodywake_cps(
    Gamr,
    vz_w,
    vr_w,
    gamw,
    vz_r,
    vr_r,
    sigr,
    vz_b,
    vr_b,
    gamb,
    panels,
    Cz_rotor,
    Omega,
    B,
    Vinf,
    Vref;
    body="duct",
)

    # - Get "surface" velocities - #

    # get induced velocities
    vz_bodywake, vr_bodywake = calculate_induced_velocities_on_bodywake(
        vz_w, vr_w, gamw, vz_r, vr_r, sigr, vz_b, vr_b, gamb, Vinf
    )

    # get "surface" velocities
    Vmat = [vz_bodywake vr_bodywake]
    vtan = [dot(v, t) for (v, t) in zip(eachrow(Vmat), panels.tangent)]

    # - Get steady pressure coefficients - #
    cp_steady = steady_cp(vtan, Vinf, Vref)

    # - Get delta cp - #
    deltacp = calculate_bodywake_delta_cp(
        Gamr, sigr, Cz_rotor, Vref, Omega, B, panels.controlpoint[2, :]; body=body
    )

    return cp_steady .+ deltacp, vtan
end

"""
Calculate dimensional and non-dimensional axial force on a single body
"""
function forces_from_pressure(cp_in, cp_out, panels; rhoinf=1.225, Vref=1.0)

    # - rename for convenience - #
    #just want x-component of normals since it's axisymmetric
    ns = panels.normal[1, :]
    #radial positions
    rs = panels.controlpoint[2, :]
    #panel lengths
    ds = panels.influence_length

    # - initialize - #
    cfx = zeros(eltype(cp_out), panels.nbodies) # axial force coefficient (all others are zero for axisymmetric case)

    # for each body
    for ib in 1:(panels.nbodies)
        # - rectangular integration due to constant panel strengths. - #
        for ip in panels.endpanelidxs[1, ib]:panels.endpanelidxs[2, ib]
            cfx[ib] += (cp_out[ip] - cp_in[ip]) * ns[ip] * ds[ip] * 2.0 * pi * rs[ip]
        end
    end

    #dimensionalize
    q = 0.5 * rhoinf * Vref^2

    #note, thrust is in negative x-direction
    return cfx .* q, cfx
end

"""
Calculate dimensional and non-dimensional axial force on a single body
"""
function forces_from_TEpanels!(
    thrust, force_coeff, cp_in, cp_out, panels; rhoinf=1.225, Vref=1.0
)

    #dimensionalize
    q = 0.5 * rhoinf * Vref^2

    for i in 1:(panels.nbodies)
        if panels.tenode[i, 2, 2] <= eps()
            # if it's the hub, don't average the first and last, but rather just the last
            cpi = cp_in[panels.endpanelidxs[i, 2]]
            cpo = cp_out[panels.endpanelidxs[i, 2]]
        else
            # if it's the duct, then average the first and last panel
            cpi =
                0.5 * (cp_in[panels.endpanelidxs[1, i]] + cp_in[panels.endpanelidxs[2, i]])
            cpo =
                0.5 *
                (cp_out[panels.endpanelidxs[1, i]] + cp_out[panels.endpanelidxs[2, i]])
        end

        r = 0.5 * sum(panels.tenode[i, :, 2])

        force_coeff[i] +=
            (cpo - cpi) *
            panels.tenormal[1, i] *
            panels.teinfluence_length[i] *
            2.0 *
            pi *
            r

        thrust[i] +=
            q *
            (cpo - cpi) *
            panels.tenormal[1, i] *
            panels.teinfluence_length[i] *
            2.0 *
            pi *
            r
    end

    return thrust, force_coeff
end

######################################################################
#                                                                    #
#                       Rotor Aero Performance                       #
#                                                                    #
######################################################################

function inviscid_rotor_thrust(Ctheta_rotor, Gamma_tilde, rotor_panel_length, rhoinf)

    # problem dimensions
    nr, nrotor = size(Gamma_tilde)

    # initialize
    dTi = similar(Gamma_tilde) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            # section thrust
            dTi[ir, irotor] =
                -rhoinf *
                Gamma_tilde[ir, irotor] *
                Ctheta_rotor[ir, irotor] *
                rotor_panel_length[ir, irotor]
        end
    end

    #sum the section thrust
    Tinv = sum(dTi; dims=1)

    return Tinv, dTi
end

function viscous_rotor_thrust(
    Cz_rotor, Cmag_rotor, B, chord, rotor_panel_length, cd, rhoinf
)

    # get dimensions
    nr, nrotor = size(Cz_rotor)

    #initialize
    dTv = similar(Cz_rotor) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            hrwc = 0.5 * rhoinf * Cmag_rotor[ir, irotor] * chord[ir, irotor]
            bdr = B[irotor] * rotor_panel_length[ir, irotor]
            dTv[ir, irotor] = -hrwc * cd[ir, irotor] * Cz_rotor[ir, irotor] * bdr
        end
    end

    Tvisc = sum(dTv; dims=1)

    return Tvisc, dTv
end

function inviscid_rotor_torque(
    Cz_rotor, rotor_panel_center, rotor_panel_length, Gamma_tilde, rhoinf
)

    # dimensions
    nr, nrotor = size(Gamma_tilde)

    # initialize
    dQi = similar(Gamma_tilde) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            rdr = rotor_panel_center[ir, irotor] * rotor_panel_length[ir, irotor]
            dQi[ir, irotor] = rhoinf * Gamma_tilde[ir, irotor] * Cz_rotor[ir, irotor] * rdr
        end
    end

    Qinv = sum(dQi; dims=1)

    return Qinv, dQi
end

function viscous_rotor_torque(
    Wtheta_rotor,
    Cmag_rotor,
    B,
    chord,
    rotor_panel_center,
    rotor_panel_length,
    cd,
    rhoinf;
    TF=nothing,
)

    # dimensions
    nr, nrotor = size(Wtheta_rotor)

    # initialize
    if isnothing(TF)
        TF = promote_type(
            eltype(Wtheta_rotor),
            eltype(Cmag_rotor),
            eltype(chord),
            eltype(rotor_panel_center),
            eltype(cd),
        )
    end

    dQv = zeros(TF, size(chord))

    for irotor in 1:nrotor
        for ir in 1:nr
            hrwc = 0.5 * rhoinf * Cmag_rotor[ir, irotor] * chord[ir, irotor]
            brdr =
                B[irotor] * rotor_panel_center[ir, irotor] * rotor_panel_length[ir, irotor]
            dQv[ir, irotor] = -hrwc * cd[ir, irotor] * Wtheta_rotor[ir, irotor] * brdr
        end
    end

    Qvisc = sum(dQv; dims=1)

    return Qvisc, dQv
end

function inviscid_rotor_power(Qinv, Omega)
    return Qinv .* Omega
end

function viscous_rotor_power(Qvisc, Omega)
    return Qvisc .* Omega
end

function get_total_efficiency(total_thrust, total_power, Vinf)
    TF = promote_type(eltype(total_thrust), eltype(total_power), eltype(Vinf))

    eta = zeros(TF, length(total_thrust))

    for i in 1:length(total_thrust)
        if Vinf <= 0.0 || total_power[i] < eps() || total_thrust[i] <= 0.0
            #do nothing, efficiency can't physically be negative or infinite.
        else
            eta[i] = total_thrust[i] * Vinf / total_power[i]
        end
    end

    return eta
end

function get_induced_efficiency(Tinv, Tduct, Pinv, Vinf)
    if Vinf <= 0.0 || Pinv <= 0.0
        return 0.0
    else
        return Vinf * (Tinv .+ Tduct) ./ Pinv
    end
end

function get_ideal_efficiency(total_thrust, rhoinf, Vinf, Rref)
    if Vinf != 0.0
        TC = total_thrust / (0.5 * rhoinf * Vinf^2 * pi * Rref^2)
        return 2.0 / (1.0 + sqrt(max(TC, -1.0) + 1.0))
    else
        return 0.0
    end
end

function tqpcoeff(thrust, torque, power, rhoinf, Omega, Rref)
    T = promote_type(eltype(thrust), eltype(torque), eltype(Omega))
    CT = zeros(T, length(Omega))
    CQ = zeros(T, length(Omega))
    CP = zeros(T, length(Omega))

    for (i, o) in enumerate(Omega)
        if isapprox(o, 0.0)
            CT[i] = CQ[i] = CP[i] = 0.0
        else
            # reference diameter
            D = 2.0 * Rref

            # rototion in rev per second
            n = o / (2.0 * pi)

            # thrust coefficient
            CT[i] = thrust[i] / (rhoinf * n^2 * D^4)

            # torque coefficient
            CQ[i] = torque[i] / (rhoinf * n^2 * D^5)

            # power coefficient
            CP[i] = power[i] / (rhoinf * n^3 * D^5)
        end
    end

    return CT, CQ, CP
end

function get_blade_loads(Wmag_rotor, beta1, cl, cd, chords, rhoinf)#, Rhub, Rtip, rotor_panel_centers,rotor_panel_lengths ,B, Omega)

    # dimensions
    nr, nrotor = size(Wmag_rotor)

    # initialize
    Np = similar(Wmag_rotor) .= 0.0
    Tp = similar(Wmag_rotor) .= 0.0

    for irotor in 1:nrotor
        for ir in 1:nr
            # rename for convenience
            cphi = cos(beta1[ir, irotor])
            sphi = sin(beta1[ir, irotor])

            # resolve lift and drag into normal and tangential coefficients
            cn = cl[ir, irotor] * cphi - cd[ir, irotor] * sphi
            ct = cl[ir, irotor] * sphi + cd[ir, irotor] * cphi

            # get the normal and tangential loads per unit length N' and T'
            Np[ir, irotor] =
                cn * 0.5 * rhoinf * Wmag_rotor[ir, irotor]^2 * chords[ir, irotor]
            Tp[ir, irotor] =
                ct * 0.5 * rhoinf * Wmag_rotor[ir, irotor]^2 * chords[ir, irotor]
        end
    end

    # Npfull = [zeros(nrotor)'; Np; zeros(nrotor)']
    # Tpfull = [zeros(nrotor)'; Tp; zeros(nrotor)']

    # TODO: consider comparing this with DFDC versions for them

    ## -- Integrate Loads to get Thrust and Torque
    # add hub/tip for complete integration.  loads go to zero at hub/tip.
    # rfull = [Rhub; rotor_panel_centers; Rtip]

    # thrust and torqe distributions
    # thrust = Npfull
    # torque = Tpfull .* rfull

    # integrate Thrust and Torque (trapezoidal)
    # T = B * fm.trapz(rfull, thrust)
    # Q = B * fm.trapz(rfull, torque)
    # - Actually use rectangle rather than trapezoid integration
    # T = B * sum(rotor_panel_lengths.*Np)
    # Q = B * sum(rotor_panel_lengths.* Tp.*rotor_panel_centers)
    # P = Q * Omega

    return Np, Tp
end
