"""
    post_process(
        solver_options,
        converged_states,
        prepost_containers,
        solve_container_caching,
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        operating_point,
        reference_parameters,
        A_bb_LU,
        airfoils,
        idmaps,
        problem_dimensions,
        multipoint_index;
        write_outputs=options.write_outputs,
        outfile=options.outfile,
        checkoutfileexists=options.checkoutfileexists,
        output_tuple_name=options.output_tuple_name,
        verbose=options.verbose,
    )

Post-process a converged nonlinear solve solution.

# Arguments
- `solver_options::SolverOptionsType` : A SolverOptionsType object (also used for dispatch)
- `converged_states::Vector{Float}` : the converged state variables
- `prepost_containers::NamedTuple` : the named tuple containing pre-allocated containers for the pre- and post-processing intermediate calculations
- `solve_container_cache::NamedTuple` : the cache and dimensions for intermediate values in the residual calculation
- `solve_parameter_cache_vector::Vector{Float}` : the applicably typed cache vector for the solve parameters
- `solve_parameter_cache_dims::NamedTuple` : the dimensions of the solver parameters
- `operating_point::OperatingPoint` : the operating point being analyzed
- `reference_parameters::ReferenceParameters` : a ReferenceParameters object
- `A_bb_LU::LinearAlgebra.LU` : LinearAlgebra LU factorization of the LHS matrix
- `airfoils::Vector{AFType}` : A matrix of airfoil types associated with each of the blade elements
- `idmaps::NamedTuple` : A named tuple containing index mapping used in bookkeeping throughout solve and post-process
- `problem_dimensions::ProblemDimensions` : A ProblemDimensions object

# Keyword Arguments
- `multipoint_index::Vector{Int}` : a one-dimensional vector containing the index of which multipoint analysis operating point is being analyzed.
- `write_outputs=options.write_outputs::Vector{Bool}` : a vector with the same length as number of multipoints indicating if the outputs should be saved.
- `outfile=options.outfile::Vector{String}` : a vector of file paths/names for where outputs should be written
- `checkoutfileexists=options.checkoutfileexists::Bool` : a flag for whether existing files should be checked for or if blind overwriting is okay.
- `output_tuple_name=options.output_tuple_name::Vector{String}` : the variable name(s) of the named tuple of outputs to be written.
- `verbose::Bool=false` : flag to print verbose statements

# Returns
- `outs::NamedTuple` : A named tuple containing all the output values including
  - `wake=(; panel_strengths=gamw)`
  - `bodies`
    - `panel_strengths`
    - `total_thrust`
    - `thrust_comp`
    - `induced_efficiency`
    - `cp_in`
    - `cp_out`
    - `cp_casing_in`
    - `cp_casing_out`
    - `zpts.casing_zpts`
    - `cp_nacelle_in`
    - `cp_nacelle_out`
    - `zpts.nacelle_zpts`
    - `cp_centerbody_in`
    - `cp_centerbody_out`
    - `zpts.centerbody_zpts`
    - `Vtot_in`
    - `Vtot_out`
    - `Vtot_prejump`
    - `vtot_body`
    - `vtot_jump`
    - `vtot_wake`
    - `vtot_rotors`
    - `Vtan_in`
    - `Vtan_out`
    - `vtan_casing_in`
    - `vtan_casing_out`
    - `vtan_nacelle_in`
    - `vtan_nacelle_out`
    - `vtan_centerbody_in`
    - `vtan_centerbody_out`
  - `rotors`
    - `circulation`
    - `panel_strengths`
    - `efficiency`
    - `inviscid_thrust`
    - `inviscid_thrust_dist`
    - `viscous_thrust`
    - `viscous_thrust_dist`
    - `thrust`
    - `CT`
    - `inviscid_torque`
    - `inviscid_torque_dist`
    - `viscous_torque`
    - `viscous_torque_dist`
    - `torque`
    - `CQ`
    - `inviscid_power`
    - `inviscid_power_dist`
    - `viscous_power`
    - `viscous_power_dist`
    - `power`
    - `CP`
    - `cl`
    - `cd`
    - `alpha`
    - `beta1`
    - `blade_normal_force_per_unit_span`
    - `blade_tangential_force_per_unit_span`
  - `totals`
    - `thrust`
    - `torque`
    - `power`
    - `CT`
    - `CQ`
    - `CP`
    - `total_efficiency`
    - `ideal_efficiency`
  - `intermediate_solve_values`
    - `vz_rotor`
    - `vtheta_rotor`
    - `Cm_wake`
    - `reynolds`
    - `mach`
    - `Cz_rotor`
    - `Ctheta_rotor`
    - `Cmag_rotor`
    - `Gamma_tilde`
    - `H_tilde`
    - `deltaGamma2`
    - `deltaH`
    - `vz_wake`
    - `vr_wake`
    - `Cm_avg`

"""
function post_process(
    solver_options,
    converged_states,
    prepost_containers,
    solve_container_caching,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
    operating_point,
    reference_parameters,
    A_bb_LU,
    airfoils,
    idmaps,
    problem_dimensions,
    multipoint_index;
    write_outputs=options.write_outputs,
    outfile=options.outfile,
    checkoutfileexists=options.checkoutfileexists,
    output_tuple_name=options.output_tuple_name,
    verbose=options.verbose,
)

    ### --- SETUP --- ###

    # - Extract Operating Point - #
    (; Vinf, asound, muinf, rhoinf, Omega) = operating_point

    # - Extract Reference Parameters - #
    (; Vref, Rref) = reference_parameters

    # - Extract PrePost Cache - #
    reset_containers!(prepost_containers; exception_keys=(:panels,:ivb))
    (;
        # stuff from pre-process
        panels,
        ivb,
        # rotor stuff
        rotor_inviscid_thrust,
        rotor_inviscid_thrust_dist,
        rotor_viscous_thrust,
        rotor_viscous_thrust_dist,
        rotor_thrust,
        rotor_inviscid_torque,
        rotor_inviscid_torque_dist,
        rotor_viscous_torque,
        rotor_viscous_torque_dist,
        rotor_torque,
        rotor_inviscid_power,
        rotor_inviscid_power_dist,
        rotor_viscous_power,
        rotor_viscous_power_dist,
        rotor_power,
        rotor_CT,
        rotor_CQ,
        rotor_CP,
        rotor_efficiency,
        induced_efficiency,
        blade_normal_force_per_unit_span,
        blade_tangential_force_per_unit_span,
        blade_loading_intermediate_containers,
        # body stuff
        zpts,
        vtan_tuple,
        cp_tuple,
        body_thrust,
        body_force_coefficient,
        # cp_tuple,
        # totals stuff
        total_thrust,
        total_torque,
        total_power,
        total_efficiency,
        ideal_efficiency,
        total_CT,
        total_CQ,
        total_CP,
    ) = prepost_containers

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
    # for combination solvers, get solver type that was last run
    if typeof(solver_options) <: CompositeSolverOptions
        idx = findlast(x -> x, (p -> p.converged[1]).(solver_options.solvers))
        sopt = solver_options.solvers[isnothing(idx) ? length(solver_options.solvers) : idx]
    elseif typeof(solver_options) <: ChainSolverOptions
        idx = findfirst(x -> x, (p -> p.converged[1]).(solver_options.solvers))
        sopt = solver_options.solvers[isnothing(idx) ? length(solver_options.solvers) : idx]
    else
        sopt = solver_options
    end

    res_vals = run_residual!(
        sopt,
        converged_states,
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
        multipoint_index
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
    rotor_panel_centers = blade_elements.rotor_panel_centers
    B = blade_elements.B
    chords = blade_elements.chords

    # - Rotor Thrust - #
    # inviscid thrust
    inviscid_rotor_thrust!(
        rotor_inviscid_thrust,
        rotor_inviscid_thrust_dist,
        Ctheta_rotor,
        Gamma_tilde,
        rotor_panel_lengths,
        rhoinf[1],
    )

    # viscous thrust
    viscous_rotor_thrust!(
        rotor_viscous_thrust,
        rotor_viscous_thrust_dist,
        Cz_rotor,
        Cmag_rotor,
        B,
        chords,
        rotor_panel_lengths,
        cd,
        rhoinf[1],
    )

    # total thrust
    rotor_thrust .= rotor_inviscid_thrust .+ rotor_viscous_thrust

    # - Rotor Torque - #
    # inviscid torque
    inviscid_rotor_torque!(
        rotor_inviscid_torque,
        rotor_inviscid_torque_dist,
        Cz_rotor,
        rotor_panel_centers,
        rotor_panel_lengths,
        Gamma_tilde,
        rhoinf[1],
    )

    # viscous torque
    viscous_rotor_torque!(
        rotor_viscous_torque,
        rotor_viscous_torque_dist,
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
    rotor_torque .= rotor_inviscid_torque .+ rotor_viscous_torque

    # - Rotor Power - #
    # inviscid power
    rotor_power!(
        rotor_inviscid_power,
        rotor_inviscid_power_dist,
        rotor_inviscid_torque,
        rotor_inviscid_torque_dist,
        Omega,
    )

    # viscous power
    rotor_power!(
        rotor_viscous_power,
        rotor_viscous_power_dist,
        rotor_viscous_torque,
        rotor_viscous_torque_dist,
        Omega,
    )

    # total power
    rotor_power .= rotor_inviscid_power .+ rotor_viscous_power

    # - Rotor Performance Coefficients - #
    tqpcoeff!(
        rotor_CT,
        rotor_CQ,
        rotor_CP,
        rotor_thrust,
        rotor_torque,
        rotor_power,
        rhoinf[1],
        Omega,
        Rref[1],
    )

    # - Rotor Efficiency - #
    get_total_efficiency!(rotor_efficiency, rotor_thrust, rotor_power, Vinf[1])

    # - Blade Loading - #
    get_blade_loads!(
        blade_normal_force_per_unit_span,
        blade_tangential_force_per_unit_span,
        Cmag_rotor,
        beta1,
        cl,
        cd,
        chords,
        rhoinf[1],
        blade_loading_intermediate_containers,
    )

    ### --- BODY OUTPUTS --- ###
    # - Surface Velocity on Bodies - #
    # TODO: update any tests for this function
    get_body_tangential_velocities!(
        vtan_tuple,
        gamb,
        gamw,
        sigr,
        ivb,
        Vinf[],
        Int(body_vortex_panels.totnode[]),
        Int(body_vortex_panels.totpanel[]),
        Int.(body_vortex_panels.nnode),
        Int.(body_vortex_panels.npanel),
        body_vortex_panels.tangent,
        body_vortex_panels.controlpoint,
        Int.(body_vortex_panels.endpanelidxs),
        idmaps.wake_node_ids_along_centerbody_wake_interface,
        idmaps.wake_node_ids_along_casing_wake_interface,
        idmaps.centerbody_panel_ids_along_centerbody_wake_interface,
        idmaps.duct_panel_ids_along_centerbody_wake_interface,
        zpts,
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
        vtan_nacelle_in,     # tangent velocity along nacelle inside of duct body
        vtan_nacelle_out,    # tangent velocity along nacelle outside of duct body
        vtan_centerbody_in,  # tangent velocity inside of centerbody
        vtan_centerbody_out, # tangent velocity outside of centerbody
    ) = vtan_tuple

    # - Surface Pressure on Bodies - #
    get_body_cps!(
        cp_tuple,
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
        zpts,
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
    forces_from_pressure!(
        body_thrust,
        body_force_coefficient,
        cp_in,
        cp_out,
        body_vortex_panels;
        rhoinf=rhoinf[1],
        Vref=Vref[1],
    )

    # add thrust from trailing edge panels on bodies
    forces_from_TEpanels!(
        body_thrust,
        body_force_coefficient,
        cp_in,
        cp_out,
        body_vortex_panels;
        rhoinf=rhoinf[1],
        Vref=Vref[1],
    )

    ### --- TOTAL OUTPUTS --- ###

    # - Total Thrust - #
    total_thrust[] = sum([rotor_inviscid_thrust'; rotor_viscous_thrust'; body_thrust])

    # - Total Torque - #
    total_torque[] = sum([rotor_inviscid_torque; rotor_viscous_torque])

    # - Total Power - #
    total_power[] = sum([rotor_inviscid_power; rotor_viscous_power])

    # - Total Efficiency - #
    get_total_efficiency!(total_efficiency, total_thrust, total_power, Vinf[])
    ideal_efficiency[] = get_ideal_efficiency(total_thrust[], rhoinf[], Vinf[], Rref[])

    # - Body Induced Efficiency - #
    get_induced_efficiency!(
        induced_efficiency,
        rotor_inviscid_thrust,
        sum(body_thrust),
        rotor_inviscid_power,
        Vinf[],
    )

    # - Total Thrust and Torque Coefficients - #
    tqpcoeff!(
        total_CT,
        total_CQ,
        total_CP,
        total_thrust,
        total_torque,
        total_power,
        rhoinf[1],
        Omega,
        Rref[1],
    )

    outs = (;
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
            zpts.casing_zpts,
            cp_nacelle_in,
            cp_nacelle_out,
            zpts.nacelle_zpts,
            cp_centerbody_in,
            cp_centerbody_out,
            zpts.centerbody_zpts,
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
#                          RESIDUAL FUNCTIONS                        #
#                                                                    #
######################################################################

"""
    run_residual!(
        solver_options::SolverOptionsType,
        converged_states,
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
        multipoint_index
    )

Run through the residual function post-convergence to save needed intermediate values for the rest of post-processing.

# Arguments
- `solver_options::SolverOptionsType` : A SolverOptionsType object (also used for dispatch)
- `converged_states::Vector{Float}` : the converged state variables
- `state_dims::NamedTuple` : a named tuple containing the sizes of the state variables
- `solve_container_cache::PreallocationTools.DiffCache` : the cache for intermediate values in the residual calculation
- `solve_container_cache_dims::NamedTuple` : the dimensions of the solve container cache
- `operating_point::OperatingPoint` : the operating point being analyzed
- `ivr::NamedTuple` : A named tuple containing arrays of induced velocities on the rotors
- `ivw::NamedTuple` : A named tuple containing arrays of induced velocities on the wake
- `linsys::NamedTuple` : A named tuple containing cacheable data for the linear system, including:
  - `A_bb::Array{Float}` : AIC (LHS) matrix for the panel method system
  - `b_bf::Array{Float}` : Initial system RHS vector based on freestrem magnitude
  - `A_br::Array{Float}` : Unit normal velocity from rotors onto body panels
  - `A_pr::Array{Float}` : Unit normal velocity from rotors onto body internal psuedo control points
  - `A_bw::Array{Float}` : Unit normal velocity from wake onto body panels
  - `A_pw::Array{Float}` : Unit normal velocity from wake onto body internal psuedo control points
  - `A_bb_LU::LinearAlgebra.LU` : LinearAlgebra LU factorization of the LHS matrix
- `blade_elements::NamedTuple` : A named tuple containing blade element information
- `wakeK::Matrix{Float}` : A matrix of precomputed geometric constants used in the calculation of the wake vortex strengths
- `idmaps::NamedTuple` : A named tuple containing index mapping used in bookkeeping throughout solve and post-process
- `multipoint_index::Vector{Int}` : a one-dimensional vector containing the index of which multipoint analysis operating point is being analyzed.

# Returns
- `res_vals::NamedTuple` : A named tuple containing the state variables and populated solve containers.
"""
function run_residual!(
solver_options::TS,
converged_states,
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
multipoint_index
) where {TS<:ExternalSolverOptions}

#=
  NOTE: we want to get all the intermediate values available to user if desired.
  The solve_containers cache will contain all the intermediate values after running the estimate states function.
=#
# - Separate out the state variables - #
vz_rotor, vtheta_rotor, Cm_wake = extract_state_variables(
    solver_options, converged_states, state_dims
)

# - Extract and Reset Cache - #
# get cache vector of correct types
solve_container_cache_vec = @views PreallocationTools.get_tmp(
    solve_container_cache, converged_states
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
converged_states,
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
multipoint_index
)

#=
  NOTE: we want to get all the intermediate values available to user if desired.
  The solve_containers cache will contain all the intermediate values after running the insides of the residual function
=#
# - Separate out the state variables - #
Gamr, sigr, gamw = extract_state_variables(solver_options, converged_states, state_dims)

# - Extract and Reset Cache - #
# get cache vector of correct types
solve_container_cache_vec = @views PreallocationTools.get_tmp(
    solve_container_cache, converged_states
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
    idmaps,
    multipoint_index;
    verbose=false,
    )

    return (; Gamr, sigr, gamw, solve_containers...)
end
