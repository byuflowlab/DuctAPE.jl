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
        boundary_layer_options,
        A_bb_LU,
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
- `BoundaryLayerOptions::BoundaryLayerOptions` : a BoundaryLayerOptions object
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
`outs::NamedTuple` : A named tuple containing all the output values including

# Extended help

Full Outputs Information:
Each for any level of these outputs with sub-levels, the higher level is a NamedTuple comprised of the of items listed below it. For example:
`outs.bodies.boundary_layers.stagnation_indices`

- `bodies::NamedTuple` : NamedTuple containing outputs related to the duct and center body
  - `panel_strengths::Vector{Float}` : body vortex panel strengths
  - `inviscid_thrust::Vector{Float}` : dimensional force in positive axial direction for duct and center body.
  - `body_force_coefficient::Vector{Float}` : force coefficients associated with dimensional inviscid thrust components
  - `viscous_drag::Vector{Float}` : dimensional force in negative axial direction for duct and center body (note: zero if boundary layer is turned off, and center body is always zero since no drag method is yet implemented for the center body)
  - `thrust_comp::Vector{Float}` : `inviscid_thrust` .- `viscous_drag`
  - `total_thrust::Float` : sum(`thrust_comp`)
  - `induced_efficiency::Vector{Float}` : body-induced propulsive efficiency
  - `cp_in::Vector{Float}` : inside pressure distribution for all bodies
  - `cp_out::Vector{Float}` : inside pressure distribution for all bodies
  - `cp_casing_in::Vector{Float}` : inside pressure distribution for duct casing
  - `cp_casing_out::Vector{Float}` : outside pressure distribution for duct casing
  - `casing_zpts::Vector{Float}` : axial component of casing control points
  - `cp_nacelle_in::Vector{Float}` : inside pressure distribution for duct nacelle
  - `cp_nacelle_out::Vector{Float}` : inside pressure distribution for duct nacelle
  - `nacelle_zpts::Vector{Float}` : axial component of nacelle control points
  - `cp_center_body_in::Vector{Float}` : inside pressure distribution for center bodies
  - `cp_center_body_out::Vector{Float}` : inside pressure distribution for center bodies
  - `center_body_zpts::Vector{Float}` : axial components of center body control points
  - `Vtot_in::Matrix{Float}` : total inner surface velocity distribution for all bodies. row 1 is vz, row 2 is vr, columns are control points.
  - `Vtot_out::Matrix{Float}` : total outer surface velocity distribution for all bodies. row 1 is vz, row 2 is vr, columns are control points.
  - `Vtot_prejump::Matrix{Float}` : total surface velocity distribution before velocity jumps terms are applied for all bodies. row 1 is vz, row 2 is vr, columns are control points.
  - `vtot_body::Vector{Float}` : body-induced velocity on the body surfaces
  - `vtot_jump::Vector{Float}` : velocity due to jump terms in Fredholm equation
  - `vtot_wake::Vector{Float}` : wake-induced velocity on the body surfaces
  - `vtot_rotors::Vector{Float}` : rotor-induced velocity on the body surfaces
  - `Vtan_in::Vector{Float}` : inner surface tangential velocity distribution for all bodies
  - `Vtan_out::Vector{Float}` : outer surface tangential velocity distribution for all bodies
  - `vtan_casing_in::Vector{Float}` : inner surface tangential velocity distribution for duct casing
  - `vtan_casing_out::Vector{Float}` : outer surface tangential velocity distribution for duct casing
  - `vtan_nacelle_in::Vector{Float}` : inner surface tangential velocity distribution for duct nacelle
  - `vtan_nacelle_out::Vector{Float}` : outer surface tangential velocity distribution for duct nacelle
  - `vtan_center_body_in::Vector{Float}` : inner surface tangential velocity distribution for center body
  - `vtan_center_body_out::Vector{Float}` : outer surface tangential velocity distribution for center bodies
  - `boundary_layers::NamedTuple` : NamedTuple containing information from the boundary layer solve (if `model_drag` in `boundary_layer_options` was set to true).
    - `stagnation_indices::Vector{Int}` : indices surrounding stagnation point
    - `upper_initial_states::Vector{Float}` : upper side initial states
    - `upper_solved_states::Matrix{Float}` : upper side solved states
    - `upper_solved_steps::Vector{Float}` : final steps associated with upper side solved states
    - `lower_initial_states::Vector{Float}` : lower side initial states
    - `lower_solved_states::Matrix{Float}` : lower side solved states
    - `lower_solved_steps::Vector{Float}` : final steps associated with lower side solved states
    - `surface_length_upper::Vector{Float}` : cumulative panel lengths on upper side
    - `surface_length_lower::Vector{Float}` : cumulative panel lengths on lower side
    - `stag_point::Float` : curve length at which stagnation point is located
    - `split_ratio::Float` : ratio of lower to total surface length
    - `separation_point_ratio_upper::Float` : ratio of upper side separation point location relative to upper side surface length
    - `separation_point_ratio_lower::Float` : ratio of lower side separation point location relative to lower side surface length
    - `cdc_upper::Float` : drag coefficient times chord length for upper side
    - `cdc_lower::Float` : drag coefficient times chord length for lower side
    - `vtdotpv::Vector{Float}` : dot product of tangential velocity and panel vector for duct
    - `boundary_layer_functions_lower::NamedTuple` : NamedTuple of functions generated for use in boundary layer solution. For Head's method these are:
      - `edge_velocity::FLOWMath.Akima` : spline of surface velocity relative to surface length
      - `edge_acceleration::FLOWMath.Akima` : spline of derivatives of `edge_velocity` relative to surface length
      - `edge_density::FLOWMath.Akima` : spline of density relative to surface length (constant for Head's method)
      - `edge_viscosity::FLOWMath.Akima` : spline of viscosity relative to surface length (constant for Head's method)
    - `boundary_layer_functions_upper::NamedTuple` : same as `boundary_layer_functions_lower` but for upper side.
- `rotors::NamedTuple` : NamedTuple of items related to rotor(s)
  - `circulation::Matrix{Float}` : blade element circulation values
  - `panel_strengths::Matrix{Float}` : balde source panel strengths
  - `efficiency::Vector{Float}` : rotor efficiency
  - `inviscid_thrust::Vector{Float}` : inviscid componenets of rotor thrust
  - `inviscid_thrust_dist::Matrix{Float}` : inviscid thrust component for each blade element
  - `viscous_thrust::Vector{Float}` : viscous componenets of rotor thrust
  - `viscous_thrust_dist::Martix{Float}` : viscous trhust component for each blade element
  - `thrust::Vector{Float}` : total rotor thrusts
  - `CT::Vector{Float}` : rotor thrust coefficients
  - `inviscid_torque::Vector{Float}` : inviscid components of rotor torque
  - `inviscid_torque_dist::Matrix{Float}` :inviscid torque component for each blade element
  - `viscous_torque::Vector{Float}` : viscous components of rotor torque
  - `viscous_torque_dist::Matrix{Float}` : viscous torque component for each blade element
  - `torque::Vector{Float}` : total rotor torques
  - `CQ::Vector{Float}` : rotor torque coefficients
  - `inviscid_power::Vector{Float}` : inviscid components of rotor power
  - `inviscid_power_dist::Matrix{Float}` : inviscid power component for each blade element
  - `viscous_power::Vector{Float}` : viscous components of rotor power
  - `viscous_power_dist::Matrix{Float}` : viscous power component for each blade element
  - `power::Vector{Float}` : total rotor powers
  - `CP::Vector{Float}` : rotor power coefficients
  - `cl::Matrix{Float}` : lift coefficient values for each blade element
  - `cd::Matrix{Float}` : drag coefficient values for each blade element
  - `alpha::Matrix{Float}` : angle of attack values for each blade element
  - `beta1::Matrix{Float}` : inflow angle values for each blade element
  - `blade_normal_force_per_unit_span::Matrix{Float}` : normal force per unit span values for each blade element
  - `blade_tangential_force_per_unit_span::Matrix{Float}` : tangential force per unit span values for each blade element
- `wake::NamedTuple` : NamedTuple containing items related to the wake
  - `panel_strengths::Vector{Float}` : wake vortex panel strengths
- `totals::NamedTuple` : NamedTuple containing total system values
  - `thrust::Vector{Float}` : total system thrust
  - `torque::Vector{Float}` : total system torque
  - `power::Vector{Float}` : total system power
  - `CT::Float` : total system thrust coefficient
  - `CQ::Float` : total system torque coefficient
  - `CP::Float` : total system power coefficient
  - `total_efficiency::Vector{Float}` : total propulsive efficiency
  - `ideal_efficiency::Vector{Float}` : ideal propulsive efficiency
- `intermediate_solve_values::NamedTuple` : NamedTuple containing items used inside the solver at their converged values.
  - `vz_rotor::Matrix{Float}` : axial velocity induced on rotor blade elements
  - `vtheta_rotor::Matrix{Float}` : swirl velocity induced on rotor blade elements
  - `reynolds::Matrix{Float}` : Reynolds numbers at each blade element
  - `mach::Matrix{Float}` : Mach numbers at each blade element
  - `Cz_rotor::Matrix{Float}` : absolute frame axial velocity at rotor blade elements
  - `Ctheta_rotor::Matrix{Float}` : absolute frame swirl velocity at rotor blade elements
  - `Cmag_rotor::Matrix{Float}` : absolute frame meridional velocity at rotor blade elements
  - `Gamma_tilde::Matrix{Float}` : net circulation of upstream and current blade elements
  - `H_tilde::Matrix{Float}` : net enthalpy of upstream and current blade elements
  - `deltaGamma2::Matrix{Float}` : squared circulation changes between adjacent blade elements
  - `deltaH::Matrix{Float}` : enthalpy changes between adjacent blade elements
  - `vz_wake::Vector{Float}` : axial velocity induced on wake control points
  - `vr_wake::Vector{Float}` : radial velocity induced on wake control oints
  - `Cm_wake::Vector{Float}` : absolute frame meridional velocity at wake control points
  - `Cm_avg::Vector{Float}` : absolute frame meridional velocity at wake panel nodes
- `reference_values::NamedTuple` : NamedTuple containing items used in computing coefficient values
  - `Vinf::Float` : Freestream velocity used in coefficient definitions
  - `Vref::Float` : Reference velocity used in coefficient definitions
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
    boundary_layer_options,
    A_bb_LU,
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
    reset_containers!(prepost_containers; exception_keys=(:panels, :ivb))
    (;
        # Stuff from Pre-process
        panels,
        ivb,
        # Rotor Stuff
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
        # Body Stuff
        zpts,
        vtan_tuple,
        cp_tuple,
        body_inviscid_thrust,
        body_viscous_drag,
        body_thrust,
        body_force_coefficient,
        # Totals Stuff
        total_thrust,
        total_torque,
        total_power,
        total_efficiency,
        ideal_efficiency,
        total_CT,
        total_CQ,
        total_CP,
        # Boundary Layer Stuff #TODO: add these to caches at some point (requires re-work of boundary layer implementation likey)
        # duct_viscous_drag,
        # boundary_layer_outputs,
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
    blade_elements = blade_elements
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
        multipoint_index,
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

    # - rename for convenience - #
    (; nbe, nrotor) = problem_dimensions
    (;rotor_panel_centers, B, chords, is_stator) = blade_elements

    ### --- BODY OUTPUTS --- ###
    # - Surface Velocity on Bodies - #
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
        idmaps.wake_node_ids_along_center_body_wake_interface,
        idmaps.wake_node_ids_along_casing_wake_interface,
        idmaps.center_body_panel_ids_along_center_body_wake_interface,
        idmaps.duct_panel_ids_along_center_body_wake_interface,
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
        vtan_center_body_in,  # tangent velocity inside of center_body
        vtan_center_body_out, # tangent velocity outside of center_body
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
        idmaps.id_of_first_center_body_panel_aft_of_each_rotor,
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
        cp_center_body_in,  # surface pressure along inside of center_body
        cp_casing_out,     # surface pressure along outside of casing
        cp_nacelle_out,    # surface pressure along outside of nacelle
        cp_center_body_out, # surface pressure along outside of center_body
    ) = cp_tuple

    # - Calculate Thrust from Bodies - #
    forces_from_pressure!(
        body_inviscid_thrust,
        body_force_coefficient,
        cp_in,
        cp_out,
        body_vortex_panels;
        rhoinf=rhoinf[1],
        Vref=Vref[1],
    )

    # add thrust from trailing edge panels on bodies
    forces_from_TEpanels!(
        body_inviscid_thrust,
        body_force_coefficient,
        cp_in,
        cp_out,
        body_vortex_panels;
        rhoinf=rhoinf[1],
        Vref=Vref[1],
    )

    # - Duct Viscous Drag - #

    if boundary_layer_options.model_drag

        #TODO; make this in place once you've finalized outputs
        body_viscous_drag[1], boundary_layer_outputs = compute_viscous_drag_duct(
            boundary_layer_options,
            Vtan_out[1:Int(body_vortex_panels.npanel[1])],
            Vtot_out[:, 1:Int(body_vortex_panels.npanel[1])],
            body_vortex_panels.controlpoint[:, 1:Int(body_vortex_panels.npanel[1])],
            body_vortex_panels.influence_length[1:Int(body_vortex_panels.npanel[1])],
            body_vortex_panels.tangent[:, 1:Int(body_vortex_panels.npanel[1])],
            # body_vortex_panels.node[2, Int(body_vortex_panels.nnode[1])],
            Rref[1],
            operating_point,
            reference_parameters;
            verbose=boundary_layer_options.verbose,
        )

    else
        body_viscous_drag = [0.0, 0.0]
        boundary_layer_outputs = nothing
    end

    # Total body thrust
    body_thrust = sum(body_inviscid_thrust) - sum(body_viscous_drag)

    ### ----- ROTOR OUTPUTS ----- ###

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

    # Apply penalty to rotor performance
    if boundary_layer_options.model_drag &&
        boundary_layer_options.apply_separation_penalty_to_rotor &&
        boundary_layer_options.separation_penalty_lower > eps()
        rotor_penalty = separation_penalty(
            boundary_layer_outputs.s_sep_lower,
            boundary_layer_outputs.lower_steps,
            boundary_layer_options.separation_allowance_lower,
            boundary_layer_options.separation_penalty_lower,
        )

        rotor_inviscid_torque .*= (1.0 .- rotor_penalty)
        rotor_viscous_torque .*= (1.0 .- rotor_penalty)
        rotor_inviscid_power .*= (1.0 .- rotor_penalty)
        rotor_viscous_power .*= (1.0 .- rotor_penalty)
        rotor_inviscid_thrust .*= (1.0 .- rotor_penalty)
        rotor_viscous_thrust .*= (1.0 .- rotor_penalty)
    end

    # - Rotor Totals - #

    # total thrust
    rotor_thrust .= rotor_inviscid_thrust .+ rotor_viscous_thrust

    # total torque
    rotor_torque .= rotor_inviscid_torque .+ rotor_viscous_torque

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
        is_stator,
        blade_loading_intermediate_containers,
    )

    ### --- TOTAL OUTPUTS --- ###

    # - Total Thrust - #
    total_thrust[] = sum([rotor_inviscid_thrust; rotor_viscous_thrust; body_thrust])

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
        Omega[1],
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
            body_force_coefficient=body_force_coefficient,
            inviscid_thrust=body_inviscid_thrust,
            viscous_drag=body_viscous_drag,
            thrust_comp=body_inviscid_thrust .- body_viscous_drag,
            total_thrust=sum(body_thrust),
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
            cp_center_body_in,
            cp_center_body_out,
            zpts.center_body_zpts,
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
            vtan_center_body_in,
            vtan_center_body_out,
            # boundary layers
            boundary_layers=boundary_layer_outputs,
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
            total_efficiency=total_efficiency,
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
        reference_values=(; Vinf=operating_point.Vinf[], Vref=reference_parameters.Vref[]),
    )

    if write_outputs
        write_data(
            outs,
            outfile;
            output_tuple_name=output_tuple_name,
            checkoutfileexists=checkoutfileexists,
        )
    end

    return deepcopy(outs)
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
    multipoint_index,
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
    multipoint_index,
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
        zeros(eltype(converged_states), 2),
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

"""
"""
function run_residual!(
    solver_options::ModCSORSolverOptions,
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
    multipoint_index,
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
    estimate_CSOR_states!(
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
