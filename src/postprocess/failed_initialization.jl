"""
    failed_initialization(
        ducted_rotor,
        operating_point,
        reference_parameters,
        prepost_containers,
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        A_bb_LU,
        idmaps,
        problem_dimensions,
        options;
        solve_container_caching=nothing,
    )

Constructs a valid, but non-converged output solution for a ducted rotor system. This is typically used when solver convergence fails but a valid output structure (with zeros or initial guesses) is still required.

# Arguments
- `ducted_rotor` : Problem definition and geometry of the ducted rotor.
- `operating_point` : Flow properties and operating conditions (e.g., velocity, density, RPM).
- `reference_parameters` : Non-dimensionalization reference values (e.g., Vref, Rref).
- `prepost_containers` : Storage containers for results and intermediate variables.
- `solve_parameter_cache_vector` : Flat cache vector containing solver-specific data.
- `solve_parameter_cache_dims` : Dimensional metadata used to interpret the solve parameter cache.
- `A_bb_LU` : LU-decomposed matrix of body-body influence coefficients.
- `idmaps` : Index maps for system matrices and solution vectors.
- `problem_dimensions` : Struct containing problem size information.
- `options` : Solver and physical modeling options.
- `solve_container_caching=nothing` : Optional preallocated cache of solve container objects.

# Returns
- `outs::NamedTuple` : Named tuple of simulation outputs with the following fields:
    - `wake` : Tuple of wake panel strengths.
    - `bodies` : Dictionary of body-related quantities, including panel strengths, thrusts, efficiencies, surface velocities/pressures, and boundary layer results (if enabled).
    - `rotors` : Rotor-related quantities including thrust, torque, power, efficiencies, blade element forces, and aerodynamic coefficients.
    - `totals` : Combined thrust, torque, power, and efficiency metrics.
    - `intermediate_solve_values` : Miscellaneous intermediate quantities computed by the solver (e.g., induced velocities, circulation deltas, etc.).
    - `reference_values` : Reference freestream and nondimensionalization values.

All quantities will be physically consistent and correctly shaped, but will not represent a converged flowfield.
"""
function failed_initialization(
    ducted_rotor,
    operating_point,
    reference_parameters,
    prepost_containers,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
    A_bb_LU,
    idmaps,
    problem_dimensions,
    options;
    solve_container_caching=nothing,
)

    ### --- SETUP --- ###

    # Set up Solve Container Cache
    if isnothing(solve_container_caching)
        solve_container_caching = allocate_solve_container_cache(
            options.solver_options, ducted_rotor.paneling_constants
        )
    end

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

    TF = eltype(total_CT)

    # - Extract Panels - #
    (; body_vortex_panels, rotor_source_panels, wake_vortex_panels) = panels

    # rename rotor panel lengths
    rotor_panel_lengths = rotor_source_panels.influence_length

    # - Extract Solve Parameter Cache - #
    solve_parameter_tuple = withdraw_solve_parameter_cache(
        options.solver_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )
    (; ivr, ivw, linsys, blade_elements, wakeK) = solve_parameter_tuple

    # put airfoils in blade elements and LU decomp into linsys
    blade_elements = blade_elements
    linsys = (; linsys..., A_bb_LU)

    # - Extract Solve Container Cache - #
    (; solve_container_cache, solve_container_cache_dims) = solve_container_caching

    # - initialize Residual to get intermediate values - #
    # for combination solvers, get solver type that was last run
    if typeof(options.solver_options) <: CompositeSolverOptions
        idx = findlast(x -> x, (p -> p.converged[1]).(options.solver_options.solvers))
        sopt = options.solver_options.solvers[if isnothing(idx)
            length(options.solver_options.solvers)
        else
            idx
        end]
    elseif typeof(options.solver_options) <: ChainSolverOptions
        idx = findfirst(x -> x, (p -> p.converged[1]).(options.solver_options.solvers))
        sopt = options.solver_options.solvers[if isnothing(idx)
            length(options.solver_options.solvers)
        else
            idx
        end]
    else
        sopt = options.solver_options
    end

    res_vals = return_failed_residual!(
        sopt,
        TF,
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
        options.multipoint_index,
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

    if options.boundary_layer_options.model_drag
        boundary_layer_outputs = return_failed_boundary_layer_outputs(
            options.boundary_layer_options,
            Vtan_out[1:Int(body_vortex_panels.npanel[1])],
            Vtot_out[:, 1:Int(body_vortex_panels.npanel[1])],
            body_vortex_panels.controlpoint[:, 1:Int(body_vortex_panels.npanel[1])],
            body_vortex_panels.influence_length[1:Int(body_vortex_panels.npanel[1])],
            body_vortex_panels.tangent[:, 1:Int(body_vortex_panels.npanel[1])],
            Rref[1],
            operating_point,
            reference_parameters,
        )

        body_viscous_drag = [0.0, 0.0]
    else
        body_viscous_drag = [0.0, 0.0]
        boundary_layer_outputs = nothing
    end

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
            thrust_comp=body_inviscid_thrust .+ body_viscous_drag,
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

    return outs
end

"""
    return_failed_residual!(
        solver_options::TS,
        TF,
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

Returns a tuple of default or fallback residual-related values for a failed solver call, specifically for use in `failed_initialization`. This ensures downstream postprocessing code can still execute by returning zeroed or nominal values consistent in shape and type.

# Arguments
- `solver_options::TS` : Specific solver configuration (must be a subtype of `ExternalSolverOptions`).
- `TF` : Floating-point type used throughout the solver (typically `Float64`).
- `state_dims` : Dimensions of the solver state vector (e.g., number of state variables per component).
- `solve_container_cache` : Preallocated container cache used to store reusable solver objects.
- `solve_container_cache_dims` : Metadata about the structure and dimensions of the solve container cache.
- `operating_point` : Flow conditions (e.g., freestream velocity, density, RPM).
- `ivr` : Rotor-related solver input variables.
- `ivw` : Wake-related solver input variables.
- `linsys` : Linear system data (e.g., LU decompositions, influence matrices).
- `blade_elements` : Blade geometry and airfoil data used in blade element calculations.
- `wakeK` : Wake kernel matrix used in unsteady wake modeling.
- `idmaps` : Identifier maps for connecting variables and system components.
- `multipoint_index` : Index to identify the current multipoint condition (for multipoint optimization).

# Returns
Named tuple with the following:
- `vz_rotor` : Axial induced velocities at rotor control points.
- `vtheta_rotor` : Tangential induced velocities at rotor control points.
- `Cm_wake` : Circulation moment (or moment coefficient) of the wake.
- Additional fields unpacked from `solve_containers`, providing zeroed or reset versions of solver outputs and state estimates.
"""
function return_failed_residual!(
    solver_options::TS,
    TF,
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

    # - Separate out the state variables - #
    vz_rotor, vtheta_rotor, Cm_wake = extract_failed_state_variables(
        solver_options, state_dims, TF
    )

    # - Extract and Reset Cache - #
    # get cache vector of correct types
    solve_container_cache_vec = @views PreallocationTools.get_tmp(
        solve_container_cache, [vz_rotor; vtheta_rotor; Cm_wake]
    )
    solve_containers = withdraw_solve_container_cache(
        solver_options, solve_container_cache_vec, solve_container_cache_dims
    )
    reset_containers!(solve_containers) #note: also zeros out state estimates

    return (; vz_rotor, vtheta_rotor, Cm_wake, solve_containers...)
end

"""
    return_failed_residual!(
        solver_options::CSORSolverOptions,
        TF,
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

Fallback residual constructor for failed CSOR solver execution.

When a solver of type `CSORSolverOptions` fails to converge, this function assembles a well-structured output with all solver state variables zero-initialized. It ensures postprocessing and analysis routines can still run without breaking due to missing data.

# Returns
A named tuple containing:
- `Gamr`: Rotor circulation (zeroed)
- `sigr`: Rotor source strengths (zeroed)
- `gamw`: Wake circulation (zeroed)
- Additional solver container fields from `solve_containers`
"""
function return_failed_residual!(
    solver_options::CSORSolverOptions,
    TF,
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

    # - Separate out the state variables - #
    Gamr, sigr, gamw = extract_failed_state_variables(solver_options, state_dims, TF)

    # - Extract and Reset Cache - #
    # get cache vector of correct types
    solve_container_cache_vec = @views PreallocationTools.get_tmp(
        solve_container_cache, [Gamr; sigr; gamw]
    )
    solve_container_cache_vec .= 0
    solve_containers = withdraw_solve_container_cache(
        solver_options, solve_container_cache_vec, solve_container_cache_dims
    )

    return (; Gamr, sigr, gamw, solve_containers...)
end

function return_failed_residual!(
    solver_options::ModCSORSolverOptions,
    TF,
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

    # - Separate out the state variables - #
    Gamr, sigr, gamw = extract_failed_state_variables(solver_options, state_dims, TF)

    # - Extract and Reset Cache - #
    # get cache vector of correct types
    solve_container_cache_vec = @views PreallocationTools.get_tmp(
        solve_container_cache, [Gamr; sigr; gamw]
    )
    solve_container_cache_vec .= 0
    solve_containers = withdraw_solve_container_cache(
        solver_options, solve_container_cache_vec, solve_container_cache_dims
    )

    return (; Gamr, sigr, gamw, solve_containers...)
end

"""
    return_failed_boundary_layer_outputs(
        boundary_layer_options::AbstractBoundaryLayerOptions,
        Vtan_duct,
        Vtot_duct,
        duct_control_points,
        duct_panel_lengths,
        duct_panel_tangents,
        rotor_tip_radius,
        operating_point,
        reference_parameters,
    )

Generates a zero-filled boundary layer output in the event of solver failure.

This fallback function ensures that the calling interface receives a complete and structurally valid boundary layer solution, even when the boundary layer integration (e.g., HEADS or Green's method) fails to converge or initialize properly.

The dimensionality of state arrays (2 or 3 columns) is determined by the boundary layer model type:
- `HeadsBoundaryLayerOptions`: 2 states (e.g., momentum thickness and shape factor)
- `GreensBoundaryLayerOptions`: 3 states (e.g., with additional energy thickness or auxiliary variables)

# Returns
A named tuple containing zeroed or default values for:
- `stagnation_indices`
- `upper_initial_states`, `upper_solved_states`, `upper_solved_steps`
- `lower_initial_states`, `lower_solved_states`, `lower_solved_steps`
- `surface_length_upper`, `surface_length_lower`
- `stag_point`, `split_ratio`
- `separation_point_ratio_upper`, `separation_point_ratio_lower`
- `cdc_upper`, `cdc_lower`
- `vtdotpv`

# Notes
- Intended for use when boundary layer models fail to initialize or run.
- Allows postprocessing to continue without raising errors.
"""
function return_failed_boundary_layer_outputs(
    boundary_layer_options::HeadsBoundaryLayerOptions,
    Vtan_duct,
    Vtot_duct,
    duct_control_points,
    duct_panel_lengths,
    duct_panel_tangents,
    rotor_tip_radius,
    operating_point,
    reference_parameters,
)
    TF = eltype(Vtan_duct)

    return (;
        stagnation_indices=[1; 1],
        upper_initial_states=zeros(TF, 2),
        upper_solved_states=zeros(TF, 2),
        upper_solved_steps=TF[0],
        lower_initial_states=zeros(TF, 2),
        lower_solved_states=zeros(TF, 2),
        lower_solved_steps=TF[0],
        surface_length_upper=zeros(TF, 2),
        surface_length_lower=zeros(TF, 2),
        stag_point=TF[0.5],
        split_ratio=TF[0.5],
        separation_point_ratio_upper=TF[0.0],
        separation_point_ratio_lower=TF[0.0],
        cdc_upper=TF[0.0],
        cdc_lower=TF[0.0],
        vtdotpv=zeros(TF, 2),
        # boundary_layer_functions_lower,
        # boundary_layer_functions_upper,
    )
end


function return_failed_boundary_layer_outputs(
    boundary_layer_options::GreensBoundaryLayerOptions,
    Vtan_duct,
    Vtot_duct,
    duct_control_points,
    duct_panel_lengths,
    duct_panel_tangents,
    rotor_tip_radius,
    operating_point,
    reference_parameters,
)
    TF = eltype(Vtan_duct)

    return (;
        stagnation_indices=[1; 1],
        upper_initial_states=zeros(TF, 3),
        upper_solved_states=zeros(TF, 3),
        upper_solved_steps=TF[0],
        lower_initial_states=zeros(TF, 3),
        lower_solved_states=zeros(TF, 3),
        lower_solved_steps=TF[0],
        surface_length_upper=zeros(TF, 2),
        surface_length_lower=zeros(TF, 2),
        stag_point=TF[0.5],
        split_ratio=TF[0.5],
        separation_point_ratio_upper=TF[0.0],
        separation_point_ratio_lower=TF[0.0],
        cdc_upper=TF[0.0],
        cdc_lower=TF[0.0],
        vtdotpv=zeros(TF, 2),
        # boundary_layer_functions_lower,
        # boundary_layer_functions_upper,
    )
end


"""
    extract_failed_state_variables(solver_options::TS, dims, TF) where {TS <: AbstractSolverOptions}

Returns zero-initialized arrays for solver state variables in the event of a failure during the solve process.

This utility provides a consistent fallback mechanism by reshaping zero arrays into the expected shapes of the solver's state variables, based on the provided dimension mappings. It supports both external and internal solver configurations.

# Arguments
- `solver_options::TS`: The solver configuration type. Determines whether external (e.g., actuator velocities) or internal (e.g., circulation) states are returned.
- `dims`: A struct containing indexing and reshaping information for all relevant state variables.
- `TF`: The element type for the arrays (typically a floating point type like `Float64`).

# Returns
A named tuple of zero-initialized and reshaped state variables:
- For `ExternalSolverOptions`:
  - `vz_rotor`: Axial velocity at the rotor.
  - `vtheta_rotor`: Tangential velocity at the rotor.
  - `Cm_wake`: Wake momentum coefficient.

- For `InternalSolverOptions`:
  - `Gamr`: Rotor-bound circulation.
  - `sigr`: Rotor-bound source strength.
  - `gamw`: Wake-bound vorticity.

# Notes
- This function is primarily used when a solve attempt fails and a placeholder state must be returned to avoid further errors.
- Shapes and indexing come from the `dims` structure, which must match the solver’s internal layout.
"""
function extract_failed_state_variables(
    solver_options::TS, dims, TF
) where {TS<:ExternalSolverOptions}

    # - Separate out - #
    vz_rotor = @views reshape(zeros(TF, length(dims.vz_rotor.index)), dims.vz_rotor.shape)
    vtheta_rotor = @views reshape(
        zeros(TF, length(dims.vtheta_rotor.index)), dims.vtheta_rotor.shape
    )
    Cm_wake = @views reshape(zeros(TF, length(dims.Cm_wake.index)), dims.Cm_wake.shape)

    return vz_rotor, vtheta_rotor, Cm_wake
end

function extract_failed_state_variables(
    solver_options::TS, dims, TF
) where {TS<:InternalSolverOptions}

    # - Separate out - #
    Gamr = @views reshape(zeros(TF, length(dims.Gamr.index)), dims.Gamr.shape)
    sigr = @views reshape(zeros(TF, length(dims.sigr.index)), dims.sigr.shape)
    gamw = @views reshape(zeros(TF, length(dims.gamw.index)), dims.gamw.shape)

    return Gamr, sigr, gamw
end
