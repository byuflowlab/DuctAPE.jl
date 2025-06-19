#---------------------------------#
#          Viscous Drag           #
#---------------------------------#

"""
    squire_young(d2, Ue, Uinf, H12)

Squire-Young formula for the viscous drag coeffiecient of one side of a 2D body.

# Arguments
- `d2::Float` : Momentum thickness at separation extrapolated back to the trailing edge (d2 = d2sep+rsep-rTE)
- `Ue::Float` : Edge velocity at separation point
- `Uinf::Float` : Freestream velocity
- `H12::Float` : Boundary layer shape factor at separation point

Note: if no separation occurs, the inputs are simply the final values for the boundary layer.

# Returns
- `cdc::Float` : viscous drag coefficient times reference chord
"""
function squire_young(d2, Ue, Uinf, H12)
    # note: formula also divides by chord, but we're going to multiply by chord later.
    return 2.0 * d2 * FLOWMath.abs_smooth(Ue / Uinf, 2.0 * eps())^((5.0 + H12) / 2.0)
end

"""
    total_viscous_drag_duct(cd_upper, cd_lower, rotor_tip_radius, Vref, rhoinf)

Calculate the total viscous drag of the duct from Squire-Young drag coefficients, integrated about exit circumference.

# Arguments
- `cdc_upper::Float` : upper side drag coefficient times refernce chord
- `cdc_lower::Float` : lower side drag coefficient times refernce chord
- `rotor_tip_radius::Float` : radius used for integrating circumferentially
- `Vref::Float` : reference velocity (Vinf)
- `rhoinf::Float` : freestream density

# Returns
- `viscous_drag::Float` : viscous drag on duct
"""
function total_viscous_drag_duct(cdc_upper, cdc_lower, rotor_tip_radius, Vref, rhoinf)
    # note: cdc's are already times chord, so no need to have a separate chord variable

    # drag per unit length
    dprime = 0.5 * rhoinf * Vref^2 * (cdc_upper + cdc_lower)

    # drag of annular airfoil
    return dprime * 2.0 * pi * rotor_tip_radius
end

"""
    compute_single_side_drag_coefficient(
        boundary_layer_options::HeadsBoundaryLayerOptions,
        steps,
        rotor_tip_radius,
        Vref,
        boundary_layer_functions,
        separation_options;
        verbose=false,
    )

Solve the integral boundary layer equations on a single side of the duct (as defined relative to the stagnation point) 
and compute the viscous drag coefficient using the Squire-Young formula. This method accounts for potential separation 
and includes optional penalties when separation occurs.

# Arguments
- `boundary_layer_options::HeadsBoundaryLayerOptions` : Configuration options for the boundary layer solver.
- `steps::Vector{Float}` : Streamwise positions along the surface for integration.
- `rotor_tip_radius::Float` : Tip radius of the rotor or casing at the trailing edge (used for scaling).
- `Vref::Vector{Float}` : One-element vector containing reference freestream velocity.
- `boundary_layer_functions::NamedTuple` : Collection of splines and functions for evaluating boundary layer properties (e.g., edge velocity, curvature).
- `separation_options::NamedTuple` : Additional configuration for separation handling (e.g., allowances and penalties).

# Keyword Arguments
- `verbose::Bool=false` : If true, print debug information during computation.

# Returns
- `cd::Float` : Viscous drag coefficient from the Squire-Young formula, with any separation penalties applied.
- `u_init` : Initial boundary layer state.
- `usol` : Solution history of boundary layer states.
- `stepsol` : Corresponding streamwise positions for the solution.
- `s_sep` : Streamwise position of separation (if detected).
- `Hsep` : Shape factor at separation (if detected).

# Notes
- The function internally calls the appropriate solver (e.g., head or green boundary layer solver) and applies separation penalties if configured.
- `Vref` is expected to be a one-element vector for consistency across the codebase.
"""
function compute_single_side_drag_coefficient(
    boundary_layer_options::HeadsBoundaryLayerOptions,
    steps,
    rotor_tip_radius,
    Vref,
    boundary_layer_functions,
    separation_options;
    verbose=false,
)
    return compute_single_side_drag_coefficient(
        initialize_head_states,
        solve_head_boundary_layer!,
        steps,
        rotor_tip_radius,
        Vref,
        boundary_layer_functions,
        (;
            boundary_layer_options.solver_type,
            boundary_layer_options.ode,
            boundary_layer_options.separation_criteria,
            boundary_layer_options.first_step_size,
            boundary_layer_options.dy_eps,
            boundary_layer_options.H1_eps,
            boundary_layer_options.H_eps,
            boundary_layer_options.terminate,
            boundary_layer_options.return_last_max_shape_factor,
            boundary_layer_options.cutoff_Hsep,
            separation_options...,
        );
        verbose=verbose,
    )
end

# function compute_single_side_drag_coefficient(
#     boundary_layer_options::GreensBoundaryLayerOptions,
#     steps,
#     rotor_tip_radius,
#     Vref,
#     boundary_layer_functions,
#     separation_options;
#     verbose=false,
# )
#     return compute_single_side_drag_coefficient(
#         initialize_green_states,
#         solve_green_boundary_layer!,
#         steps,
#         rotor_tip_radius,
#         Vref,
#         boundary_layer_functions,
#         (;
#             boundary_layer_options...,
#             boundary_layer_options.ode,
#             boundary_layer_options.separation_criteria,
#             boundary_layer_options.first_step_size,
#             boundary_layer_options.lambda,
#             boundary_layer_options.longitudinal_curvature,
#             boundary_layer_options.lateral_strain,
#             boundary_layer_options.dilation,
#             boundary_layer_options.terminate,
#             boundary_layer_options.return_last_max_shape_factor,
#             separation_options...,
#         );
#         verbose=verbose,
#     )
# end

"""
    compute_single_side_drag_coefficient(
        initialize_states,
        solve_boundary_layer,
        steps,
        rotor_tip_radius,
        Vref,
        boundary_layer_functions,
        single_side_boundary_layer_options;
        verbose=false,
    )

Generic boundary layer drag coefficient calculator for one side of a duct. Uses provided initialization 
and solver functions to compute the boundary layer solution and applies the Squire-Young formula 
(with optional separation penalties) to compute the viscous drag coefficient.

# Arguments
- `initialize_states::Function` : Function to initialize the boundary layer state at the starting location.
- `solve_boundary_layer::Function` : Solver function that integrates the boundary layer equations along `steps`.
- `steps::Vector{Float}` : Streamwise positions along the surface.
- `rotor_tip_radius::Float` : Tip radius at the trailing edge, used for scaling.
- `Vref::Vector{Float}` : One-element vector containing the freestream velocity magnitude.
- `boundary_layer_functions::NamedTuple` : Contains splines/functions for edge velocity, curvature, etc.
- `single_side_boundary_layer_options::NamedTuple` : Options for the solver, including ODE type, separation criteria, penalties, etc.

# Keyword Arguments
- `verbose::Bool=false` : If true, enables debug printing during boundary layer computation.

# Returns
- `cd::Float` : Viscous drag coefficient from the Squire-Young formula, with separation penalties if applicable.
- `u_init` : Initial boundary layer state vector.
- `usol` : Solution history (boundary layer state at each step).
- `stepsol` : Streamwise positions corresponding to `usol`.
- `s_sep` : Streamwise location of separation (if detected).
- `Hsep` : Shape factor at separation (if detected).

# Notes
- The drag coefficient is augmented by a penalty if separation occurs and penalties are configured.
- `Vref` is assumed to be a one-element vector for consistency with other API components.
- This form is intended for use inside higher-level solver-specific wrappers (e.g., head or green boundary layer solvers).
"""
function compute_single_side_drag_coefficient(
    initialize_states,
    solve_boundary_layer,
    steps,
    rotor_tip_radius,
    Vref,
    boundary_layer_functions,
    single_side_boundary_layer_options;
    verbose=false,
)
    u_init = initialize_states(boundary_layer_functions, steps[1]; verbose=verbose)

    usep, Hsep, s_sep, usol, stepsol = solve_boundary_layer(
        single_side_boundary_layer_options.solver_type,
        single_side_boundary_layer_options.ode,
        u_init,
        steps,
        (; boundary_layer_functions..., single_side_boundary_layer_options...);
        verbose=verbose,
    )

    # printdebug("Hsep:", Hsep)
    # printdebug("δ_2:", usep[1])

    cdsqy = squire_young(
        usep[1], boundary_layer_functions.edge_velocity(s_sep), Vref[], Hsep
    )

    if single_side_boundary_layer_options.separation_penalty < eps()
        cd = cdsqy
    else
        cdadd = separation_penalty(
            s_sep,
            steps,
            single_side_boundary_layer_options.separation_allowance,
            single_side_boundary_layer_options.separation_penalty,
        )

        cdsqy += cdadd
    end

    return cdsqy, u_init, usol, stepsol, s_sep, Hsep
end

"""
    compute_viscous_drag_duct(
        boundary_layer_options::BoundaryLayerOptions,
        Vtan_duct,
        cp_duct,
        duct_panel_lengths,
        rotor_tip_radius,
        operating_point,
        reference_parameters,
    )

Determine total, dimensional viscous drag on the duct.

# Arguments
- `boundary_layer_options::BoundaryLayerOptions` : BoundaryLayerOptions object, used for dispatch as well
- `Vtan_duct::Vector{Float}` : tangential velocity magnitudes for the entire duct
- `duct_panel_lengths::Vector{Float}` : panel lengths for the entire duct
- `rotor_tip_radius::Float` : radius at duct trailing edge (casing side)
- `operating_point::OperatingPoint` : OperatingPoint object
- `reference_parameters::ReferenceParameters` : ReferenceParameters object

# Returns
- `duct_viscous_drag::Float` : total viscous drag of duct
- `boundary_layer_outputs::NamedTuple` : named tuple of various boundary layer related outputs:
  - `stagnation_indices`
  - `upper_solved_states`
  - `upper_solved_steps`
  - `lower_solved_states`
  - `lower_solved_steps`
  - `surface_length_upper`
  - `surface_length_lower`
  - `split_ratio`
  - `separation_point_ratio_upper`
  - `separation_point_ratio_lower`
  - `cdc_upper`
  - `cdc_lower`
"""
function compute_viscous_drag_duct(
    boundary_layer_options::HeadsBoundaryLayerOptions,
    Vtan_duct,
    Vtot_duct,
    duct_control_points,
    duct_panel_lengths,
    duct_panel_tangents,
    rotor_tip_radius,
    operating_point,
    reference_parameters;
    verbose=false,
)
    return compute_viscous_drag_duct(
        setup_boundary_layer_functions_head,
        boundary_layer_options,
        Vtan_duct,
        Vtot_duct,
        duct_control_points,
        duct_panel_lengths,
        duct_panel_tangents,
        rotor_tip_radius,
        operating_point,
        reference_parameters;
        verbose=verbose,
    )
end

# function compute_viscous_drag_duct(
#     boundary_layer_options::GreensBoundaryLayerOptions,
#     Vtan_duct,
#     Vtot_duct,
#     duct_control_points,
#     duct_panel_lengths,
#     duct_panel_tangents,
#     rotor_tip_radius,
#     operating_point,
#     reference_parameters;
#     verbose=false,
# )
#     return compute_viscous_drag_duct(
#         setup_boundary_layer_functions_green,
#         boundary_layer_options,
#         Vtan_duct,
#         Vtot_duct,
#         duct_control_points,
#         duct_panel_lengths,
#         duct_panel_tangents,
#         rotor_tip_radius,
#         operating_point,
#         reference_parameters;
#         verbose=verbose,
#     )
# end

"""
    compute_viscous_drag_duct(
        setup_boundary_layer_functions,
        boundary_layer_options,
        Vtan_duct,
        Vtot_duct,
        duct_control_points,
        duct_panel_lengths,
        duct_panel_tangents,
        rotor_tip_radius,
        operating_point,
        reference_parameters;
        verbose=false,
    )

Compute the total viscous drag for a duct by solving the boundary layer on both sides of the stagnation point 
and integrating the drag contribution using Squire-Young-based formulas with optional separation penalties.

# Arguments
- `setup_boundary_layer_functions::Function` : Function that builds splines and helper functions for the boundary layer solver.
- `boundary_layer_options::BoundaryLayerOptions` : Options controlling solver type, step sizes, separation criteria, penalties, etc.
- `Vtan_duct::Vector{Float}` : Tangential velocity at duct control points.
- `Vtot_duct::Vector{Float}` : Total velocity at duct control points.
- `duct_control_points::Matrix{Float}` : Coordinates of control points on the duct.
- `duct_panel_lengths::Vector{Float}` : Lengths of duct surface panels.
- `duct_panel_tangents::Matrix{Float}` : Tangent vectors for duct surface panels.
- `rotor_tip_radius::Float` : Tip radius of the rotor (used for scaling viscous drag).
- `operating_point::NamedTuple` : Contains operating condition data, including freestream density `rhoinf`.
- `reference_parameters::NamedTuple` : Contains reference freestream velocity `Vref`.

# Keyword Arguments
- `verbose::Bool=false` : If true, enables diagnostic printing during boundary layer computations.

# Returns
- `total_drag::Float` : Total viscous drag (dimensional).
- `details::NamedTuple` : Contains detailed solution data, including:
  - `stagnation_indices` : Indices of stagnation points (upper, lower).
  - `stag_point` : Coordinates of stagnation point.
  - `split_ratio` : Surface distance split ratio at stagnation.
  - `surface_length_upper` / `surface_length_lower` : Arc lengths on upper and lower sides.
  - `upper_steps` / `lower_steps` : Integration steps used.
  - `cdc_upper` / `cdc_lower` : Drag coefficients for upper and lower sides.
  - `upper_initial_states` / `lower_initial_states` : Initial boundary layer states.
  - `upper_solved_states` / `lower_solved_states` : Solution state history.
  - `upper_solved_steps` / `lower_solved_steps` : Solution step positions.
  - `Hsep_upper` / `Hsep_lower` : Shape factor at separation.
  - `s_sep_upper` / `s_sep_lower` : Streamwise position of separation.
  - `separation_point_ratio_upper` / `separation_point_ratio_lower` : Separation location as fraction of total length.
  - `vtdotpv` : Velocity tangent dot products at panel midpoints.
  - `boundary_layer_functions_upper` / `boundary_layer_functions_lower` : Functions used for boundary layer evaluation.

# Notes
- The function automatically detects the stagnation point and splits the duct surface for integration.
- Handles cases where no upper side exists (e.g., open ducts or numerical artifacts).
- The final drag is dimensionalized using `rotor_tip_radius`, `Vref`, and `rhoinf`.
"""
function compute_viscous_drag_duct(
    setup_boundary_layer_functions,
    boundary_layer_options,
    Vtan_duct,
    Vtot_duct,
    duct_control_points,
    duct_panel_lengths,
    duct_panel_tangents,
    rotor_tip_radius,
    operating_point,
    reference_parameters;
    verbose=false,
)
    # find stagnation point
    s_upper, s_lower, stag_ids, stag_point, split_ratio, dots = split_at_stagnation_point(
        duct_panel_lengths,
        duct_panel_tangents,
        Vtot_duct,
        Vtan_duct,
        boundary_layer_options.first_step_size,
    )

    # set up boundary layer solve parameters
    if !isnothing(s_upper)
        boundary_layer_functions_upper = setup_boundary_layer_functions(
            s_upper,
            [0.0; Vtan_duct[stag_ids[2]:end]],
            duct_control_points,
            operating_point,
            boundary_layer_options;
            verbose=verbose,
        )
    end

    # set up boundary layer solve parameters
    boundary_layer_functions_lower = setup_boundary_layer_functions(
        s_lower,
        [0.0; Vtan_duct[stag_ids[1]:-1:1]],
        duct_control_points,
        operating_point,
        boundary_layer_options;
        verbose=verbose,
    )

    # - Set integration steps - #
    if !isnothing(s_upper)
        # upper side
        if boundary_layer_options.solver_type == RK()
            if isnothing(boundary_layer_options.upper_step_size)
                upper_steps =
                    set_boundary_layer_steps(
                        boundary_layer_options.n_steps,
                        boundary_layer_options.first_step_size,
                        s_upper[end] - boundary_layer_options.offset,
                    ) .+ boundary_layer_options.offset
            else
                upper_steps =
                    range(
                        0.0,
                        s_upper[end] - boundary_layer_options.offset;
                        step=boundary_layer_options.upper_step_size,
                    ) .+ boundary_layer_options.offset
            end
        else
            upper_steps =
                range(
                    0.0,
                    s_upper[end] - boundary_layer_options.offset;
                    length=length(s_upper),
                ) .+ boundary_layer_options.offset
        end
    end

    # lower side
    if boundary_layer_options.solver_type == RK()
        if isnothing(boundary_layer_options.lower_step_size)
            lower_steps =
                set_boundary_layer_steps(
                    boundary_layer_options.n_steps,
                    boundary_layer_options.first_step_size,
                    s_lower[end] - boundary_layer_options.offset,
                ) .+ boundary_layer_options.offset
        else
            lower_steps =
                range(
                    0.0,
                    s_lower[end] - boundary_layer_options.offset;
                    step=boundary_layer_options.lower_step_size,
                ) .+ boundary_layer_options.offset
        end
    else
        lower_steps =
            range(
                0.0, s_lower[end] - boundary_layer_options.offset; length=length(s_lower)
            ) .+ boundary_layer_options.offset
    end

    # - Get drag coeffients - #

    if !isnothing(s_upper)
        cdc_upper, u_init_upper, usol_upper, stepsol_upper, s_sep_upper, Hsep_upper = compute_single_side_drag_coefficient(
            boundary_layer_options,
            upper_steps,
            rotor_tip_radius,
            reference_parameters.Vref,
            boundary_layer_functions_upper,
            (;
                separation_allowance=boundary_layer_options.separation_allowance_upper,
                separation_penalty=boundary_layer_options.separation_penalty_upper,
            );
            verbose=verbose,
        )
    else
        cdc_upper = 0.0
        u_init_upper = 0.0
        boundary_layer_functions_upper = (;)
        usol_upper = -ones(eltype(Vtan_duct), 3)
        stepsol_upper = -1
        s_sep_upper = -1.0
        Hsep_upper = 0.0
    end

    cdc_lower, u_init_lower, usol_lower, stepsol_lower, s_sep_lower, Hsep_lower = compute_single_side_drag_coefficient(
        boundary_layer_options,
        lower_steps,
        rotor_tip_radius,
        reference_parameters.Vref,
        boundary_layer_functions_lower,
        (;
            separation_allowance=boundary_layer_options.separation_allowance_lower,
            separation_penalty=boundary_layer_options.separation_penalty_lower,
        );
        verbose=verbose,
    )

    # Calculate total viscous drag
    total_drag = total_viscous_drag_duct(
        cdc_upper,
        cdc_lower,
        rotor_tip_radius,
        reference_parameters.Vref[],
        operating_point.rhoinf[],
    )

    # Return
    return total_drag,
    (;
        stagnation_indices=stag_ids,
        upper_initial_states=u_init_upper,
        Hsep_upper,
        upper_solved_states=usol_upper,
        upper_solved_steps=stepsol_upper,
        Hsep_lower,
        lower_initial_states=u_init_lower,
        lower_solved_states=usol_lower,
        lower_solved_steps=stepsol_lower,
        surface_length_upper=s_upper,
        surface_length_lower=s_lower,
        stag_point,
        split_ratio,
        s_sep_upper,
        s_sep_lower,
        upper_steps = !isnothing(s_upper) ? upper_steps : [],
        lower_steps,
        separation_point_ratio_upper= !isnothing(s_upper) ? s_sep_upper/upper_steps[end] : 1.0,
        separation_point_ratio_lower=s_sep_lower/lower_steps[end],
        cdc_upper=cdc_upper,
        cdc_lower=cdc_lower,
        vtdotpv=dots,
        boundary_layer_functions_lower,
        boundary_layer_functions_upper,
    )
end

#---------------------------------#
#    SCHLICHTING APPROXIMATION    #
#---------------------------------#

"""
    compute_viscous_drag_duct_schlichting(
        vtan_duct_TE, duct_chord, TE_radius, operating_point
    )

Computes Schlichting approximation of skin friction drag dimensionalized to drag per unit length using duct chord, and using the trailing edge circuference as the total length.

# Arguments
- `vtan_duct_TE::Float` : tangential velocity at the duct trailing edge
- `duct_chord::Float` : length of duct
- `TE_radius::Float` : radius of the trailing edge point
- `operating_point::OperatingPoint` : OperatingPoint object

# Returns
- `drag::Float` : Estimated viscous drag force (dimensional, per unit length).
- `Cf::Vector{Float}` : Skin friction coefficients used in the calculation.
"""
function compute_viscous_drag_duct_schlichting(
    vtan_duct_TE, duct_chord, TE_radius, operating_point
)
    (; rhoinf, muinf) = operating_point
    q = 0.5 * rhoinf[] * vtan_duct_TE .^ 2

    Cf =
        0.074 ./
        calculate_Re.(Ref(rhoinf[]), vtan_duct_TE, Ref(duct_chord), Ref(muinf[])) .^ 0.2

    return (Cf[1] * q[1] + Cf[2] * q[2]) * duct_chord * 2.0 * pi * TE_radius, Cf
end
