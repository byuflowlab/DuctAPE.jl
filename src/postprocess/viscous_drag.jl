#---------------------------------#
#          Viscous Drag           #
#---------------------------------#

"""
    squire_young(d2, Ue, Uinf, H12)

Squire-Young formula for the viscous drag coeffiecient of one side of a 2D body.

# Arguments:
- `d2::Float` : Momentum thickness at separation extrapolated back to the trailing edge (d2 = d2sep+rsep-rTE)
- `Ue::Float` : Edge velocity at separation point
- `Uinf::Float` : Freestream velocity
- `H12::Float` : Boundary layer shape factor at separation point

Note: if no separation occurs, the inputs are simply the final values for the boundary layer.

# Returns:
- `cdc::Float` : viscous drag coefficient times reference chord
"""
function squire_young(d2, Ue, Uinf, H12)
    # note: formula also divides by chord, but we're going to multiply by chord later.
    return 2.0 * d2 * FLOWMath.abs_smooth(Ue / Uinf, 2.0 * eps())^((5.0 + H12) / 2.0)
end

"""
    total_viscous_drag_duct(cd_upper, cd_lower, exit_radius, Vref, rhoinf)

Calculate the total viscous drag of the duct from Squire-Young drag coefficients, integrated about exit circumference.

# Arguments:
- `cdc_upper::Float` : upper side drag coefficient times refernce chord
- `cdc_lower::Float` : lower side drag coefficient times refernce chord
- `exit_radius::Float` : radius used for integrating circumferentially
- `Vref::Float` : reference velocity (Vinf)
- `rhoinf::Float` : freestream density

# Returns:
- `viscous_drag::Float` : viscous drag on duct
"""
function total_viscous_drag_duct(cdc_upper, cdc_lower, exit_radius, Vref, rhoinf)
    # note: cdc's are already times chord, so no need to have a separate chord variable

    # drag per unit length
    dprime = 0.5 * rhoinf * Vref^2 * (cdc_upper + cdc_lower)

    # drag of annular airfoil
    return dprime * 2.0 * pi * exit_radius
end

#---------------------------------#
#          HEAD'S METHOD          #
#---------------------------------#

"""
    compute_single_side_drag_coefficient_head(
        steps,
        exit_radius,
        operating_point,
        boundary_layer_options;
        verbose=false,
    )

Solve integral boundary layer and obtain viscous drag coefficient from Squire-Young formula for one side of the duct (side being defined as portion of surface on once side of the stagnation point)

# Arguments:
- `steps::Vector{Float}` : positions along surface for integration
- `exit_radius::Float` : radius at duct trailing edge (casing side)
- `operating_point::Float` : OperatingPoint object
- `boundary_layer_functions::NamedTuple` : Various Akima splines and other functions for boundary layer values
- `boundary_layer_options::BoundaryLayerOptions` : BoundaryLayerOptions object

# Returns:
- `cd::Float` : viscous drag coefficient
"""
function compute_single_side_drag_coefficient_head(
    steps,
    exit_radius,
    operating_point,
    boundary_layer_functions,
    boundary_layer_options;
    verbose=false,
)
    (; edge_density, edge_velocity, edge_viscosity) = boundary_layer_functions

    # verbose && printdebug("initial velocity: ", edge_velocity(steps[1]))
    # verbose && printdebug("initial density: ", edge_density(steps[1]))
    # verbose && printdebug("initial viscosity: ", edge_viscosity(steps[1]))

    # - Initialize Boundary Layer States - #
    H10 = 10.6
    d20 =
        0.036 * steps[1] / calculate_Re(
            edge_density(steps[1]),
            edge_velocity(steps[1]),
            steps[1],
            edge_viscosity(steps[1]),
        )
    initial_states = [d20; H10]
    H0 = 1.28

    usep, Hsep, s_sep, sepid, usol, stepsol = solve_head_boundary_layer!(
        boundary_layer_residual_head,
        boundary_layer_options.rk,
        [initial_states, H0],
        steps,
        (; boundary_layer_functions..., boundary_layer_options.separation_criteria);
        verbose=false,
    )

    # verbose && println(s_sep / steps[end])

    cdsq = squire_young(
        usep[1], boundary_layer_functions.edge_velocity(s_sep), operating_point.Vinf[], Hsep
    )

    cd = FLOWMath.ksmin([boundary_layer_options.separation_penalty; cdsq])

    cdadd = FLOWMath.ksmax(
        [
            0.0
            FLOWMath.linear(
                [0.0; steps[end - boundary_layer_options.separation_allowance]],
                [boundary_layer_options.separation_penalty; 0.0],
                s_sep,
            )
        ],
        100,
    )

    cd += cdadd

    return cd, usol, stepsol, s_sep / steps[end]
end

"""
    compute_viscous_drag_duct(
        boundary_layer_options::BoundaryLayerOptions,
        Vtan_duct,
        cp_duct,
        duct_panel_lengths,
        exit_radius,
        operating_point,
    )

Determine total, dimensional viscous drag on the duct.

# Arguments:
- `boundary_layer_options::BoundaryLayerOptions` : BoundaryLayerOptions object, used for dispatch as well
- `Vtan_duct::Vector{Float}` : tangential velocity magnitudes for the entire duct
- `duct_panel_lengths::Vector{Float}` : panel lengths for the entire duct
- `exit_radius::Float` : radius at duct trailing edge (casing side)
- `operating_point::Float` : OperatingPoint object

# Returns:
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
"""
function compute_viscous_drag_duct(
    boundary_layer_options::HeadsBoundaryLayerOptions,
    Vtan_duct,
    Vtot_duct,
    duct_panel_lengths,
    duct_panel_tangents,
    exit_radius,
    operating_point;
    verbose=false,
)

    # find stagnation point
    s_upper, s_lower, stag_ids, split_ratio = split_at_stagnation_point(
        duct_panel_lengths, duct_panel_tangents, Vtot_duct
    )

    # set up boundary layer solve parameters
    if !isnothing(s_upper)
        boundary_layer_functions_upper = setup_boundary_layer_functions_head(
            s_upper,
            [0.0; Vtan_duct[stag_ids[2]:end]],
            operating_point,
            boundary_layer_options;
            verbose=verbose,
        )
    end

    # set up boundary layer solve parameters
    boundary_layer_functions_lower = setup_boundary_layer_functions_head(
        s_lower,
        [0.0; Vtan_duct[stag_ids[1]:-1:1]],
        operating_point,
        boundary_layer_options;
        verbose=verbose,
    )

    # - Set integration steps - #
    if !isnothing(s_upper)
        # upper side
        upper_steps =
            set_boundary_layer_steps(
                boundary_layer_options.n_steps,
                boundary_layer_options.first_step_size,
                s_upper[end] - boundary_layer_options.offset,
            ) .+ boundary_layer_options.offset
    end

    # lower side
    lower_steps =
        set_boundary_layer_steps(
            boundary_layer_options.n_steps,
            boundary_layer_options.first_step_size,
            s_lower[end] - boundary_layer_options.offset,
        ) .+ boundary_layer_options.offset

    # - Get drag coeffients - #

    if !isnothing(s_upper)
        cdc_upper, usol_upper, stepsol_upper, s_sep_upper = compute_single_side_drag_coefficient_head(
            upper_steps,
            exit_radius,
            operating_point,
            boundary_layer_functions_upper,
            (;
                boundary_layer_options.rk,
                boundary_layer_options.separation_criteria,
                separation_allowance=boundary_layer_options.separation_allowance_upper,
                separation_penalty=boundary_layer_options.separation_penalty_upper,
            );
            verbose=verbose,
        )
    else
        cdc_upper = 0.0
        usol_upper = -ones(eltype(Vtan_duct), 3)
        stepsol_upper = -1
        s_sep_upper = -1.0
    end

    cdc_lower, usol_lower, stepsol_lower, s_sep_lower = compute_single_side_drag_coefficient_head(
        lower_steps,
        exit_radius,
        operating_point,
        boundary_layer_functions_lower,
        (;
            boundary_layer_options.rk,
            boundary_layer_options.separation_criteria,
            separation_allowance=boundary_layer_options.separation_allowance_lower,
            separation_penalty=boundary_layer_options.separation_penalty_lower,
        );
        verbose=verbose,
    )

    # Calculate total viscous drag
    total_drag = total_viscous_drag_duct(
        cdc_upper, cdc_lower, exit_radius, operating_point.Vinf[], operating_point.rhoinf[]
    )

    # Return
    return total_drag,
    (;
        stagnation_indices=stag_ids,
        upper_solved_states=usol_upper,
        upper_solved_steps=stepsol_upper,
        lower_solved_states=usol_lower,
        lower_solved_steps=stepsol_lower,
        surface_length_upper=s_upper,
        surface_length_lower=s_lower,
        split_ratio,
        separation_point_ratio_upper=s_sep_upper,
        separation_point_ratio_lower=s_sep_lower,
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

# Arguments:
- `vtan_duct_TE::Float` : tangential velocity at the duct trailing edge
- `duct_chord::Float` : length of duct
- `TE_radius::Float` : radius of the trailing edge point
- `operating_point::OperatingPoint` : OperatingPoint object
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
