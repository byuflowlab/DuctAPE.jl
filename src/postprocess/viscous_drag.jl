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
    return 2.0 * d2 * (Ue / Uinf)^((5.0 + H12) / 2.0)
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

    verbose && printdebug("initial velocity: ", edge_velocity(steps[1]))
    verbose && printdebug("initial density: ", edge_density(steps[1]))
    verbose && printdebug("initial viscosity: ", edge_viscosity(steps[1]))

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

    usep, Hsep, s_sep, sepid = solve_head_boundary_layer!(
        boundary_layer_residual_head,
        boundary_layer_options.rk,
        [initial_states, H0],
        steps,
        (; boundary_layer_functions..., boundary_layer_options.separation_criteria);
        verbose=verbose,
    )

    println(s_sep)

    # println("s_sep - steps[end]: ", abs(s_sep - steps[end]))

    cdsq = squire_young(
        usep[1], boundary_layer_functions.edge_velocity(s_sep), operating_point.Vinf[], Hsep
    )

    # println("cdsq: ", cdsq)

    cd = FLOWMath.ksmin([boundary_layer_options.separation_penalty; cdsq])

    # println("cd nosep: ", cd)

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

    # println("cdadd: ", cdadd)

    # cd += cdadd

    # println("cd w/sep: ", cd)
    # println()

    return cd
end

"""
    compute_viscous_drag_duct(
        boundary_layer_options::BoundaryLayerOptions,
        vtan_duct,
        cp_duct,
        duct_control_points,
        duct_panel_lengths,
        exit_radius,
        operating_point,
    )

Determine total, dimensional viscous drag on the duct.

# Arguments:
- `boundary_layer_options::BoundaryLayerOptions` : BoundaryLayerOptions object, used for dispatch as well
- `vtan_duct::Vector{Float}` : tangential velocity magnitudes for the entire duct
- `duct_control_points::Matrix{Float}` : control point positions for the entire duct
- `duct_panel_lengths::Vector{Float}` : panel lengths for the entire duct
- `exit_radius::Float` : radius at duct trailing edge (casing side)
- `operating_point::Float` : OperatingPoint object

# Returns:
- `duct_viscous_drag::Float` : total viscous drag of duct
"""
function compute_viscous_drag_duct(
    boundary_layer_options::HeadsBoundaryLayerOptions,
    vtan_duct,
    sid,
    duct_control_points,
    duct_panel_lengths,
    exit_radius,
    operating_point;
    verbose=false,
)

    # find stagnation point
    s_upper, s_lower = split_at_stagnation_point(duct_panel_lengths, sid)

    # set up boundary layer solve parameters
    boundary_layer_functions_upper = setup_boundary_layer_functions_head(
        s_upper,
        vtan_duct[sid:end],
        duct_control_points[:, sid:end],
        operating_point,
        boundary_layer_options;
        verbose=verbose,
    )

    # set up boundary layer solve parameters
    boundary_layer_functions_lower = setup_boundary_layer_functions_head(
        s_lower,
        vtan_duct[sid:-1:1],
        duct_control_points[:, sid:-1:1],
        operating_point,
        boundary_layer_options;
        verbose=verbose,
    )

    # - Set integration steps - #
    # upper side
    upper_steps =
        set_boundary_layer_steps(
            boundary_layer_options.n_steps,
            boundary_layer_options.first_step_size,
            s_upper[end] - boundary_layer_options.offset,
        ) .+ boundary_layer_options.offset

    # lower side
    lower_steps =
        set_boundary_layer_steps(
            boundary_layer_options.n_steps,
            boundary_layer_options.first_step_size,
            s_lower[end] - boundary_layer_options.offset,
        ) .+ boundary_layer_options.offset

    # - Get drag coeffients - #

    # upper side
    print("upper: ")

    cdc_upper = compute_single_side_drag_coefficient_head(
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

    # lower side
    print("lower: ")

    cdc_lower = compute_single_side_drag_coefficient_head(
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

    # Return total viscous drag
    return total_viscous_drag_duct(
        cdc_upper, cdc_lower, exit_radius, operating_point.Vinf[], operating_point.rhoinf[]
    )
end

#---------------------------------#
#  GREEN's LAG-ENTRAINMENT METHOD #
#---------------------------------#

"""
    compute_single_side_drag_coefficient_green(
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
function compute_single_side_drag_coefficient_green(
    steps,
    exit_radius,
    operating_point,
    boundary_layer_functions,
    boundary_layer_options;
    verbose=false,
)

    # - Initialize Boundary Layer States - #
    initial_states, Cf0, H12_0 = initialize_turbulent_boundary_layer_states_green(
        steps[1],
        boundary_layer_functions.r_coords(steps[1]),
        boundary_layer_functions.edge_velocity(steps[1]),
        boundary_layer_functions.edge_mach(steps[1]),
        boundary_layer_functions.edge_density(steps[1]),
        boundary_layer_functions.edge_viscosity(steps[1]);
        verbose=verbose,
    )

    usep, H12sep, s_sep, sepid = solve_green_boundary_layer!(
        boundary_layer_residual_green,
        boundary_layer_options.rk,
        [initial_states, Cf0, H12_0],
        steps,
        (;
            boundary_layer_functions...,
            boundary_layer_options.lambda,
            boundary_layer_options.longitudinal_curvature,
            boundary_layer_options.lateral_strain,
            boundary_layer_options.dilation,
        );
        verbose=verbose,
    )

    if sepid[] < length(steps)
        println(
            "separated at sepid=$(sepid[]) of $(length(steps)): $(steps[sepid[]]) of $(maximum(steps)-minimum(steps))",
        )
    end

    println("s_sep - steps[end]: ", abs(s_sep - steps[end]))
    cdsq = squire_young(
        usep[1],
        boundary_layer_functions.edge_velocity(s_sep),
        operating_point.Vinf[],
        H12sep,
    )

    println("cdsq: ", cdsq)

    cd = FLOWMath.ksmin([
        2.0
        cdsq
    ])

    println("cd nosep: ", cd)

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

    println("cdadd: ", cdadd)

    cd += cdadd

    println("cd w/sep: ", cd)
    println()

    return cd
end

# has a docstring for the other dispatch above
function compute_viscous_drag_duct(
    boundary_layer_options::GreensBoundaryLayerOptions,
    vtan_duct,
    # cp_duct,
    sid,
    duct_control_points,
    duct_panel_lengths,
    exit_radius,
    operating_point;
    verbose=false,
)

    # # find stagnation point
    # s, s_stagnation, lower_length, upper_length = split_at_stagnation_point(
    #     duct_panel_lengths, cp_duct
    # )

    # find stagnation point
    s_upper, s_lower = split_at_stagnation_point(duct_panel_lengths, sid)

    # set up boundary layer solve parameters
    boundary_layer_functions_upper = setup_boundary_layer_functions_green(
        s_upper,
        vtan_duct[sid:end],
        duct_control_points[:, sid:end],
        operating_point,
        boundary_layer_options;
        verbose=verbose,
    )

    # set up boundary layer solve parameters
    boundary_layer_functions_lower = setup_boundary_layer_functions_green(
        s_lower,
        vtan_duct[sid:-1:1],
        duct_control_points[:, sid:-1:1],
        operating_point,
        boundary_layer_options;
        verbose=verbose,
    )

    # - Set integration steps - #
    # upper side
    upper_steps =
        set_boundary_layer_steps(
            boundary_layer_options.n_steps,
            boundary_layer_options.first_step_size,
            s_upper[end] - boundary_layer_options.offset,
        ) .+ boundary_layer_options.offset

    # lower side
    lower_steps =
        set_boundary_layer_steps(
            boundary_layer_options.n_steps,
            boundary_layer_options.first_step_size,
            s_lower[end] - boundary_layer_options.offset,
        ) .+ boundary_layer_options.offset

    # - Get drag coeffients - #

    # upper side
    cdc_upper = compute_single_side_drag_coefficient_green(
        upper_steps,
        exit_radius,
        operating_point,
        boundary_layer_functions_upper,
        (;
            boundary_layer_options.lambda,
            boundary_layer_options.longitudinal_curvature,
            boundary_layer_options.lateral_strain,
            boundary_layer_options.dilation,
            boundary_layer_options.rk,
            separation_allowance=boundary_layer_options.separation_allowance_upper,
            separation_penalty=boundary_layer_options.separation_penalty_upper,
        );
        verbose=verbose,
    )

    # lower side
    cdc_lower = compute_single_side_drag_coefficient_green(
        lower_steps,
        exit_radius,
        operating_point,
        boundary_layer_functions_lower,
        (;
            boundary_layer_options.lambda,
            boundary_layer_options.longitudinal_curvature,
            boundary_layer_options.lateral_strain,
            boundary_layer_options.dilation,
            boundary_layer_options.rk,
            separation_allowance=boundary_layer_options.separation_allowance_lower,
            separation_penalty=boundary_layer_options.separation_penalty_lower,
        );
        verbose=verbose,
    )

    # Return total viscous drag
    return total_viscous_drag_duct(
        cdc_upper, cdc_lower, exit_radius, operating_point.Vinf[], operating_point.rhoinf[]
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
