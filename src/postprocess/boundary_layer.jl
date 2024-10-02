#---------------------------------#
#         Setup functions         #
#---------------------------------#

"""
    split_at_stagnation_point(cp_duct)

Find the stagnation point as the point with local cp closest to 1.0.

# Arguments:
- `cp_duct::Vector{Float}` : Vector of surface pressure coefficients on duct from casing trailing edge clockwise to nacelle trailing edge.

# Returns:
- `upper_boundary_layer_indices::UnitRange{Int}` : range of panel indices for the upper boundary layer.
- `lower_boundary_layer_indices::StepRange{Int}` : range of panel indices for the lower boundary layer; in reverse order, and therefore front to back.
"""
function split_at_stagnation_point(cp_duct)
    _, stagnation_id = findmin(abs.(1.0 .- cp_duct))

    return stagnation_id:length(cp_duct), stagnation_id:-1:1
end

"""
    calc_radius_of_curvature(s, xy, ss)

Determine the radius of curvature.

(see https://en.wikipedia.org/wiki/Radius_of_curvature)

# Arugments:
- `s::Vector{Float}` : vector of surface lengths.
- `controlpoint::Matrix{Float} : control point positions, axial in first row, radial in second row
- `ss::Float` : position along surface length to find radius of curvature.

# Return:
- `radius_of_curvature::Float` : Radius of curvature at point `ss` along surface.
"""
function calc_radius_of_curvature(s, controlpoint, ss)

    # spline x and y coordinates with respect to arc length
    x_of_s = FLOWMath.Akima(s, controlpoint[1, :])
    y_of_s = FLOWMath.Akima(s, controlpoint[2, :])

    #get first and second derivatives of coordinates with respect to arc length at integration step locations
    xdot = FLOWMath.derivative.(Ref(x_of_s), ss)
    xddot = FLOWMath.second_derivative.(Ref(x_of_s), ss)
    ydot = FLOWMath.derivative.(Ref(y_of_s), ss)
    yddot = FLOWMath.second_derivative.(Ref(y_of_s), ss)

    # assemble the numerator and denominator of the radius of curvature expression
    num = (xdot .^ 2 .+ ydot .^ 2) .^ (1.5)
    den = xdot .* yddot .- ydot .* xddot

    # convex positive
    return num ./ den
end

"""
    bl_step_fun(n, m, p)

Function used in determining step sizes for boundary layer calculation. f(n) = m*n^p

Given a number of steps, n ∈ [1:N], provides the cumulative step lengths according to the power, p, and the multiplicative factor, m; where p determined from the `set_bl_steps` functions.
"""
function bl_step_fun(n, m, p)
    return m .* n .^ p
end

"""
    set_bl_steps(N::Int, first_step_size, total_length)

Sets boundary layer steps based on desired number of steps (must be an Integer), an initial step size, and the total cumulative length of the steps.

# Arguments:
- `N::Int` : Number of steps to take
- `first_step_size::Float` : size of first step (which is `m` in `bl_step_fun`)
- `total_length::Float` : total surface length to divide up.

# Returns:
- `steps::Vector{Float}` : steps along surface length satisfying the equation: f(n) = m*n^p with the condition that `m` is the first step size and f(N) = `total_length`
"""
function set_bl_steps(N::Int, first_step_size, total_length)
    function f(p)
        return total_length - first_step_size * N^p
    end
    p = find_zero(f, 1.0)

    return bl_step_fun(1:N, first_step_size, p)
end

"""
    arc_lengths_from_panel_lengths(panel_lengths, bl_ids)

Cumulative sum of panel lengths for the given section of surface associated with the upper or lower boundary layer.

# Arguments:
- `panel_lengths::Vector{Float}` : vector of panel lengths (called influence_length in body_vortex_panels).
- `bl_ids::UnitRange{Int}` : range of indices over which to determine the surface length

# Returns:
- `s::Vector{Float}` : cumulative sum of panel lengths between control points in the given index range, starting from zero.
"""
function arc_lengths_from_panel_lengths(panel_lengths, bl_ids)
    if bl_ids[1] < bl_ids[end]
        bl_range = (bl_ids[1] + 1):(bl_ids[end])

        return cumsum(
            [0.0; [0.5 * (panel_lengths[i] + panel_lengths[i - 1]) for i in bl_range]]
        )
    else
        bl_range = (bl_ids[1] - 1):-1:bl_ids[end]

        return cumsum(
            [0.0; [0.5 * (panel_lengths[i] + panel_lengths[i + 1]) for i in bl_range]]
        )
    end
end

"""
    setup_bl(
        vtan_duct,
        controlpoint,
        panel_lengths,
        bl_ids,
        operating_point,
        boundary_layer_options;
        verbose=false
    )

# Arguments:
- `panel_lengths::Vector{Float}` : vector of panel lengths for the duct body
- `controlpoint::Matrix{Float}` : Control point coordinates along the duct surface
- `panel_lengths::Matrix{Float}` : Panel lengths along the duct surface
- `bl_ids::UnitRange` : range of indices for which the
- `operating_point::OperatingPoint` : OperatingPoint object
- `boundary_layer_options::BoundaryLayerOptions` : BoundaryLayerOptions object

# Returns:
- `steps::Vector{Float}` : steps for ODE integration
- `boundary_layer_parameters::NamedTuple` : Namped Tuple containing boundary layer solver parameters:
  - `edge_velocity::FLOWMath.Akima` : spline of edge velocities relative to surface length
  - `edge_mach::FLOWMath.Akima` : spline of edge Mach numbers relative to surface length
  - `edge_acceleration::FLOWMath.Akima` : spline of edge acceleration (dUe/ds) relative to surface length
  - `edge_density::FLOWMath.Akima` : spline of edge density relative to surface length
  - `edge_viscosity::FLOWMath.Akima` : spline of edge viscosity relative to surface length
  - `r_coords::FLOWMath.Akima` : spline radial coordaintes relative to surface length
  - `radial_derivative::FLOWMath.Akima` : spline of dr/ds relative to surface length
  - `radius_of_curvature::FLOWMath.Akima` : spline of radius of curvature relative to surface length
"""
function setup_bl(
    vtan_duct,
    controlpoint,
    panel_lengths,
    bl_ids,
    operating_point,
    boundary_layer_options;
    verbose=false,
)

    # Arc lengths from panel lengths
    s = arc_lengths_from_panel_lengths(panel_lengths, bl_ids)

    # # Get integration steps
    # steps = range(
    #     boundary_layer_options.offset,
    #     s[end] - 0.5 * sum(panel_lengths[bl_ids[[1:1; end:end]]]) -
    #     boundary_layer_options.offset;
    #     step=boundary_layer_options.first_step_size,
    # )

    steps =
        set_bl_steps(
            boundary_layer_options.n_steps,
            boundary_layer_options.first_step_size,
            # total length = sum of all panel lengths - half the first and last panel lengths - offset for starting point
            sum(panel_lengths[bl_ids]) - boundary_layer_options.offset -
            0.5 * sum(panel_lengths[bl_ids[[1:1; end:end]]]),
        ) .+ boundary_layer_options.offset

    # Edge Velocities
    edge_velocity = FLOWMath.Akima(s, vtan_duct[bl_ids])

    # Edge Accelerations (dUe/ds)
    edge_acceleration(ss) = FLOWMath.derivative.(Ref(edge_velocity), ss)

    # r's
    r_coords = FLOWMath.Akima(s, controlpoint[2, bl_ids])

    # dr/ds
    radial_derivative(ss) = FLOWMath.derivative.(Ref(r_coords), ss)

    # radii of curvature
    radius_of_curvature(ss) = calc_radius_of_curvature(s, controlpoint[:, bl_ids], ss)

    # local mach number
    edge_mach(ss) = calc_mach(edge_velocity(ss), operating_point.asound[])

    # local density
    Pe(ss) = static_pressure(operating_point.Ptot[], edge_mach(ss))
    edge_density(ss) = static_density(Pe(ss), operating_point.asound[])

    # local viscosity
    Te(ss) = static_temperature.(operating_point.Ttot[], edge_mach(ss))
    edge_viscosity(ss) = sutherlands_law(Te(ss))

    return steps,
    (;
        edge_velocity,
        edge_mach,
        edge_acceleration,
        edge_density,
        edge_viscosity,
        r_coords,
        radial_derivative,
        radius_of_curvature,
        verbose,
    )
end

#---------------------------------#
#  state initilization functions  #
#---------------------------------#

"""
    d2_init(s_init, Rex)

Initialize momentum thickness state with flat plate Schlichting solution

# Arguments:
- `x::Float` : length
- `Rex::Float` : Reynolds number at the given length

# Returns:
- `d2::Float` : momentum thickness
"""
function d2_init(x, Rex)
    return 0.036 * x / Rex^0.2
end

"""
Rename of `calc_H12bar0` used for initialization to avoid confusion
"""
function H12bar_init(Cf0, M)
    return calc_H12bar0(Cf0, M)
end

"""
Wrapper of `calc_CEeq` used for initialization to avoid confusion (`lambda` set to 1)
"""
function CE_init(Ctaueq0, M, Cf0)
    return calc_CEeq(Ctaueq0, M, 1, Cf0)
end

"""
    initialize_turbulent_boundary_layer_states(
        s_init, r_init, Ue, M, rhoe, mue; verbose=false
    )

Initialize the states for the turbulent boundary layer solve.

# Arguments:
- `s_init::Float` : surface length starting point
- `r_init::Float` : surface radial position at starting point
- `Ue::Float` : edge velocity at starting point
- `M::Float` : edge Mach at starting point
- `rhoe::Float` : edge density at starting point
- `mue::Float` : edge dynamic viscosity at starting point

# Returns:
- `initial_states::Vector{Float}` : Initial values for the states: r*delta2, Hbar12, CE
- `Cf_init::Vector{Float}` : initial value for friction coefficietn Cf (used for determining separateion)
- `H12_init::Vector{Float}` : initial value for shape factor H12 (used in viscous drag calculation)
"""
function initialize_turbulent_boundary_layer_states(
    s_init, r_init, Ue, M, rhoe, mue; verbose=false
)
    verbose && println("INITIALIZATION")
    # - Initialize States - #
    # initialize momentum thickness (d2) using flat plate schlichting model
    Rex = calc_Re(rhoe, Ue, s_init, mue)
    verbose && printdebug("Rex: ", Rex)
    d2 = d2_init(s_init, Rex)
    verbose && printdebug("d2: ", d2)

    # initialize H12bar (compressible shape factor) using _0 equations
    Red2 = calc_Re(rhoe, Ue, d2, mue)
    verbose && printdebug("Red2: ", Red2)
    Cf0 = calc_Cf0(Red2, M)
    verbose && printdebug("Cf0: ", Cf0)
    H12bar0 = H12bar_init(Cf0, M)
    verbose && printdebug("H12bar0: ", H12bar0)

    # initialize the entrainment coefficient (CE) using equilibrium equations
    H1 = calc_H1(H12bar0)
    verbose && printdebug("H1: ", H1)
    H12 = calc_H12(H12bar0, M)
    verbose && printdebug("H12: ", H12)
    Cf = calc_Cf(H12bar0, H12bar0, Cf0)
    verbose && printdebug("Cf: ", Cf)
    d2dUedsUeeq0 = calc_d2dUedsUeeq0(H12bar0, H12, Cf, M)
    verbose && printdebug("d2dUedsUeeq0: ", d2dUedsUeeq0)
    CEeq0 = calc_CEeq0(H1, H12, Cf, d2dUedsUeeq0)
    verbose && printdebug("CEeq0: ", CEeq0)
    Ctaueq0 = calc_Ctaueq0(CEeq0, Cf0, M)
    verbose && printdebug("Ctaueq0: ", Ctaueq0)
    CE = CE_init(Ctaueq0, M, Cf0)
    verbose && printdebug("CE: ", CE)

    # initial states
    return [r_init * d2, H12bar0, CE], Cf, H12
end

#---------------------------------#
#     Intermediate functions      #
#---------------------------------#

"""
    Fc(M) = sqrt(1.0 + 0.2 * M^2)
"""
function Fc(M)
    return sqrt(1.0 + 0.2 * M^2)
end

"""
    FR(M) = 1.0 + 0.056 * M^2
"""
function FR(M)
    return 1.0 + 0.056 * M^2
end

"""
    calc_Re(ρ, V, L, μ)

Calculate Reynolds number
"""
function calc_Re(rhoe, Ue, d2, mue)
    return rhoe * Ue * d2 / mue
end

"""
    calc_Cf0(Red2, M; hardness=50)

Calculate Cf₀ value based on momentum thickness Reynolds number (Red2) and edge mach number (M).
"""
function calc_Cf0(Red2, M; hardness=50)
    return (0.01013 / (log(10, FR(M) * Red2) - 1.02) - 0.00075) / Fc(M)
end

"""
    calc_H12bar0(Cf0, M)
"""
function calc_H12bar0(Cf0, M)
    return 1.0 / (1.0 - 6.55 * sqrt(Cf0 / 2.0 * (1.0 + 0.04 * M^2)))
end

"""
    calc_Cf(H12bar, H12bar0, Cf0)

Calculate friction coefficient.
"""
function calc_Cf(H12bar, H12bar0, Cf0)
    return Cf0 * (0.9 / (H12bar / H12bar0 - 0.4) - 0.5)
end

"""
    calc_H12(H12bar, M, Pr=1.0)
"""
function calc_H12(H12bar, M, Pr=1.0)
    return (H12bar + 1.0) * (1.0 + Pr^(1.0 / 3.0) * M^2 / 5) - 1.0
end

"""
    calc_Ctau(CE, Cf0, M)
"""
function calc_Ctau(CE, Cf0, M)
    return (0.024 * CE + 1.2 * CE^2 + 0.32 * Cf0) * (1.0 + 0.1 * M^2)
end

"""
    calc_F(CE, Cf0)
"""
function calc_F(CE, Cf0)
    return (0.02 * CE + CE^2 + 0.8 * Cf0 / 3) / (0.01 + CE)
end

"""
    calc_H1(H12bar)
"""
function calc_H1(H12bar)
    return 3.15 + 1.72 / (H12bar - 1.0) - 0.01 * (H12bar - 1.0)^2
end

"""
    calc_dH12bardH1(H12bar)
"""
function calc_dH12bardH1(H12bar)
    return -(H12bar - 1.0)^2 / (1.72 + 0.02 * (H12bar - 1.0)^3)
end

"""
    calc_richardson_number(H12bar, d2, H12, H1, R)
"""
function calc_richardson_number(H12bar, d2, H12, H1, R)
    return 2.0 * d2 / (3.0 * R) * (H12 + H1) * (H1 / H12bar + 0.3)
end

"""
    beta(Ri; hardness=100.0)

Sigmoind blended version of piecewise function:
    | 7.0 if Ri > 0
β = |
    | 4.5 if Ri < 0

# Arguments:
- `Ri::float` : Richardson number

# Keyword Arguments:
- `hardness::float` : hardness factor for sigmoid blend

# Returns:
- `beta::float` : factor used in secondary influence from longitudinal curvature.
"""
function beta(Ri, hardness=100.0)
    return FLOWMath.sigmoid_blend(4.5, 7.0, Ri, 0.0, hardness)
end

"""
    longitudinal_curvature_influence(M, Ri)
"""
function longitudinal_curvature_influence(M, Ri)
    return 1 + beta(Ri) * (1.0 + M^2 / 5) * Ri
end

"""
    lateral_strain_influence(H12bar, d2, H12, H1, r, drds)
"""
function lateral_strain_influence(H12bar, d2, H12, H1, r, drds)
    return 1.0 - 7.0 / 3.0 * (H1 / H12bar + 1.0) * (H12 + H1) * d2 * drds / r
end

"""
    dilation_influence(H12bar, d2, H12, H1, M, Ue, dUeds)
"""
function dilation_influence(H12bar, d2, H12, H1, M, Ue, dUeds)
    return 1.0 + 7.0 / 3.0 * M^2 * (H1 / H12bar + 1.0) * (H12 + H1) * d2 * dUeds / Ue
end

"""
    calc_lambda(args...)

λ = *(args...)
"""
function calc_lambda(args...)
    return *(args...)
end

"""
    calc_d2dUedsUeeq0(H12bar, H12, Cf, M)
"""
function calc_d2dUedsUeeq0(H12bar, H12, Cf, M)
    return 1.25 * (Cf / 2.0 - ((H12bar - 1.0) / (6.432 * H12bar))^2 / (1.0 + 0.04 * M^2)) /
           H12
end

"""
    calc_CEeq0(H1, H12, Cf, d2dUedsUeeq0)

Note: applies smooth-max to makes sure this value stays non-negative
"""
function calc_CEeq0(H1, H12, Cf, d2dUedsUeeq0)
    # apply smooth-max to make sure we stay non-negative.
    return FLOWMath.ksmax([0.0; H1 * (Cf / 2.0 - (H12 + 1.0) * d2dUedsUeeq0)])
end

"""
    calc_Ctaueq0(CEeq0, Cf0, M)
"""
function calc_Ctaueq0(CEeq0, Cf0, M)
    return (0.24 * CEeq0 + 1.2 * CEeq0^2 + 0.32 * Cf0) * (1.0 + 0.1 * M^2)
end

"""
    calc_CEeq(Ctaueq0, M, lambda, Cf0)
"""
function calc_CEeq(Ctaueq0, M, lambda, Cf0)
    c = Ctaueq0 / ((1.0 + 0.1 * M^2) * lambda^2) - 0.32 * Cf0
    return sqrt(c / 1.2 + 0.0001) - 0.01
end

"""
    calc_d2dUedsUeeq(H1, H12, Cf, CEeq)
"""
function calc_d2dUedsUeeq(H1, H12, Cf, CEeq)
    return (Cf / 2.0 - CEeq / H1) / (H12 + 1.0)
end

#---------------------------------#
#       Residual functions        #
#---------------------------------#

"""
    boundary_layer_residual(y, s, parameters)

Out-of-place version of `boundary_layer_residual!`
"""
function boundary_layer_residual(y, s, parameters)
    dy = zeros(length(y))
    return boundary_layer_residual!(dy, y, s, parameters)
end

"""
    boundary_layer_residual!(dy, y, s, parameters)

Calculate dy give the current states, y, the input position, s, and various parameters.
"""
function boundary_layer_residual!(dy, y, s, parameters)
    (;
        r_coords,
        edge_velocity,
        edge_density,
        edge_viscosity,
        edge_mach,
        edge_acceleration,
        radius_of_curvature,
        radial_derivative,
        verbose,
    ) = parameters

    # - rename for convenience - #
    verbose && println("   Inputs:")
    r = r_coords(s)
    verbose && printdebug("r: ", r)
    Ue = edge_velocity(s)
    verbose && printdebug("Ue: ", Ue)
    rhoe = edge_density(s)
    verbose && printdebug("rhoe: ", rhoe)
    mue = edge_viscosity(s)
    verbose && printdebug("mue: ", mue)
    M = edge_mach(s)
    verbose && printdebug("M: ", M)
    dUeds = edge_acceleration(s)
    verbose && printdebug("dUeds: ", dUeds)
    R = radius_of_curvature(s)
    verbose && printdebug("R: ", R)
    drds = radial_derivative(s)
    verbose && printdebug("drds: ", drds)

    # - unpack variables - #
    rd2, H12bar, CE = y
    verbose && println("   States:")
    verbose && printdebug("rd2: ", rd2)
    verbose && printdebug("H12bar: ", H12bar)
    verbose && printdebug("CE: ", CE)

    # - calculate intermediate variables - #
    verbose && println("   Intermediates:")
    d2 = rd2 / r
    verbose && printdebug("d2: ", d2)
    Red2 = calc_Re(rhoe, Ue, d2, mue)
    verbose && printdebug("Red2: ", Red2)

    Cf0 = calc_Cf0(Red2, M)
    verbose && printdebug("Cf0: ", Cf0)

    H12bar0 = calc_H12bar0(Cf0, M)
    verbose && printdebug("H12bar0: ", H12bar0)

    Cf = calc_Cf(H12bar, H12bar0, Cf0)
    verbose && printdebug("Cf: ", Cf)

    H12 = calc_H12(H12bar, M)
    verbose && printdebug("H12: ", H12)

    dH12bardH1 = calc_dH12bardH1(H12bar) #yes, this should be negative according to example
    verbose && printdebug("dH12bardH1: ", dH12bardH1)
    H1 = calc_H1(H12bar)
    verbose && printdebug("H1: ", H1)

    F = calc_F(CE, Cf0)
    verbose && printdebug("F: ", F)

    d2dUedsUeeq0 = calc_d2dUedsUeeq0(H12bar, H12, Cf, M)
    verbose && printdebug("d2dUedsUeeq0: ", d2dUedsUeeq0)

    CEeq0 = calc_CEeq0(H1, H12, Cf, d2dUedsUeeq0)
    verbose && printdebug("CEeq0: ", CEeq0)

    Ctaueq0 = calc_Ctaueq0(CEeq0, Cf0, M)
    verbose && printdebug("Ctaueq0: ", Ctaueq0)

    Ctau = calc_Ctau(CE, Cf0, M)
    verbose && printdebug("Ctau: ", Ctau)

    if parameters.lambda
        if parameters.longitudinal_curvature
            Ri = calc_richardson_number(H12bar, d2, H12, H1, r)
            l1 = longitudinal_curvature_influence(M, Ri)
        else
            l1 = 1
        end

        if parameters.lateral_strain
            l2 = lateral_strain_influence(H12bar, d2, H12, H1, r, drds)
        else
            l2 = 1
        end

        if parameters.dilation
            l3 = dilation_influence(H12bar, d2, H12, H1, M, Ue, dUeds)
        else
            l3 = 1
        end

        lambda = calc_lambda(l1, l2, l3)
    else
        lambda = 1
    end
    verbose && printdebug("lambda: ", lambda)

    CEeq = calc_CEeq(Ctaueq0, M, lambda, Cf0)
    verbose && printdebug("CEeq: ", CEeq)

    d2dUedsUeeq = calc_d2dUedsUeeq(H1, H12, Cf, CEeq)
    verbose && printdebug("d2dUedsUeeq: ", d2dUedsUeeq)

    # - system of equations - #

    # momentum
    dy[1] = r * Cf / 2.0 - (H12 + 2.0 - M^2) * rd2 * dUeds / Ue

    # entrainment
    dy[2] = dH12bardH1 * (CE - H1 * (Cf / 2.0 - (H12 + 1.0) * d2 * dUeds / Ue)) / d2

    # lag
    dy[3] =
        F / d2 * (
            2.8 / (H12 + H1) * (sqrt(Ctaueq0) - lambda * sqrt(Ctau)) + d2dUedsUeeq -
            d2 * dUeds / Ue * (1.0 + 0.075 * M^2 * ((1.0 + 0.2 * M^2) / (1.0 + 0.1 * M^2)))
        )
    verbose && println("   Residuals:")
    verbose && printdebug("d(rd2)/ds: ", dy[1])
    verbose && printdebug("d(H12bar)/ds: ", dy[2])
    verbose && printdebug("d(CE)/ds: ", dy[3])

    return dy, Cf, H12
end

#---------------------------------#
#        Solver functions         #
#---------------------------------#

"""
    RK2(f, y, s, ds, parameters)

2nd Order Runge-Kutta integration scheme.

# Arguments:
- `f::function_handle` : residual function for integration
- `y::Vector{Float}` : current states
- `s::Float` : current position
- `ds::Float` : step size
- `parameters::NamedTuple` : BoundaryLayerOptions and other various parameters
"""
function RK2(f, y, s, ds, parameters)
    parameters.verbose && println("  1st call:")
    k1, _, _ = f(y, s, parameters)
    parameters.verbose && println("  2nd call:")
    k2, Cfnext, H12next = f(y .+ (ds / 2) .* k1, s + (ds / 2), parameters)
    unext = @. y + k2 * ds
    return unext, Cfnext, H12next
end

"""
    RK4(f, y, s, ds, parameters)

4th Order Runge-Kutta integration scheme.

# Arguments:
- `f::function_handle` : residual function for integration
- `y::Vector{Float}` : current states
- `s::Float` : current position
- `ds::Float` : step size
- `parameters::NamedTuple` : BoundaryLayerOptions and other various parameters
"""
function RK4(f, y, s, ds, parameters)
    k1, Cf1, H12_1 = f(y, s, parameters)
    k2, Cf2, H12_2 = f(y .+ (ds / 2) .* k1, s + (ds / 2), parameters)
    k3, Cf3, H12_3 = f(y .+ (ds / 2) .* k2, s + (ds / 2), parameters)
    k4, Cf4, H12_4 = f(y .+ ds .* k3, s + ds, parameters)
    unext = @. y + (k1 + k2 * 2 + k3 * 2 + k4) * ds / 6
    Cfnext = (Cf1 + Cf2 * 2 + Cf3 * 2 + Cf4) / 6
    H12next = (H12_1 + H12_2 * 2 + H12_3 * 2 + H12_4) / 6
    return unext, Cfnext, H12next
end

"""
    solve_turbulent_boundary_layer_rk!(f, rk, rk, initial_states, steps, parameters; verbose=false)

Integrate the turbulent boundary layer using a Runge-Kutta method.

# Arguments:
- `f::function_handle` : Governing residual equations to integrate
- `rk::function_handle` : Runge-Kutta method to use (RK2 or RK4)
- `initial_states::Float` : initial states
- `steps::Vector{Float}` : steps for integration
- `parameters::NamedTuple` : boundary layer solve options and other parameters
"""
function solve_turbulent_boundary_layer_rk!(f, rk, initial_states, steps, parameters; verbose=false)

    # Unpack States and variables for viscous drag
    u0, Cf0, H12_0 = rk, initial_states

    # Initilization separate flags and outputs
    sep = [false]
    sepid = [1]

    # Allocate intermediate states and outputs
    us = zeros(eltype(u0), length(u0), length(steps))
    Cfs = zeros(eltype(u0), length(steps))
    H12s = zeros(eltype(u0), length(steps))

    # Initialize intermediate states and outputs
    us[:, 1] = u0
    Cfs[1] = Cf0
    H12s[1] = H12_0

    # Take Runge-Kutta Steps until either the end of the surface or separation is detected
    for i in 1:(length(steps) - 1)
        if verbose
            println()
            println(i, " of $(length(steps)-1)")
            println("  s = ", steps[i])
            println("  ds = ", steps[i + 1] - steps[i])
            println()
        end

        # take step
        us[:, i + 1], Cfs[i + 1], H12s[i + 1] = rk(
            f, us[:, i], steps[i], steps[i + 1] - steps[i], parameters
        )

        sepid[1] = i
        if Cfs[i + 1] <= 0.0
            sep[1] = true
            break
        end
    end

    if sep[1] == true
        u1sep = FLOWMath.akima(
            Cfs[(sepid[] - 1):sepid[]], us[1, (sepid[] - 1):sepid[]], 0.0
        )
        u2sep = FLOWMath.akima(
            Cfs[(sepid[] - 1):sepid[]], us[2, (sepid[] - 1):sepid[]], 0.0
        )
        u3sep = FLOWMath.akima(
            Cfs[(sepid[] - 1):sepid[]], us[3, (sepid[] - 1):sepid[]], 0.0
        )
        usep = [u1sep; u2sep; u3sep]
        H12sep = FLOWMath.akima(
            Cfs[(sepid[] - 1):sepid[]], H12s[(sepid[] - 1):sepid[]], 0.0
        )
        ssep = FLOWMath.akima(Cfs[(sepid[] - 1):sepid[]], steps[(sepid[] - 1):sepid[]], 0.0)
    else
        usep = us[:, end]
        H12sep = H12s[end]
        ssep = steps[end]
    end

    # return states at separate, and separation state, and separation postition
    return us, usep, Cfs, H12sep, ssep, sepid[1]
end

# """
# """
# function solve_turbulent_boundary_layer(f!, u0, xspan, parameters)
#     prob = DifferentialEquations.ODEProblem(f!, u0, xspan, parameters)
#     # alg = DifferentialEquations.radau()
#     alg = DifferentialEquations.Tsit5()
#     sol = DifferentialEquations.solve(prob, alg)
#     return sol
# end

#---------------------------------#
#          Viscous Drag           #
#---------------------------------#

# TODO: figure out how to get dimensional drag from squire young.
# should we be dividing by c? (probably)
"""
"""
function squire_young(d2, Ue, Uinf, H12, c)
    return 2.0 * d2 / c * (Ue / Uinf)^((5.0 + H12) / 2.0)
end

#---------------------------------#
#        Overall Functions        #
#---------------------------------#

"""
"""
function compute_viscous_drag_duct()

    # Find Stagnation Point
    upper_bl_ids, lower_bl_ids = split_at_stagnation_point(cp_duct)

    # - Upper Side - #
    us, usep, Cfs, H12sep, ssep, sepid = compute_viscous_drag_single_side(
        vtan_duct,
        controlpoint,
        influence_length,
        upper_bl_ids,
        operating_point,
        boundary_layer_options,
    )

    # - Lower Side - #
    us, usep, Cfs, H12sep, ssep, sepid = compute_viscous_drag_single_side(
        vtan_duct,
        controlpoint,
        influence_length,
        lower_bl_ids,
        operating_point,
        boundary_layer_options,
    )

    return nothing
end

"""
"""
function compute_viscous_drag_single_side(
    vtan_duct,
    controlpoint,
    influence_length,
    bl_ids,
    operating_point,
    boundary_layer_options;
    verbose=false,
)

    # - Set up boundary layer solve parameters - #
    steps, boundary_layer_parameters = setup_bl(
        vtan_duct,
        controlpoint,
        influence_length,
        bl_ids,
        operating_point,
        boundary_layer_options;
        verbose=verbose,
    )

    # - Initialize Boundary Layer States - #
    initial_states, Cf0, H12_0 = initialize_turbulent_boundary_layer_states
    verbose = false(
        steps[1],
        boundary_layer_parameters.r_coords(steps[1]),
        boundary_layer_parameters.edge_velocity(steps[1]),
        boundary_layer_parameters.edge_mach(steps[1]),
        boundary_layer_parameters.edge_density(steps[1]),
        boundary_layer_parameters.edge_viscosity(steps[1]);
        verbose=verbose,
    )

    return solve_turbulent_boundary_layer_rk!(
        boundary_layer_residual,
        boundary_layer_options.rk,
        [initial_states, Cf0, H12_0],
        steps,
        parameters,
    )
end
