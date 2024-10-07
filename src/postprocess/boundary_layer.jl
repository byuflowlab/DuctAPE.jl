#---------------------------------#
#         Setup functions         #
#---------------------------------#

"""
    arc_lengths_from_panel_lengths(duct_panel_lengths, bl_ids)

Cumulative sum of panel lengths for the given section of surface associated with the upper or lower boundary layer.

# Arguments:
- `duct_panel_lengths::Vector{Float}` : vector of panel lengths (called influence_length in body_vortex_panels) associated with the duct (casing + nacelle).

# Returns:
- `s::Vector{Float}` : cumulative sum of panel lengths between control points in the given index range, starting from zero.
"""
function arc_lengths_from_panel_lengths(duct_panel_lengths)
    return cumsum(
        [
            0.0
            [
                0.5 * (duct_panel_lengths[i] + duct_panel_lengths[i - 1]) for
                i in 2:length(duct_panel_lengths)
            ]
        ],
    )
end

"""
   split_at_stagnation_point(duct_panel_lengths, cp_duct)

Find the stagnation point as the point with local cp closest to 1.0.

# Arguments:
- `duct_panel_lengths::Vector{Float}` : Vector of panel lengths for the duct from casing trailing edge clockwise to nacelle trailing edge.
- `cp_duct::Vector{Float}` : Vector of surface pressure coefficients on duct from casing trailing edge clockwise to nacelle trailing edge.

# Returns:
- `s::Vector{Float}` : cumulative sum of distance along panels between control points starting at 0.
- `s_stagnation::Float` : surface length at which the stagnation point occurs, calculated using extrapolation technique described below
- `lower_length::Float` : length from stagnation point to casing trailing edge control point
- `upper_length::Float` : length from stagnation point to nacelle trailing edge control point


To determine point of intersection between extrapolated lines, we do the following:

Determine the point-slope form of the lines extrapolated from either side of the panel with a pressure coefficient closest to 1.0, and set up equations to find a point somewhere along each line

intersection    Known points on    unknown    slopes of
   point         each line         factors    each line
    P1 =       [s[i-2]; cp[i-2]] + x1    *   [ds1; dcp1]
    P2 =       [s[i+1]; cp[i+1]] + x2    *   [ds2; dcp2]

where if we set P1=P2, in other words, to be the point the lines intersect, we can solve for the unknown factors that yield the point of intersection.

We therefore set P1=P2 and assemble a system of linear equations, then plug one of the solved facotrs back into its associated equation to obtain the s location of intersection of the extrapolated lines and take this to be the stagnation point, regardless of what the pressure coefficient ends up being.
"""
function split_at_stagnation_point(duct_panel_lengths, cp_duct)

    # get surface length along entire surface
    s = arc_lengths_from_panel_lengths(duct_panel_lengths)

    # find index of cp closest to 1.0
    _, spi = findmin(abs.(1.0 .- cp_duct))

    # - Extrapolate to determine smooth stagnation point along surface length - #
        # set up system of equations to determine intersection of extrapolated lines on either side of the index closest to cp = 1.0
    A = [
        s[spi - 1]-s[spi - 2] s[spi + 1]-s[spi + 2]
        cp_duct[spi - 1]-cp_duct[spi - 2] cp_duct[spi + 1]-cp_duct[spi + 2]
    ]
    b = [s[spi + 1] - s[spi - 2]; cp_duct[spi + 1] - cp_duct[spi - 2]]

    # note: this should not happen, but if it does, here's how to treat it.
    # if det(A) < eps()
    #     # singular matrix, lines are parallel
    #     # take stagnation point to be panel center
    #     s_stagnation = s[spi]
    # else
    # solve linear system for unknown factors
    x = ImplicitAD.linear_solve(A, b)

    # get stagnation point from one of the equations (choose one arbitrarily)
    s_stagnation = s[spi - 2] + x[1] * (s[spi - 1] - s[spi - 2])
    # end

    # lower length from zero to stagnation point
    lower_length = s_stagnation

    # upper length from end of surface back to stagnation point
    upper_length = maximum(s) - s_stagnation

    return s, s_stagnation, lower_length, upper_length
end

"""
    bl_step_fun(n, m, p)

Function used in determining step sizes for boundary layer calculation. f(n) = m*n^p

Given a number of steps, n ∈ [1:N], provides the cumulative step lengths according to the power, p, and the multiplicative factor, m; where p determined from the `set_boundary_layer_steps` functions.
"""
function bl_step_fun(n, m, p)
    return m .* n .^ p
end

"""
    set_boundary_layer_steps(N::Int, first_step_size, total_length)

Sets boundary layer steps based on desired number of steps (must be an Integer), an initial step size, and the total cumulative length of the steps.

# Arguments:
- `N::Int` : Number of steps to take
- `first_step_size::Float` : size of first step (which is `m` in `bl_step_fun`)
- `total_length::Float` : total surface length to divide up.

# Returns:
- `steps::Vector{Float}` : steps along surface length satisfying the equation: f(n) = m*n^p with the condition that `m` is the first step size and f(N) = `total_length`
"""
function set_boundary_layer_steps(N::Int, first_step_size, total_length)
    function f(p)
        return total_length - first_step_size * N^p
    end
    p = find_zero(f, 1.0)

    return bl_step_fun(1:N, first_step_size, p)
end

"""
    calculate_radius_of_curvature(s, xy, ss)

Determine the radius of curvature.

(see https://en.wikipedia.org/wiki/Radius_of_curvature)

# Arugments:
- `s::Vector{Float}` : vector of surface lengths.
- `controlpoint::Matrix{Float} : control point positions, axial in first row, radial in second row
- `ss::Float` : position along surface length to find radius of curvature.

# Return:
- `radius_of_curvature::Float` : Radius of curvature at point `ss` along surface.
"""
function calculate_radius_of_curvature(s, controlpoint, ss)

    # spline x and y coordinates with respect to arc length
    x_of_s = Akima_smooth(s, controlpoint[1, :])
    y_of_s = Akima_smooth(s, controlpoint[2, :])

    #get first and second derivatives of coordinates with respect to arc length at integration step locations
    xdot = FLOWMath.derivative.(Ref(x_of_s), ss)
    xdotsp = Akima_smooth(s,FLOWMath.derivative.(Ref(x_of_s), s))
    xddot = FLOWMath.derivative.(Ref(xdotsp), ss)
    # xddot = FLOWMath.second_derivative.(Ref(x_of_s), ss)
    ydot = FLOWMath.derivative.(Ref(y_of_s), ss)
    ydotsp = Akima_smooth(s,FLOWMath.derivative.(Ref(y_of_s), s))
    yddot = FLOWMath.derivative.(Ref(ydotsp), ss)
    # yddot = FLOWMath.second_derivative.(Ref(y_of_s), ss)

    # assemble the numerator and denominator of the radius of curvature expression
    num = (xdot .^ 2 .+ ydot .^ 2) .^ (1.5)
    den = xdot .* yddot .- ydot .* xddot

    # convex positive
    return num ./ den
end

"""
    setup_boundary_layer_functions(
        s,
        vtan_duct,
        duct_control_points,
        operating_point,
        boundary_layer_options;
        verbose=false
    )

# Arguments:
- `s::Vector{Float}` : cumulative sum of panel lengths between control points in the given index range, starting from zero.
- `vtan_duct::Vector{Float}` : tangential velocity magnitudes for the entire duct
- `duct_control_points::Matrix{Float}` : Control point coordinates along the duct surface
- `operating_point::OperatingPoint` : OperatingPoint object
- `boundary_layer_options::BoundaryLayerOptions` : BoundaryLayerOptions object

# Returns:
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
function setup_boundary_layer_functions(
    s,
    vtan_duct,
    duct_control_points,
    operating_point,
    boundary_layer_options;
    verbose=false,
)

    # Edge Velocities
    # note: the ss that get's passed in is the ss in the full surface length, so in practice this will be stagnation s ± boundary layer step s
    edge_velocity = Akima_smooth(s, vtan_duct)

    # Edge Accelerations (dUe/ds)
    edge_acceleration(ss) = FLOWMath.derivative.(Ref(edge_velocity), ss)

    # r's
    r_coords = Akima_smooth(s, duct_control_points[2, :])

    # dr/ds
    radial_derivative(ss) = FLOWMath.derivative.(Ref(r_coords), ss)

    # radii of curvature
    radius_of_curvature(ss) = calculate_radius_of_curvature(s, duct_control_points, ss)

    # local mach number
    edge_mach(ss) = calculate_mach(edge_velocity(ss), operating_point.asound[])

    # local density
    Pe(ss) = static_pressure(operating_point.Ptot[], edge_mach(ss))
    edge_density(ss) = static_density(Pe(ss), operating_point.asound[])

    # local viscosity
    Te(ss) = static_temperature.(operating_point.Ttot[], edge_mach(ss))
    edge_viscosity(ss) = sutherlands_law(Te(ss))

    return (;
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
Rename of `calculate_H12bar0` used for initialization to avoid confusion
"""
function H12bar_init(Cf0, M)
    return calculate_H12bar0(Cf0, M)
end

"""
Wrapper of `calculate_CEeq` used for initialization to avoid confusion (`lambda` set to 1)
"""
function CE_init(Ctaueq0, M, Cf0)
    return calculate_CEeq(Ctaueq0, M, 1, Cf0)
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
    Rex = calculate_Re(rhoe, Ue, s_init, mue)
    verbose && printdebug("Rex: ", Rex)
    d2 = d2_init(s_init, Rex)
    verbose && printdebug("d2: ", d2)

    # initialize H12bar (compressible shape factor) using _0 equations
    Red2 = calculate_Re(rhoe, Ue, d2, mue)
    verbose && printdebug("Red2: ", Red2)
    Cf0 = calculate_Cf0(Red2, M)
    verbose && printdebug("Cf0: ", Cf0)
    H12bar0 = H12bar_init(Cf0, M)
    verbose && printdebug("H12bar0: ", H12bar0)

    # initialize the entrainment coefficient (CE) using equilibrium equations
    H1 = calculate_H1(H12bar0)
    verbose && printdebug("H1: ", H1)
    H12 = calculate_H12(H12bar0, M)
    verbose && printdebug("H12: ", H12)
    Cf = calculate_Cf(H12bar0, H12bar0, Cf0)
    verbose && printdebug("Cf: ", Cf)
    d2dUedsUeeq0 = calculate_d2dUedsUeeq0(H12bar0, H12, Cf, M)
    verbose && printdebug("d2dUedsUeeq0: ", d2dUedsUeeq0)
    CEeq0 = calculate_CEeq0(H1, H12, Cf, d2dUedsUeeq0)
    verbose && printdebug("CEeq0: ", CEeq0)
    Ctaueq0 = calculate_Ctaueq0(CEeq0, Cf0, M)
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
    calculate_Re(ρ, V, L, μ)

Calculate Reynolds number
"""
function calculate_Re(rhoe, Ue, d2, mue)
    return rhoe * Ue * d2 / mue
end

"""
    calculate_Cf0(Red2, M)

Calculate Cf₀ value based on momentum thickness Reynolds number (Red2) and edge mach number (M).
"""
function calculate_Cf0(Red2, M)
    # if (0.01013 / (log(10, FR(M) * Red2) - 1.02) - 0.00075) / Fc(M) < 0.0
    # println("Red2: ", Red2)
    # println("FR(M): ", FR(M))
    # println("Log: ", log(10, FR(M) * Red2))
    # println("Log-1.02: ", log(10, FR(M) * Red2)-1.02)
    # println("Cf0: ", (0.01013 / (log(10, FR(M) * Red2) - 1.02) - 0.00075) / Fc(M))
# end
return FLOWMath.ksmax([0.0;(0.01013 / (log(10, FLOWMath.ksmax([1.0;FR(M) * Red2])) - 1.02) - 0.00075) / Fc(M)])
end

"""
    calculate_H12bar0(Cf0, M)
"""
function calculate_H12bar0(Cf0, M)
    return 1.0 / (1.0 - 6.55 * sqrt(Cf0 / 2.0 * (1.0 + 0.04 * M^2)))
end

"""
    calculate_Cf(H12bar, H12bar0, Cf0)

Calculate friction coefficient.
"""
function calculate_Cf(H12bar, H12bar0, Cf0)
    return Cf0 * (0.9 / (H12bar / H12bar0 - 0.4) - 0.5)
end

"""
    calculate_H12(H12bar, M, Pr=1.0)
"""
function calculate_H12(H12bar, M, Pr=1.0)
    return (H12bar + 1.0) * (1.0 + Pr^(1.0 / 3.0) * M^2 / 5) - 1.0
end

"""
    calculate_Ctau(CE, Cf0, M)
"""
function calculate_Ctau(CE, Cf0, M)
    return (0.024 * CE + 1.2 * CE^2 + 0.32 * Cf0) * (1.0 + 0.1 * M^2)
end

"""
    calculate_F(CE, Cf0)
"""
function calculate_F(CE, Cf0)
    return (0.02 * CE + CE^2 + 0.8 * Cf0 / 3) / (0.01 + CE)
end

"""
    calculate_H1(H12bar)
"""
function calculate_H1(H12bar)
    return 3.15 + 1.72 / (H12bar - 1.0) - 0.01 * (H12bar - 1.0)^2
end

"""
    calculate_dH12bardH1(H12bar)
"""
function calculate_dH12bardH1(H12bar)
    return -(H12bar - 1.0)^2 / (1.72 + 0.02 * (H12bar - 1.0)^3)
end

"""
    calculate_richardson_number(H12bar, d2, H12, H1, R)
"""
function calculate_richardson_number(H12bar, d2, H12, H1, R)
    return 2.0 * d2 / (3.0 * R) * (H12 + H1) * (H1 / H12bar + 0.3)
end

"""
    calcualte_beta(Ri; hardness=100.0)

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
function calcualte_beta(Ri, hardness=100.0)
    return FLOWMath.sigmoid_blend(4.5, 7.0, Ri, 0.0, hardness)
end

"""
    longitudinal_curvature_influence(M, Ri)
"""
function longitudinal_curvature_influence(M, Ri)
    return 1 + calculate_beta(Ri) * (1.0 + M^2 / 5) * Ri
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
    calculate_lambda(args...)

λ = *(args...)
"""
function calculate_lambda(args...)
    return *(args...)
end

"""
    calculate_d2dUedsUeeq0(H12bar, H12, Cf, M)
"""
function calculate_d2dUedsUeeq0(H12bar, H12, Cf, M)
    return 1.25 * (Cf / 2.0 - ((H12bar - 1.0) / (6.432 * H12bar))^2 / (1.0 + 0.04 * M^2)) /
           H12
end

"""
    calculate_CEeq0(H1, H12, Cf, d2dUedsUeeq0)

Note: applies smooth-max to makes sure this value stays non-negative
"""
function calculate_CEeq0(H1, H12, Cf, d2dUedsUeeq0)
    # apply smooth-max to make sure we stay non-negative.
    return FLOWMath.ksmax([0.0; H1 * (Cf / 2.0 - (H12 + 1.0) * d2dUedsUeeq0)])
end

"""
    calculate_Ctaueq0(CEeq0, Cf0, M)
"""
function calculate_Ctaueq0(CEeq0, Cf0, M)
    return (0.24 * CEeq0 + 1.2 * CEeq0^2 + 0.32 * Cf0) * (1.0 + 0.1 * M^2)
end

"""
    calculate_CEeq(Ctaueq0, M, lambda, Cf0)
"""
function calculate_CEeq(Ctaueq0, M, lambda, Cf0)
    c = FLOWMath.ksmax([0.0;Ctaueq0 / ((1.0 + 0.1 * M^2) * lambda^2) - 0.32 * Cf0])
    return sqrt(c / 1.2 + 0.0001) - 0.01
end

"""
    calculate_d2dUedsUeeq(H1, H12, Cf, CEeq)
"""
function calculate_d2dUedsUeeq(H1, H12, Cf, CEeq)
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
    Red2 = calculate_Re(rhoe, Ue, d2, mue)
    verbose && printdebug("Red2: ", Red2)

    Cf0 = calculate_Cf0(Red2, M)
    verbose && printdebug("Cf0: ", Cf0)

    H12bar0 = calculate_H12bar0(Cf0, M)
    verbose && printdebug("H12bar0: ", H12bar0)

    Cf = calculate_Cf(H12bar, H12bar0, Cf0)
    verbose && printdebug("Cf: ", Cf)

    H12 = calculate_H12(H12bar, M)
    verbose && printdebug("H12: ", H12)

    dH12bardH1 = calculate_dH12bardH1(H12bar) #yes, this should be negative according to example
    verbose && printdebug("dH12bardH1: ", dH12bardH1)
    H1 = calculate_H1(H12bar)
    verbose && printdebug("H1: ", H1)

    F = calculate_F(CE, Cf0)
    verbose && printdebug("F: ", F)

    d2dUedsUeeq0 = calculate_d2dUedsUeeq0(H12bar, H12, Cf, M)
    verbose && printdebug("d2dUedsUeeq0: ", d2dUedsUeeq0)

    CEeq0 = calculate_CEeq0(H1, H12, Cf, d2dUedsUeeq0)
    verbose && printdebug("CEeq0: ", CEeq0)

    Ctaueq0 = calculate_Ctaueq0(CEeq0, Cf0, M)
    verbose && printdebug("Ctaueq0: ", Ctaueq0)

    Ctau = calculate_Ctau(CE, Cf0, M)
    verbose && printdebug("Ctau: ", Ctau)

    if parameters.lambda
        if parameters.longitudinal_curvature
            Ri = calculate_richardson_number(H12bar, d2, H12, H1, r)
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

        lambda = calculate_lambda(l1, l2, l3)
    else
        lambda = 1
    end
    verbose && printdebug("lambda: ", lambda)

    CEeq = calculate_CEeq(Ctaueq0, M, lambda, Cf0)
    verbose && printdebug("CEeq: ", CEeq)

    d2dUedsUeeq = calculate_d2dUedsUeeq(H1, H12, Cf, CEeq)
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
    solve_turbulent_boundary_layer_rk!(f, rk, initial_states, steps, parameters; verbose=false)

Integrate the turbulent boundary layer using a Runge-Kutta method.

# Arguments:
- `f::function_handle` : Governing residual equations to integrate
- `rk::function_handle` : Runge-Kutta method to use (RK2 or RK4)
- `initial_states::Float` : initial states
- `steps::Vector{Float}` : steps for integration
- `parameters::NamedTuple` : boundary layer solve options and other parameters
"""
function solve_turbulent_boundary_layer_rk!(
    f, rk, initial_states, steps, parameters; verbose=false
)

    # Unpack States and variables for viscous drag
    u0, Cf0, H12_0 = initial_states

    # Initilization separate flags and outputs
    sep = [false]
    sepid = [1]

    # Allocate intermediate states and outputs
    # TODO; put these in a cache that gets updated in place.
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

        sepid[1] = i + 1
        if Cfs[i + 1] <= 0.0
            sep[1] = true
            break
        end
    end


    if sep[1] == true

        # - Interpolate to find actual s_sep - #

        usep = [
            FLOWMath.linear(Cfs[(sepid[] - 1):sepid[]], us[1, (sepid[] - 1):sepid[]], 0.0)
            FLOWMath.linear(Cfs[(sepid[] - 1):sepid[]], us[2, (sepid[] - 1):sepid[]], 0.0)
            FLOWMath.linear(Cfs[(sepid[] - 1):sepid[]], us[3, (sepid[] - 1):sepid[]], 0.0)
        ]

        H12sep = FLOWMath.linear(
            Cfs[(sepid[] - 1):sepid[]], H12s[(sepid[] - 1):sepid[]], 0.0
        )

        s_sep = FLOWMath.linear(Cfs[(sepid[] - 1):sepid[]], steps[(sepid[] - 1):sepid[]], 0.0)
    else
        usep = us[:, end]
        H12sep = H12s[end]
        s_sep = steps[end]
    end

    # return states at separate, and separation shape factor, and surface length at separation
    return usep, H12sep, s_sep
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
function total_viscous_drag_duct(cd_upper, cd_lower, exit_radius, Vref, rhoinf)
    # note: cd's are already times chord, so no need to have a separate chord variable

    # drag per unit length
    dprime = 0.5 * rhoinf * Vref^2 * (cd_upper + cd_lower)

    # drag of annular airfoil
    return dprime * 2.0 * pi * exit_radius
end

#---------------------------------#
#        Overall Functions        #
#---------------------------------#

"""
    compute_single_side_drag_coefficient(
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
function compute_single_side_drag_coefficient(
    steps,
    exit_radius,
    operating_point,
    boundary_layer_functions,
    boundary_layer_options;
    verbose=false,
)

    # - Initialize Boundary Layer States - #
    initial_states, Cf0, H12_0 = initialize_turbulent_boundary_layer_states(
        steps[1],
        boundary_layer_functions.r_coords(steps[1]),
        boundary_layer_functions.edge_velocity(steps[1]),
        boundary_layer_functions.edge_mach(steps[1]),
        boundary_layer_functions.edge_density(steps[1]),
        boundary_layer_functions.edge_viscosity(steps[1]);
        verbose=verbose,
    )

    usep, H12sep, s_sep = solve_turbulent_boundary_layer_rk!(
        boundary_layer_residual,
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

    cd = squire_young(
        abs(
            usep[1] / boundary_layer_functions.r_coords(s_sep) +
            boundary_layer_functions.r_coords(steps[end]) - exit_radius,
        ),
        boundary_layer_functions.edge_velocity(s_sep),
        operating_point.Vinf[],
        H12sep,
    )

    return cd
end

"""
    compute_viscous_drag_duct(
        vtan_duct,
        cp_duct,
        duct_control_points,
        duct_panel_lengths,
        exit_radius,
        operating_point,
        boundary_layer_options,
    )

Determine total, dimensional viscous drag on the duct.

# Arguments:
- `vtan_duct::Vector{Float}` : tangential velocity magnitudes for the entire duct
- `duct_control_points::Matrix{Float}` : control point positions for the entire duct
- `duct_panel_lengths::Vector{Float}` : panel lengths for the entire duct
- `exit_radius::Float` : radius at duct trailing edge (casing side)
- `operating_point::Float` : OperatingPoint object
- `boundary_layer_options::NamedTuple` : BoundaryLayerOptions object

# Returns:
- `duct_viscous_drag::Float` : total viscous drag of duct
"""
function compute_viscous_drag_duct(
    vtan_duct,
    cp_duct,
    duct_control_points,
    duct_panel_lengths,
    exit_radius,
    operating_point,
    boundary_layer_options;
    verbose=false,
)

    # find stagnation point
    s, s_stagnation, lower_length, upper_length = split_at_stagnation_point(
        duct_panel_lengths, cp_duct
    )

    # set up boundary layer solve parameters
    boundary_layer_functions = setup_boundary_layer_functions(
        s,
        vtan_duct,
        duct_control_points,
        operating_point,
        boundary_layer_options;
        verbose=verbose,
    )

    # - Set integration steps - #

    # upper side
    upper_steps =
        s_stagnation .+ set_boundary_layer_steps(
            boundary_layer_options.n_steps,
            boundary_layer_options.first_step_size,
            upper_length - boundary_layer_options.offset,
        ) .+ boundary_layer_options.offset

    # lower side
    lower_steps =
        s_stagnation .- set_boundary_layer_steps(
            boundary_layer_options.n_steps,
            boundary_layer_options.first_step_size,
            lower_length - boundary_layer_options.offset,
        ) .+ boundary_layer_options.offset

    # - Get drag coeffients - #

    # upper side
    cdc_upper = compute_single_side_drag_coefficient(
        upper_steps,
        exit_radius,
        operating_point,
        boundary_layer_functions,
        boundary_layer_options;
        verbose=verbose,
    )

    # lower side
    cdc_lower = compute_single_side_drag_coefficient(
        lower_steps,
        exit_radius,
        operating_point,
        boundary_layer_functions,
        boundary_layer_options;
        verbose=verbose,
    )

    # Return total viscous drag
    return total_viscous_drag_duct(
        cdc_upper, cdc_lower, exit_radius, operating_point.Vinf[], operating_point.rhoinf[]
    )
end
