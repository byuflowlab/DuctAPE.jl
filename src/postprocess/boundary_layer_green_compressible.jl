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
    x_of_s = smooth_Akima(s, controlpoint[1, :])
    y_of_s = smooth_Akima(s, controlpoint[2, :])

    #get first and second derivatives of coordinates with respect to arc length at integration step locations
    xdot = FLOWMath.derivative.(Ref(x_of_s), ss)
    xdotsp = smooth_Akima(s, FLOWMath.derivative.(Ref(x_of_s), s))
    xddot = FLOWMath.derivative.(Ref(xdotsp), ss)

    ydot = FLOWMath.derivative.(Ref(y_of_s), ss)
    ydotsp = smooth_Akima(s, FLOWMath.derivative.(Ref(y_of_s), s))
    yddot = FLOWMath.derivative.(Ref(ydotsp), ss)

    # assemble the numerator and denominator of the radius of curvature expression
    num = (xdot .^ 2 .+ ydot .^ 2) .^ (1.5)
    den = xdot .* yddot .- ydot .* xddot

    # convex positive
    if abs(den) < eps()
        return 0.0
    else
        return num ./ den
    end
end

"""
    setup_boundary_layer_functions_green(
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
function setup_boundary_layer_functions_green(
    s,
    vtan_duct,
    duct_control_points,
    operating_point,
    boundary_layer_options;
    verbose=false,
)

    # Edge Velocities
    # note: the ss that get's passed in is the ss in the full surface length, so in practice this will be stagnation s ± boundary layer step s
    edge_velocity = smooth_Akima(s, vtan_duct)

    # Edge Accelerations (dUe/ds)
    edge_acceleration(ss) = FLOWMath.derivative.(Ref(edge_velocity), ss)

    # r's
    r_coords = smooth_Akima(s, duct_control_points[2, :])

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
    function Te(ss)
        return convert_temperature_to_kelvin.(
            Ref(operating_point.units),
            static_temperature.(Ref(operating_point.Ttot[]), edge_mach(ss)),
        )
    end
    function edge_viscosity(ss)
        return convert_viscosity.(Ref(operating_point.units), sutherlands_law(Te(ss)))
    end

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
    initialize_turbulent_boundary_layer_states_green(
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
function initialize_turbulent_boundary_layer_states_green(
    s_init, r_init, Ue, M, rhoe, mue; verbose=false
)
    verbose && println("INITIALIZATION")

    # println("s_init: ", s_init)
    # println("r_init: ", r_init)
    # println("Ue: ", Ue)
    # println("M: ", M)
    # println("rho: ", rhoe)
    # println("mue: ", mue)

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

    return (0.01013 / (log(10, FR(M) * Red2) - 1.02) - 0.00075) / Fc(M)
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
    calculate_beta(Ri; hardness=100.0)

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
function calculate_beta(Ri, hardness=100.0)
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
    return  H1 * (Cf / 2.0 - (H12 + 1.0) * d2dUedsUeeq0)
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
    c = Ctaueq0 / ((1.0 + 0.1 * M^2) * lambda^2) - 0.32 * Cf0
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
    boundary_layer_residual_green(y, s, parameters)

Out-of-place version of `boundary_layer_residual_green!`
"""
function boundary_layer_residual_green(y, s, parameters)
    dy = similar(y) .= 0
    return boundary_layer_residual_green!(dy, y, s, parameters)
end

"""
    boundary_layer_residual_green!(dy, y, s, parameters)

Calculate dy give the current states, y, the input position, s, and various parameters.
"""
function boundary_layer_residual_green!(dy, y, s, parameters)
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

    return dy, [Cf; H12]
end

"""
    solve_green_boundary_layer!(f, rk, initial_states, steps, parameters; verbose=false)

Integrate the turbulent boundary layer using a Runge-Kutta method.

# Arguments:
- `f::function_handle` : Governing residual equations to integrate
- `rk::function_handle` : Runge-Kutta method to use (RK2 or RK4)
- `initial_states::Float` : initial states
- `steps::Vector{Float}` : steps for integration
- `parameters::NamedTuple` : boundary layer solve options and other parameters
"""
function solve_green_boundary_layer!(
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
            println("  ds = ", abs(steps[i + 1] - steps[i]))
            println()
        end

        # take step
        us[:, i + 1], aux = rk(
            f, us[:, i], steps[i], abs(steps[i + 1] - steps[i]), parameters
        )

        Cfs[i + 1], H12s[i + 1] = aux

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

        s_sep = FLOWMath.linear(
            Cfs[(sepid[] - 1):sepid[]], steps[(sepid[] - 1):sepid[]], 0.0
        )
    else
        usep = us[:, end]
        H12sep = H12s[end]
        s_sep = steps[end]
    end

    # return states at separate, and separation shape factor, and surface length at separation
    return usep, H12sep, s_sep, sepid
end
