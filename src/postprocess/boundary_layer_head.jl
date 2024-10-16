"""
    setup_boundary_layer_functions_head(
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
  - `edge_acceleration::FLOWMath.Akima` : spline of edge acceleration (dUe/ds) relative to surface length
  - `edge_density::FLOWMath.Akima` : spline of edge density relative to surface length
  - `edge_viscosity::FLOWMath.Akima` : spline of edge viscosity relative to surface length
"""
function setup_boundary_layer_functions_head(
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

    # local density
    edge_mach(ss) = calculate_mach(edge_velocity(ss), operating_point.asound[])
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

    # r's
    r_coords = Akima_smooth(s, duct_control_points[2, :])

    return (;
        edge_velocity, edge_acceleration, edge_density, edge_viscosity, r_coords, verbose
    )
end

"""
"""
function calculate_H(H1)

    # get each side of the piecewise equation
    hgeq = 0.86 * (H1 - 3.3)^(-0.777) + 1.1
    hlt = 1.1538 * (H1 - 3.3)^(-0.325) + 0.6778

    # blend the pieces smoothly
    return FLOWMath.sigmoid_blend(hlt, hgeq, H1, H1)
end

"""
"""
function calculate_cf(H, Red2)
    return 0.245 * 10^(-0.678 * H) * Red2^(-0.268)
end

"""
"""
function boundary_layer_residual_head(y, s, parameters)
    dy = similar(y) .= 0
    return boundary_layer_residual_head!(dy, y, s, parameters)
end

"""
    boundary_layer_residual_head!(dy, y, s, parameters)

Calculate dy give the current states, y, the input position, s, and various parameters.
"""
function boundary_layer_residual_head!(dy, y, s, parameters)

    # - unpack parameters - #
    (; verbose) = parameters

    # - unpack variables - #
    d2, H1 = y
    verbose && printdebug("d2: ", d2)
    verbose && printdebug("H1: ", H1)

    # limit H1 to be greater than 3.3
    H1lim = FLOWMath.ksmax([H1; 3.3 + 1e-2])
    verbose && printdebug("H1lim: ", H1lim)

    # - unpack variables - #
    (; edge_velocity, edge_acceleration, edge_density, edge_viscosity) = parameters

    # - Intermediate Calculations - #
    # determine dUedx
    dUedx = edge_acceleration(s)
    verbose && printdebug("dUedx: ", dUedx)

    # determine H
    H = calculate_H(H1lim)
    verbose && printdebug("H: ", H)

    # determine local edge velocity
    Ue = edge_velocity(s)
    verbose && printdebug("Ue: ", Ue)

    # determine momentum thickness reynolds number
    Red2 = calculate_Re(edge_density(s), Ue, d2, edge_viscosity(s))
    verbose && printdebug("Red2: ", Red2)

    # determine cf
    cf = calculate_cf(H, Red2)
    verbose && printdebug("cf: ", cf)

    # - "Residuals" - #
    # dδ₂/dx
    dy[1] = cf / 2.0 - dUedx * d2 / Ue * (H + 2.0)
    verbose && printdebug("dy[1]: ", dy[1])

    # dH₁/dx
    dy[2] = 0.0306 / d2 * (H1lim - 3.0)^(-0.6169) - dUedx * H1lim / Ue - dy[1] * H1lim / d2
    verbose && printdebug("dy[2]: ", dy[2])

    return dy, H
end

"""
    solve_head_boundary_layer!(f, rk, initial_states, steps, parameters; verbose=false)

Integrate the turbulent boundary layer using a Runge-Kutta method.

# Arguments:
- `f::function_handle` : Governing residual equations to integrate
- `rk::function_handle` : Runge-Kutta method to use (RK2 or RK4)
- `initial_states::Float` : initial states
- `steps::Vector{Float}` : steps for integration
- `parameters::NamedTuple` : boundary layer solve options and other parameters
"""
function solve_head_boundary_layer!(f, rk, initial_states, steps, parameters; verbose=false)

    # Unpack States and variables for viscous drag
    u0, H0 = initial_states

    # Initilization separate flags and outputs
    sep = [false]
    sepid = [1]

    # Allocate intermediate states and outputs
    # TODO; put these in a cache that gets updated in place.
    us = zeros(eltype(u0), length(u0), length(steps))
    Hs = zeros(eltype(u0), length(steps))

    # Initialize intermediate states and outputs
    us[:, 1] = u0
    Hs[1] = H0

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
        us[:, i + 1], Hs[i + 1] = rk(
            f, us[:, i], steps[i], abs(steps[i + 1] - steps[i]), parameters
        )

        sepid[1] = i + 1
        if Hs[i + 1] >= 3.0
            sep[1] = true
            break
        end
    end

    if sep[1] == true

        # - Interpolate to find actual s_sep - #

        usep = [
            FLOWMath.linear(Hs[(sepid[] - 1):sepid[]], us[1, (sepid[] - 1):sepid[]], 3.0)
            FLOWMath.linear(Hs[(sepid[] - 1):sepid[]], us[2, (sepid[] - 1):sepid[]], 3.0)
        ]

        Hsep = FLOWMath.linear(Hs[(sepid[] - 1):sepid[]], Hs[(sepid[] - 1):sepid[]], 3.0)

        s_sep = FLOWMath.linear(
            Hs[(sepid[] - 1):sepid[]], steps[(sepid[] - 1):sepid[]], 3.0
        )
    else
        usep = us[:, end]
        Hsep = Hs[end]
        s_sep = steps[end]
    end

    # return states at separate, and separation shape factor, and surface length at separation
    return usep, Hsep, s_sep, sepid
end
