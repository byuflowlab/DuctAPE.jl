"""
    setup_boundary_layer_functions_head(
        s,
        vtan_duct,
        duct_control_points,
        operating_point,
        boundary_layer_options;
        verbose=false
    )

# Arguments
- `s::Vector{Float}` : cumulative sum of panel lengths between control points in the given index range, starting from zero.
- `vtan_duct::Vector{Float}` : tangential velocity magnitudes for the entire duct
- `duct_control_points::Matrix{Float}` : Control point coordinates along the duct surface
- `operating_point::OperatingPoint` : OperatingPoint object
- `boundary_layer_options::BoundaryLayerOptions` : BoundaryLayerOptions object

# Returns
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
    edge_velocity = smooth_Akima(s, vtan_duct)

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

    # # r's
    # r_coords = smooth_Akima(s, duct_control_points[2, :])

    return (; edge_velocity, edge_acceleration, edge_density, edge_viscosity, verbose)
end

function regularization(fx, x; sigma=1e-4)
    if x < sigma
        return fx + (x - sigma)^2
    else
        return fx
    end
end

"""
    calculate_H(H1)

Calculate the value of the shape factor used in Head's method.
"""
function calculate_H(H1; eps=0.0)

    # get each side of the piecewise equation
    # hgeq = 0.86 * ((H1 - 3.3)^(-0.777)) + 1.1
    # hgeq = 0.86 * 1.0 / ((H1 - 3.3)^(0.777) + eps) + 1.1
    hgeq = 0.86 * 1.0 / (regularization((H1 - 3.3)^(0.777), H1; sigma=eps)) + 1.1
    # hlt = 1.1538 * ((H1 - 3.3)^(-0.326)) + 0.6778
    # hlt = 1.1538 * 1.0 / ((H1 - 3.3)^(0.326) + eps) + 0.6778
    hlt = 1.1538 * 1.0 / (regularization((H1 - 3.3)^(0.326), H1; sigma=eps)) + 0.6778
    # hlt = FLOWMath.ksmin([3.1; 1.1538 * (H1 - 3.3)^(-0.326) + 0.6778], 5)

    # blend the pieces smoothly
    return FLOWMath.sigmoid_blend(hlt, hgeq, H1, 5.3)
end

"""
    limH1(H1)

Returns a limited H1 to avoid undefined behavior
"""
function limH1(H1; eps=1e-4)
    return FLOWMath.ksmax([H1; 3.3 + eps], 25)
    # return H1
end

"""
    calculate_cf(H, Red2)

Calculate the skin friction coefficient used in Head's method
"""
function calculate_cf(H, Red2)
    return 0.246 * 10^(-0.678 * H) * Red2^(-0.268)
end

"""
    boundary_layer_residual_head(y, parameters, s)

Out of place residual function for Head's method.
"""
function boundary_layer_residual_head(y, parameters, s; debug=false)
    dy = similar(y) .= 0
    return boundary_layer_residual_head!(dy, y, parameters, s; debug=debug)
end

"""
    boundary_layer_residual_head!(dy, y, parameters, s)

Calculate dy given the current states, y, the input position, s, and various parameters.
"""
function boundary_layer_residual_head!(dy, y, parameters, s; debug=false)

    # - unpack parameters - #
    (;
        edge_velocity,
        edge_acceleration,
        edge_density,
        edge_viscosity,
        verbose,
        dy_eps,
        H1_eps,
        H_eps,
    ) = parameters

    # - unpack variables - #
    d2, H1 = y
    verbose && printdebug("d2: ", d2)
    verbose && printdebug("H1: ", H1)

    # limit H1 to be greater than 3.3(+1e-4)
    H1lim = limH1(H1; eps=H1_eps)
    verbose && printdebug("H1lim: ", H1lim)

    # - Intermediate Calculations - #

    # determine H
    H = calculate_H(H1lim; eps=H_eps)
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

    # determine dUedx
    dUedx = edge_acceleration(s)
    verbose && printdebug("dUedx: ", dUedx)

    # - "Residuals" - #
    # dδ₂/ds
    dy[1] = cf / 2.0 - dUedx * d2 / Ue * (H + 2.0)
    verbose && printdebug("dy[1]: ", dy[1])

    # dH₁/ds
    # dy[2] = 0.0306 / d2 * (H1lim - 3.0)^(-0.6169) - dUedx * H1lim / Ue - dy[1] * H1lim / d2
    # dy[2] =
    #     0.0306 / d2 * 1.0 / ((H1lim - 3.0)^(0.6169) + dy_eps) - dUedx * H1lim / Ue -
    #     dy[1] * H1lim / d2
    dy[2] =
        0.0306 / d2 * 1.0 / (regularization((H1lim - 3.0)^(0.6169), H1lim; sigma=dy_eps)) -
        dUedx * H1lim / Ue - dy[1] * H1lim / d2
    verbose && printdebug("dy[2]: ", dy[2])

    if debug
        return dy[1], dy[2], H1lim, H, Ue, dUedx, Red2, cf
    else
        return dy
    end
end

"""
    initialize_head_states(boundary_layer_functions, s_init; verbose=false)

Initialize the boundary layer state variables at the start of a surface coordinate.

# Arguments
- `boundary_layer_functions::Tuple`: A tuple containing functions `(edge_density, edge_velocity, edge_viscosity)` that return properties at a given coordinate.
- `s_init::Float64`: Initial position along the surface where initialization is performed.

# Keyword Arguments
- `verbose::Bool=false`: If true, enables verbose output for debugging.

# Returns
- `initial_states::Vector{Float64}`: Initial boundary layer state vector `[d20, H10]`, where:
  - `d20` is the initial boundary layer momentum thickness Reynolds number estimate.
  - `H10` is the initial shape factor guess.
- `H0::Float64`: Calculated boundary layer shape parameter `H` corresponding to `H10`.
"""
function initialize_head_states(boundary_layer_functions, s_init; verbose=false)
    (; edge_density, edge_velocity, edge_viscosity) = boundary_layer_functions

    # - Initialize Boundary Layer States - #
    H10 = 10.6
    d20 =
        0.036 * s_init /
        calculate_Re(
            edge_density(s_init), edge_velocity(s_init), s_init, edge_viscosity(s_init)
        )^0.2

    initial_states = [d20; H10]
    H0 = calculate_H(H10)

    return initial_states, H0
end

"""
    find_last_max_H(usol, stepsol)

Find the last maximum value of the shape factor `H` along a solution path and corresponding position.

# Arguments
- `usol::Matrix{Float64}`: Solution array where the second row corresponds to a variable related to `H` (e.g., `limH1`).
- `stepsol::Vector{Float64}`: Independent variable vector corresponding to positions along the solution path.

# Returns
- `usep::Vector{Float64}`: State vector at the last maximum of `H`.
- `Hsep::Float64`: Value of the shape factor `H` at the last maximum.
- `s_sep::Float64`: Position along the surface where the last maximum of `H` occurs.

# Notes
- Computes `H` from the solution using `calculate_H` and `limH1`.
- If `H` decreases at the end, uses a smoothed Akima spline to find the exact last maximum by identifying where the derivative changes sign.
- Uses root-finding via `Roots.find_zero` to find the zero derivative point precisely.
- If no sign change is found, the last point is taken as the maximum.
- The state vector at the maximum is interpolated/smoothed accordingly.
"""
function find_last_max_H(usol, stepsol)
    Hsol = calculate_H.(limH1.(usol[2, :]))
    # check if H drops at the end
    if Hsol[end] < Hsol[end - 1]
        # if it drops, spline H and s
        hsp = smooth_Akima(stepsol, Hsol)
        # determine where the sign of the derivative changes,
        zidx = findlast(
            x -> sign(x) != sign(FLOWMath.derivative(hsp, stepsol[end])),
            FLOWMath.derivative.(Ref(hsp), stepsol),
        )
        if isnothing(zidx)
            # if it doesn't change, return the final point.
            usep = usol[:, end]
            Hsep = Hsol[end]
            s_sep = stepsol[end]
        else
            # if it does change, zero find the zero derivative point.
            maxwrap(x) = FLOWMath.derivative.(Ref(hsp), x)
            bracket = [stepsol[max(1, zidx - 1)]; stepsol[min(zidx + 1, length(stepsol))]]
            if sign(maxwrap(bracket[1])) == sign(maxwrap(bracket[2]))
                s_sep = Roots.find_zero(maxwrap, stepsol[zidx])
            else
                s_sep = Roots.find_zero(maxwrap, bracket)
            end

            Hsep = hsp(s_sep)
            usep = [smooth_akima(stepsol, usol[1, :], s_sep); Hsep]
        end

    else
        usep = usol[:, end]
        Hsep = Hsol[end]
        s_sep = stepsol[end]
    end

    return usep, Hsep, s_sep
end

"""
    solve_head_boundary_layer!(f, ode, initial_states, steps, parameters; verbose=false)

Integrate the turbulent boundary layer using a Runge-Kutta method.

# Arguments
- `f::function_handle` : Governing residual equations to integrate
- `ode::function_handle` : ODE method to use (RK2 or RK4)
- `initial_states::Float` : initial states
- `steps::Vector{Float}` : steps for integration
- `parameters::NamedTuple` : boundary layer solve options and other parameters

# Returns
- `usep::Vector{Float64}`: State vector at the separation point or last step
- `Hsep::Float64`: Shape factor value at the separation or last maximum
- `s_sep::Float64`: Surface coordinate corresponding to separation or last maximum
- `us::Matrix{Float64}`: Solution states over all integration steps
- `steps::Vector{Float64}`: Integration step positions along the surface
"""
function solve_head_boundary_layer!(
    ::RK, ode, initial_states, steps, parameters; verbose=false
)
    f = boundary_layer_residual_head

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
        us[:, i + 1] = ode(f, us[:, i], steps[i], abs(steps[i + 1] - steps[i]), parameters)

        Hs[i + 1] = calculate_H(limH1(us[2, i + 1]))

        if Hs[i + 1] >= parameters.separation_criteria
            sepid[1] = i
            sep[1] = true
            if parameters.terminate
                break
            end
        end
    end

    Hsep = min(parameters.separation_criteria, maximum(Hs))

    if Hsep == parameters.separation_criteria
        # println("separated")

        # - Interpolate to find actual s_sep - #
        #spline H and steps from solution, get step at H=3
        sepwrap(x) = parameters.separation_criteria - smooth_akima(steps, Hs, x)
        hid = findfirst(x -> x > parameters.separation_criteria, Hs)
        s_sep = Roots.find_zero(sepwrap, [steps[hid - 1]; steps[hid]])
        # spline states and steps, get states at step for H=3
        usep = [smooth_akima(steps, u, s_sep) for u in eachrow(us)]

    else
        # println("not separated")
        if parameters.return_last_max_shape_factor
            # println("return last max")
            usep, Hsep, s_sep = find_last_max_H(us, steps)
        else
            # println("return last")
            usep = us[:, end]
            Hsep = calculate_H(limH1(usep[2]))
            s_sep = steps[end]
        end
    end

    # return states at separate, and separation shape factor, and surface length at separation
    if parameters.cutoff_Hsep
        Hsep = FLOWMath.ksmin([Hsep; parameters.separation_criteria], 50)
    end

    return usep, Hsep, s_sep, us, steps
end

"""
    update_Cf_head(state, step, parameters)

Calculate the skin friction coefficient (`Cf`) at a given boundary layer state and location.

# Arguments
- `state::Vector{Float64}`: Current boundary layer state vector `[d2, H1]`.
- `step::Float64`: The current position or step along the surface.
- `parameters::NamedTuple`: Parameters containing edge velocity, density, viscosity functions and tolerances:
  - `edge_velocity::Function`: Function returning edge velocity at given step.
  - `edge_density::Function`: Function returning edge density at given step.
  - `edge_viscosity::Function`: Function returning edge viscosity at given step.
  - `H1_eps::Float64`: Small epsilon used in limiting `H1`.
  - `H_eps::Float64`: Small epsilon used in calculating `H`.

# Returns
- `cf::Float64`: Calculated skin friction coefficient at the given state and position.
"""
function update_Cf_head(state, step, parameters)

    # Unpack State
    d2, H1 = state

    # Unpack parameters
    (; edge_velocity, edge_density, edge_viscosity, H1_eps, H_eps) = parameters

    Red2 = calculate_Re(edge_density(step), edge_velocity(step), d2, edge_viscosity(step))

    H = calculate_H(limH1(H1; eps=H1_eps); eps=H_eps)

    cf = calculate_cf(H, Red2)

    return cf
end

"""
    solve_head_boundary_layer!(::DiffEq, ode, initial_states, steps, parameters; verbose=false)

Integrate the turbulent boundary layer using a Runge-Kutta method.

# Arguments:
- `f::function_handle` : Governing residual equations to integrate
- `ode::function_handle` : ODE method to use (one of the DifferentialEquations.jl options)
- `initial_states::Float` : initial states
- `steps::Vector{Float}` : steps for integration
- `parameters::NamedTuple` : boundary layer solve options and other parameters

# Returns
- `usep::Vector{Float64}`: State vector at the separation point or last step
- `Hsep::Float64`: Shape factor at separation or last maximum
- `s_sep::Float64`: Surface coordinate corresponding to separation or last maximum
- `usol::Matrix{Float64}`: Solution states over all integration steps
- `stepsol::Vector{Float64}`: Integration step positions from the solver
"""
function solve_head_boundary_layer!(
    ::DiffEq, ode, initial_states, steps, parameters; verbose=false
)
    f = boundary_layer_residual_head!

    prob = ODEProblem(f, initial_states[1], [steps[1], steps[end]], parameters)

    if parameters.terminate
        # set up separation termination conditions
        function condition(u, t, integrator)
            return calculate_H(limH1(u[2])) - parameters.separation_criteria
            return update_Cf_head(u, t, parameters) + 1e-2
        end
        function affect!(integrator)
            return terminate!(integrator)
        end
        cb = ContinuousCallback(condition, affect!)

    else
        # TODO: figure out what default continuous callback is and use that here
        cb = nothing
    end

    sol = diffeq_solve(
        prob,
        ode();
        callback=cb,
        abstol=eps(),
        dtmax=1e-4, # TODO: put this in options and pass into parameters
        # alg_hints=[:stiff],
        # dt=parameters.first_step_size,
        verbose=verbose,
    )

    usol = reduce(hcat, (sol.u))
    stepsol = sol.t

    Hsol = calculate_H.(limH1.(usol[2, :]))
    Hsep = min(parameters.separation_criteria, maximum(Hsol))

    if Hsep == parameters.separation_criteria
        # println("separated")

        #spline H and steps from solution, get step at H=3
        sepwrap(x) = parameters.separation_criteria - smooth_akima(stepsol, Hsol, x)
        hid = findfirst(x -> x > parameters.separation_criteria, Hsol)
        s_sep = Roots.find_zero(sepwrap, [stepsol[hid - 1]; stepsol[hid]])
        # spline states and steps, get states at step for H=3
        usep = [smooth_akima(stepsol, u, s_sep) for u in eachrow(usol)]
    else
        # println("not separated")
        if parameters.return_last_max_shape_factor
            # println("return last max value")
            usep, Hsep, s_sep = find_last_max_H(usol, stepsol)
        else
            # println("return last value")
            usep = usol[:, end]
            Hsep = calculate_H(limH1(usep[2]))
            s_sep = stepsol[end]
        end
    end

    # return states at separate, and separation shape factor, and surface length at separation
    return usep, Hsep, s_sep, usol, stepsol
end
