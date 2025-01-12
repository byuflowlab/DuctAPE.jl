#---------------------------------#
#           Runge-Kutta           #
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
    k1 = f(y, parameters, s)
    parameters.verbose && println("  2nd call:")
    k2 = f(y .+ (ds / 2) .* k1, parameters, s + (ds / 2))
    unext = @. y + k2 * ds
    return unext
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
    k1 = f(y, parameters, s)
    k2 = f(y .+ (ds / 2) .* k1, parameters, s + (ds / 2))
    k3 = f(y .+ (ds / 2) .* k2, parameters, s + (ds / 2))
    k4 = f(y .+ ds .* k3, parameters, s + ds)
    unext = @. y + (k1 + k2 * 2 + k3 * 2 + k4) * ds / 6
    return unext
end
