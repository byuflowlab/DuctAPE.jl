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
    k1, _ = f(y, s, parameters)
    parameters.verbose && println("  2nd call:")
    k2, aux = f(y .+ (ds / 2) .* k1, s + (ds / 2), parameters)
    unext = @. y + k2 * ds
    return unext, aux
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
    k1, aux1 = f(y, s, parameters)
    k2, aux2 = f(y .+ (ds / 2) .* k1, s + (ds / 2), parameters)
    k3, aux3 = f(y .+ (ds / 2) .* k2, s + (ds / 2), parameters)
    k4, aux4 = f(y .+ ds .* k3, s + ds, parameters)
    uaux = @. y + (k1 + k2 * 2 + k3 * 2 + k4) * ds / 6
    aux = [(aux1[i] + aux2[i] * 2 + aux3[i] * 2 + aux4[i]) / 6 for i in length(aux1)]
    return unext, aux
end

