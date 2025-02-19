# - Functions to check if the input is a scalar - #
import Base.BroadcastStyle
"""
    isscalar(x::T) where {T} = isscalar(T)
    isscalar(::Type{T}) where {T} = BroadcastStyle(T) isa Broadcast.DefaultArrayStyle{0}

Determines if the input is a scalar. Note that `Base.BroadcastStyle` is imported.
"""
isscalar(x::T) where {T} = isscalar(T)
isscalar(::Type{T}) where {T} = BroadcastStyle(T) isa Broadcast.DefaultArrayStyle{0}

"""
    printval(text, val)

Used for debugging; prints values of `val`, preceeded by `text` for floats, strings, ints and ForwardDiff duals without the user having to know if duals are being used or not.
"""
function printval(text, val)
    if eltype(val) == Float64
        println(text, val, " (Float)")
    elseif eltype(val) <: Int
        println(text, val, " (Int)")
    elseif eltype(val) <: String
        println(text, val, " (String)")
    else
        if typeof(val) <: AbstractMatrix || typeof(val) <: AbstractVector
            println(text, (p -> p.value).(val), " (ForDiff)")
        else
            println(text, val.value, " (ForDiff)")
        end
    end
    return nothing
end

"""
    printdebug(variable_name, variable, nspaces=4)

Formatted printing when you have lots of debugging to do and want things to line up nicely.
"""
function printdebug(variable_name, variable, nspaces=4; colw=16)
    if length(variable) == 1
        s = @sprintf "%*s %-*s %2.18f" nspaces " " colw variable_name ForwardDiff.value(
            variable[1]
        )
        println(rstrip(s, ['0', ' ']))
    else
        @printf "%*s %-*s\n" nspaces " " colw variable_name
        for v in variable
            s = @sprintf "%*s %-*s %2.18f" (nspaces) " " colw "" ForwardDiff.value(v)
            println(rstrip(s, ['0', ' ']))
        end
    end
    return nothing
end

"""
    dot(A, B) = sum(a * b for (a, b) in zip(A, B))

A faster dot product.
"""
dot(A, B) = sum(a * b for (a, b) in zip(A, B))

"""
    norm(A) = sqrt(mapreduce(x -> x^2, +, A))

A faster 2-norm.
"""
norm(A) = sqrt(mapreduce(x -> x^2, +, A))

"""
    cross2mag(A, B) = A[1] * B[2] - A[2] * B[1]

2D "cross product" magnitude
"""
cross2mag(A, B) = A[1] * B[2] - A[2] * B[1]

"""
    linear_transform(range1, range2, values)

Linear transfrom of values from range `(source_range[1], source_range[end])` to `(target_range[1], target_range[end])`

# Arguments
- `source_range::Vector{Float}` : range values come from (can also be a Tuple)
- `target_range::Vector{Float}` : range onto which we are transforming (can also be a Tuple)
- `source_values::Array{Float}` : array of source values to transform

# Returns
 - `target_values::Array{Float}` : array of transformed sourcevalues onto target range
"""
function linear_transform(source_range, target_range, source_values)
    return target_range[1] .+
           (target_range[end] - target_range[1]) .* (source_values .- source_range[1]) /
           (source_range[end] - source_range[1])
end

"""
    extract_primals!(Avalue, A::AbstractMatrix{T}) where {T}

Extracts primals of A and places them in Avalue.
"""
function extract_primals!(Avalue, A::AbstractMatrix{T}) where {T}
    if T <: ForwardDiff.Dual #|| T<:ReverseDiff.TrackedReal  # Automatic differentiation case

        # Extract primal values of A
        # value = T<:ForwardDiff.Dual ? ForwardDiff.value : ReverseDiff.value
        value = ForwardDiff.value
        map!(value, Avalue, A)

    else                                # Normal case
        # Deep copy A
        Avalue .= A
    end

    return Avalue
end

"""
    lfs(shape)

Determines length from shape (output of `size` function).
"""
function lfs(shape)
    if length(shape) == 1
        return shape[1]
    else
        return *(shape...)
    end
end

"""
    reset_containers!(containers; exception_keys=[])

Resets all fields (not incluing any contained in exception keys) of containers---which must be arrays, structs of arrays, or tuples of arrays---to zeros.
"""
function reset_containers!(c; exception_keys=[])
    if typeof(c) <: AbstractArray
        #do nothing if it's a string
        (eltype(c) == String) || (c .= 0)
    else
        for f in fieldnames(typeof(c))
            if !(f in exception_keys)
                cp = getfield(c, f)
                if typeof(cp) <: AbstractArray
                    if eltype(cp) <: Tuple
                        for i in 1:length(cp[1])
                            for j in eachindex(cp)
                                cp[j][i] .= 0.0
                            end
                        end
                    else
                        #do nothing if it's a string
                        (eltype(cp) == String) || (cp .= 0)
                    end
                else
                    reset_containers!(cp; exception_keys=exception_keys)
                end
            end
        end
    end

    return c
end

"""
    promote_ducted_rotor_type(ducted_rotor, operating_point)

Convenience function for promoting types based on any potential elements of the ducted_rotor object dependent on optimization design variables.

# Arguments
- `ducted_rotor::DuctedRotor` : the ducted_rotor input
- `operating_point::OperatingPoint` : the operating_point input

# Returns
- `TP::Type` : the promoted type
"""
function promote_ducted_rotor_type(d, o)
    return promote_type(
        eltype(d.duct_coordinates),
        eltype(d.centerbody_coordinates),
        eltype(o.Vinf),
        eltype(o.rhoinf),
        eltype(o.muinf),
        eltype(o.asound),
        eltype(o.Omega),
        eltype(d.rotor.B),
        eltype(d.rotor.rotorzloc),
        eltype(d.rotor.r),
        eltype(d.rotor.Rhub),
        eltype(d.rotor.Rtip),
        eltype(d.rotor.chords),
        eltype(d.rotor.twists),
    )
end

"""
    cache_dims!(total_length, l, s)

A function that returns a named tuple containing an index range and shape and increases `total_length` by `l`.

This function is used heavily in the cache allocation functions for setting up the dimension maps used to access the vectorized caches.

# Arguments
- `total_length::Vector{Int}` : single element vector containing the current total length of the eventual cache vector. Modified in place.
- `l::Int` : total length of the object in question
- `s::Int` : size of the object in question

# Returns
- `dims::NamedTuple` : A named tuple containing `index` and `shape` fields
"""
function cache_dims!(total_length, l, s)
    dims = (; index=(total_length[] + 1):(total_length[] + l), shape=s)
    total_length[] += l
    return dims
end

"""
    smooth_akima(x, y, xpt; delta=2.0 * eps(), eps=eps())

Wrapper for FLOWMath.akima with different optional argument values.
"""
function smooth_akima(x, y, xpt; delta=2.0 * eps(), eps=eps())
    return FLOWMath.akima(x, y, xpt, delta, eps)
end

"""
    smooth_Akima(x, y; delta=2.0 * eps(), eps=eps())

Wrapper for FLOWMath.Akima with different optional argument values.
"""
function smooth_Akima(x, y; delta=2.0 * eps(), eps=eps())
    return FLOWMath.Akima(x, y, delta, eps)
end

"""
"""
function smooth_abs(x; eps=eps())
    return FLOWMath.abs_smooth.(x, Ref(eps))
end
