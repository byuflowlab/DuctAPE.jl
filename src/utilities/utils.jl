"""
allows us to get rotor parameters easier
"""
function Base.getproperty(obj::AbstractVector{<:NamedTuple}, sym::Symbol)
    return getfield.(obj, sym)
end

# - Function for adding in xlocations - #
# from https://stackoverflow.com/questions/25678112/insert-item-into-a-sorted-list-with-julia-with-and-without-duplicates
function insert_and_dedup!(v, x)
    for i in 1:length(x)
        # find ranges and replace with discrete values (thus deleting duplicates if present)
        v = (splice!(v, searchsorted(v, x[i]), x[i]); v)
    end
end

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

# dot product
dot(A, B) = sum(a * b for (a, b) in zip(A, B))
# norm of vector
norm(A) = sqrt(mapreduce(x -> x^2, +, A))
# 2D "cross product" magnitude
cross2mag(A, B) = A[1] * B[2] - A[2] * B[1]

"""
    linear_transform(range1, range2, values)

Linear transfrom of values from range (source_range[1], raend) to (target_range[1], target_range[end])

# Arguments:
- `source_range::Vector{Float{` : range values come from
- `target_range::Vector{Float}` : range onto which we are transforming
- `source_values::Array{Float}` : array of source_values to transform

# Returns:
 - `target_values::Array{Float}` : array of transformed source_values onto target range
"""
function linear_transform(source_range, target_range, source_values)
    return target_range[1] .+
           (target_range[end] - target_range[1]) .* (source_values .- source_range[1]) /
           (source_range[end] - source_range[1])
end

"""
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
length from size
move to utilities
"""
function lfs(shape)
    if length(shape) == 1
        return shape[1]
    else
        return *(shape...)
    end
end

"""
note: containers must be Arrays, structs of arrays, or tuples of arrays
move to utilities
"""
function reset_containers!(c)
    if typeof(c) <: AbstractArray
        #do nothing if it's a string
        (eltype(c) == String) || (c .= 0)
    else
        for p in propertynames(c)
            cp = getfield(c, p)
            if typeof(cp) <: AbstractArray
                #do nothing if it's a string
                (eltype(cp) == String) || (cp .= 0)
            else
                reset_containers!(cp)
            end
        end
    end

    return c
end
