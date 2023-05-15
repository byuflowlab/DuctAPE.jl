"""
allows us to get rotor parameters easier
"""
function Base.getproperty(obj::AbstractVector{<:NamedTuple}, sym::Symbol)
    return getfield.(obj, sym)
end

"""
cosine spacing, but also scales and transforms
"""
function scaled_cosine_spacing(N, scale, transform; mypi=pi)
    return transform .+ scale * [0.5 * (1 - cos(mypi * (i - 1) / (N - 1))) for i in 1:N]
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
    if eltype(val) != Float64
        println(text, val.value)
    else
        println(text, val)
    end
    return nothing
end
