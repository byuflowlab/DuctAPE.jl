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
move to utils.jl
"""
function split_bodies(vec, panels; duct=true)
    # get type of vector for consistent outputs
    TF = eltype(vec)

    #check if duct is used
    if !duct
        #hub only
        return TF[], TF[], vec, TF[], TF[], panels.panel_center[:, 1]
    else
        # get duct leading edge index. assumes duct comes first in vector
        _, leidx = findmin(panels[1].panel_center[:, 1])
        ndpan = length(panels[1].panel_center[:, 1])

        if length(panels) > 1
            #duct and hub
            return vec[1:leidx],
            vec[(leidx + 1):ndpan],
            vec[(ndpan + 1):end],
            panels[1].panel_center[1:leidx, 1],
            panels[1].panel_center[(leidx + 1):ndpan, 1],
            panels[2].panel_center[:, 1]
        else
            #duct only
            return vec[1:leidx],
            vec[(leidx + 1):ndpan],
            TF[],
            panels[1].panel_center[1:leidx, 1],
            panels[1].panel_center[(leidx + 1):ndpan, 1],
            TF[]
        end
    end

    # shouldn't get to this point...
    return nothing
end

# dot product
dot(A, B) = sum(a * b for (a, b) in zip(A, B))
# norm of vector
norm(A) = sqrt(mapreduce(x -> x^2, +, A))
