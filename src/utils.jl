"""
allows us to get rotor parameters easier
"""
function Base.getproperty(obj::AbstractVector{<:NamedTuple}, sym::Symbol)
    return getfield.(obj, sym)
end

"""
    function cosine_spacing(N=180)

Calculate N cosine spaced points.

**Arguments:**
 - `N::Int` : Number of points

"""
function cosine_spacing(N=80)
    return scaled_cosine_spacing(N, 1.0, 0.0; mypi=pi)
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
"""
function split_bodies(vec, body_panels; duct=true, hub=true)
    # get type of vector for consistent outputs
    TF = eltype(vec)

    #check if duct is used
    if !duct
        #hub only
        return TF[], TF[], vec, TF[], TF[], body_panels.controlpoint[:, 1]
    else
        # split duct into inner and outer
        ndpan = body_panels.endpointidxs[1, 2]
        # get duct leading edge index. assumes duct comes first in vector
        _, leidx = findmin(body_panels.controlpoint[1:ndpan, 1])
        if !hub
            #duct only
            return vec[1:leidx,:],
            vec[(leidx + 1):ndpan,:],
            TF[],
            body_panels.controlpoint[1:leidx, 1],
            body_panels.controlpoint[(leidx + 1):ndpan, 1],
            TF[]
        else
            #duct and hub
            return vec[1:leidx,:],
            vec[(leidx + 1):ndpan,:],
            vec[(ndpan + 1):end,:],
            body_panels.controlpoint[1:leidx, 1],
            body_panels.controlpoint[(leidx + 1):ndpan, 1],
            body_panels.controlpoint[(ndpan + 1):end, 1]
        end
    end

    # shouldn't get to this point...
    return nothing
end

# dot product
dot(A, B) = sum(a * b for (a, b) in zip(A, B))
# norm of vector
norm(A) = sqrt(mapreduce(x -> x^2, +, A))
# 2D "cross product" magnitude
cross2mag(A,B) = A[1]*B[2] - A[2]*B[1]

"""
    repanel_airfoil(x,y;N)
    repanel_airfoil(xy;N)
Takes x and y coordinates of an airfoil  and uses a cosine spaced akima spline to fill in the gaps
**Arguments**
- `x::Vector{Float64}` : vector containing the x coordinates of the airfoil
- `y::Vector{Float64}` : vector containing the y components of the airfoil
- `xy::Array{Float64,2}` : Array of x and y coordinates with X in the first column and y in the 2nd
**Keyword Arguements**
- `N::Int` : Number of data points to be returned after repaneling. Will only return odd numbers, if N is even, N+1 points will be returned.
**Returns**
- `xreturn::Vector{Float64}` : Repaneled, cosine spaced x corrdinates of the airfoil
- `yreturn::Vector{Float64}` : y coordinates of the repaneled airfoil obtained using an akima spline
- `xyreturn::Array{Float64}` : If the coordinates were input as an array, this will be returned with x in the 1st column and y in the 2nd.
"""
function repanel_airfoil(x, y; N=160, normalize=true)
    @assert length(x) == length(y) "X and Y vectors must be the same length"

    if normalize
        scale = 1.0
        le = 0.0
    else
        scale = maximum(x) - minimum(x)
        le = minimum(x)
    end
    #First normalize the airfoil to between 0 and 1
    normalize_airfoil!(x, y)

    #let's figure out the cosine spacing.
    npoints = ceil(Int, N / 2)
    akimax = cosine_spacing(npoints)

    #now we split the top and bottom of the airfoil
    x1, x2, y1, y2 = split_upper_lower(x, y)

    #Now check and see which x and y need to be reversed
    #x has to be ascending (0-->1)

    if x1[1] > x1[end]
        x1 = reverse(x1)
        y1 = reverse(y1)
    end

    if x2[1] > x2[end]
        x2 = reverse(x2)
        y2 = reverse(y2)
    end

    #do the akima spline
    akimay1 = FLOWMath.akima(x1, y1, akimax)
    akimay2 = FLOWMath.akima(x2, y2, akimax)

    #figure out which spline is on top
    if maximum(akimay1) > maximum(akimay2) #then akimay1 is on top so I need to reverse akimay2
        yreturn = [reverse(akimay2); akimay1[2:end]]

    else #otherwise akimay2 is on top so I need to reverse akimay1
        yreturn = [reverse(akimay1); akimay2[2:end]]
    end

    xreturn = [reverse(akimax); akimax[2:end]]

    return xreturn, yreturn
end

function repanel_airfoil(coordinates; N=160, normalize=true)
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    xpane, ypane = repanel_airfoil(x, y; N=N, normalize=normalize)

    return [xpane ypane]
end

function repanel_revolution(coordinates; N=160, normalize=true)
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    xpane, ypane = repanel_revolution(x, y; N=N, normalize=normalize)

    return [xpane ypane]
end

function repanel_revolution(x, y; N=160, normalize=true)
    @assert length(x) == length(y) "X and Y vectors must be the same length"

    if normalize
        scale = 1.0
        le = 0.0
    else
        scale = maximum(x) - minimum(x)
        le = minimum(x)
    end
    #First normalize the airfoil to between 0 and 1
    normalize_airfoil!(x, y)

    #let's figure out the cosine spacing.
    npoints = ceil(Int, N)
    akimax = cosine_spacing(npoints)

    #Now check and see which x and y need to be reversed
    #x has to be ascending (0-->1)

    if x[1] > x[end]
        x = reverse(x)
        y = reverse(y)
    end

    #do the akima spline
    akimay = fm.akima(x, y, akimax)

    return akimax * scale .+ le, akimay * scale
end

"""
    split_upper_lower(x, z)

Split the upper and lower halves of the airfoil coordinates. Assumes odd number of coordinates (leading edge repeated).

**Arguments:**
 - `x::Array{Float}` : Array of x coordinates
 - `z::Array{Float}` : Array of z coordinates

**Returns:**
 - `xu::Array{Float}` : Array of upper half of x coordinates
 - `xl::Array{Float}` : Array of lower half of x coordinates
 - `zu::Array{Float}` : Array of upper half of z coordinates
 - `zl::Array{Float}` : Array of lower half of z coordinates

"""
function split_upper_lower(x, z)

    # get half length of geometry coordinates
    _, N = findmin(x)

    return x[1:N], x[N:end], z[1:N], z[N:end]
end

"""
    normalize_airfoil!(x, z)

Normalize airfoil to unit chord and shift leading edge to zero. Adjusts coordinates in place.

**Arguments:**
 - `x::Array{Float}` : Array of x coordinates
 - `z::Array{Float}` : Array of z coordinates

"""
function normalize_airfoil!(x, z)
    chord = maximum(x) - minimum(x) #get current chord length
    x .-= minimum(x) #shift to zero
    x ./= chord #normalize chord
    z ./= chord #scale z coordinates to match

    return nothing
end

"""
finds local maxima on given vector, x
"""
function findlocalmax(x::AbstractVector{TF}) where {TF}
    inds = Int[]
    if length(x) > 1
        if x[1] > x[2]
            push!(inds, 1)
        end
        for i in 2:(length(x) - 1)
            if x[i - 1] < x[i] > x[i + 1]
                push!(inds, i)
            end
        end
        if x[end] > x[end - 1]
            push!(inds, length(x))
        end
    end
    return inds
end

"""
finds local minima on given vector, x
"""
function findlocalmin(x::AbstractVector{TF}) where {TF}
    inds = Int[]
    if length(x) > 1
        if x[2] > x[1]
            push!(inds, 1)
        end
        for i in 2:(length(x) - 1)
            if x[i + 1] > x[i] < x[i - 1]
                push!(inds, i)
            end
        end
        if x[end - 1] > x[end]
            push!(inds, length(x))
        end
    end
    return inds
end

"""
    linear_transform(range1, range2, values)

Linear transfrom of values from range (source_range[1], raend) to (target_range[1], target_range[end])

**Arguments:**
- `source_range::Vector{Float{` : range values come from
- `target_range::Vector{Float}` : range onto which we are transforming
- `source_values::Array{Float}` : array of source_values to transform

**Returns:**
 - `target_values::Array{Float}` : array of transformed source_values onto target range
"""
function linear_transform(source_range, target_range, source_values)
    return target_range[1] .+
           (target_range[end] - target_range[1]) .* (source_values .- source_range[1]) /
           (source_range[end] - source_range[1])
end
