"""
allows us to get rotor parameters easier
"""
function Base.getproperty(obj::AbstractVector{<:NamedTuple}, sym::Symbol)
    return getfield.(obj, sym)
end
