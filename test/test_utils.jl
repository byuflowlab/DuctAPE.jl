"""
check monotinicity of array
"""
function ismonotonic(A::AbstractArray, column::Int, cmp=<)
    current = A[begin, column] # begin instead of 1
    for i in axes(A, 1)[2:end] # skip the first element
        newval = A[i, column] # don't use copy here
        cmp(newval, current) && return false
        current = newval
    end
    return true
end

#---------------------------------#
#    Compare Structs and Tuples   #
#---------------------------------#
function compare_structs(obj1, obj2)
    getfieldnames(obj) = fieldnames(typeof(obj))

    names = getfieldnames(obj1)
    flags = Vector{Bool}(undef, length(names))

    if !isa(obj1, typeof(obj2))
        @warn("compare_fieldnames(): the objects are not of the same type.")
        flags .= false
        return flags
    end

    for i in eachindex(names)
        i1 = getproperty(obj1, names[i])
        i2 = getproperty(obj2, names[i])
        if isapprox(i1, i2)
            flags[i] = true
        else
            flags[i] = false
        end
    end
    return flags
end

function compare_namedtuples(tuple1, tuple2)
    flags = Bool[]
    for field in keys(tuple1)
        try
            append!(
                flags,
                compare_structs(getproperty(tuple1, field), getproperty(tuple2, field)),
            )
        catch
            if typeof(getproperty(tuple1, field)) <: NamedTuple
                append!(
                    flags,
                    compare_namedtuples(
                        getproperty(tuple1, field), getproperty(tuple2, field)
                    ),
                )
            elseif typeof(getproperty(tuple1, field)) <: AbstractArray
                for i in 1:length(getproperty(tuple1, field))
                    try
                        append!(
                            flags,
                            compare_structs(
                                getproperty(tuple1, field)[i], getproperty(tuple2, field)[i]
                            ),
                        )
                    catch
                        if typeof(getproperty(tuple1, field)[i]) <: NamedTuple
                            append!(
                                flags,
                                compare_namedtuples(
                                    getproperty(tuple1, field)[i],
                                    getproperty(tuple2, field)[i],
                                ),
                            )
                        else
                            append!(
                                flags,
                                isapprox(
                                    getproperty(tuple1, field)[i],
                                    getproperty(tuple2, field)[i],
                                ),
                            )
                        end #if named tuple
                    end
                end #for entries in vector
            else
                append!(
                    flags, isapprox(getproperty(tuple1, field), getproperty(tuple2, field))
                )
            end #if/else
        end
    end #for fields
    return flags
end

