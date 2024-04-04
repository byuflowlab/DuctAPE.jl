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

#-----------------------------------------------------#
#    Compare Structs and Tuples and Arrays of such    #
#-----------------------------------------------------#

getfieldnames(obj) = fieldnames(typeof(obj))

function compare_structs(
    obj1, obj2, names; current=nothing, parent="", tab="", verbose=false
)
    if verbose
        if isnothing(current)
            println("Comparing Structs:")
        else
            if current == ""
                println(tab * "$(parent)")
            else
                println(tab * "$(parent).$(current)")
            end
        end
    end

    flags = Bool[]

    for name in names
        if typeof(getproperty(obj1, name)) <: NamedTuple
            append!(
                flags,
                compare_namedtuples(
                    getproperty(obj1, name),
                    getproperty(obj2, name);
                    current=name,
                    parent=parent,
                    tab=tab * "  ",
                    verbose=verbose,
                ),
            )
        elseif typeof(getproperty(obj1, name)) <: AbstractArray
            append!(
                flags,
                compare_arrays(
                    getproperty(obj1, name),
                    getproperty(obj2, name);
                    current=name,
                    parent=parent,
                    tab=tab * "  ",
                    verbose=verbose,
                ),
            )
        else
            try
                subnames = getfieldnames(getproperty(obj1, name))
                if isempty(subnames)
                    throw(subnames)
                end
                append!(
                    flags,
                    compare_structs(
                        getproperty(obj1, name),
                        getproperty(obj2, name),
                        subnames;
                        current="",
                        parent=parent * ".$(name)",
                        tab=tab * "  ",
                        verbose=verbose,
                    ),
                )
            catch
                if verbose
                    println(tab * "  $(parent).$(name)")
                end
                append!(flags, isapprox(getproperty(obj1, name), getproperty(obj2, name)))
            end #if/else
        end #catch
    end
    return flags
end

function compare_namedtuples(
    tuple1, tuple2; current=nothing, parent="", tab="", verbose=false
)
    if verbose
        if isnothing(current)
            println("Comparing NamedTuples:")
        else
            if current == ""
                println(tab * "$(parent)")
            else
                println(tab * "$(parent).$(current)")
            end
        end
    end

    flags = Bool[]
    for field in keys(tuple1)
        if typeof(getproperty(tuple1, field)) <: NamedTuple
            append!(
                flags,
                compare_namedtuples(
                    getproperty(tuple1, field),
                    getproperty(tuple2, field);
                    current="",
                    parent=parent * ".$(field)",
                    tab=tab * "  ",
                    verbose=verbose,
                ),
            )
        elseif typeof(getproperty(tuple1, field)) <: AbstractArray
            append!(
                flags,
                compare_arrays(
                    getproperty(tuple1, field),
                    getproperty(tuple2, field);
                    current="",
                    parent=parent * ".$(field)",
                    tab=tab * "  ",
                    verbose=verbose,
                ),
            )
        else
            try
                names = getfieldnames(getproperty(tuple1, field))
                if isempty(names)
                    throw(names)
                end
                append!(
                    flags,
                    compare_structs(
                        getproperty(tuple1, field),
                        getproperty(tuple2, field),
                        names;
                        current="",
                        parent=parent * ".$(field)",
                        tab=tab * "  ",
                        verbose=verbose,
                    ),
                )
            catch
                if verbose
                    println(tab * "  $(parent).$(field)")
                end
                append!(
                    flags, isapprox(getproperty(tuple1, field), getproperty(tuple2, field))
                )
            end #if/else
        end #catch
    end #for fields
    return flags
end

function compare_arrays(arr1, arr2; current=nothing, parent="", tab="", verbose=false)
    if verbose
        if isnothing(current)
            println("Comparing Arrays:")
        else
            if current == ""
                println(tab * "$(parent)")
            else
                println(tab * "$(parent).$(current)")
            end
        end
    end
    flags = Bool[]
    if eltype(arr1) <: NamedTuple
        for i in 1:length(arr1)
            append!(
                flags,
                compare_namedtuples(
                    arr1[i],
                    arr2[i];
                    current="",
                    parent=parent * "$(current)[$i]",
                    tab=tab * "  ",
                    verbose=verbose,
                ),
            )
        end
    elseif eltype(arr1) <: AbstractArray
        for i in 1:length(arr1)
            append!(
                flags,
                compare_arrays(
                    arr1[i],
                    arr2[i];
                    current="",
                    parent=parent * "$(current)[$i]",
                    tab=tab * "  ",
                    verbose=verbose,
                ),
            )
        end
    else
        try
            names = getfieldnames(arr1[1])
            if isempty(names)
                throw(names)
            end
            for i in 1:length(arr1)
                append!(
                    flags,
                    compare_structs(
                        arr1[i],
                        arr2[i],
                        names;
                        current="",
                        parent=parent * "$(current)[$i]",
                        tab=tab * "  ",
                        verbose=verbose,
                    ),
                )
            end
        catch
            append!(flags, isapprox(arr1, arr2; atol=1e-3))
        end #if named tuple
    end
    return flags
end
