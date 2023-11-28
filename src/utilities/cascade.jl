#=
Types and functions related to cascade data not included in CCBlade

Formatting copied directly from CCBlade source code.
=#

using Printf: @printf

"""
    parsefile(filename, radians, solidity)

Cascade version of `parsefile` function from CCBlade. Assumes stagger is given before reynolds and Mach number, and solidity is given after
"""
function parsecascadefile(filename, radians)
    inflow = Float64[]
    cl = Float64[]
    cd = Float64[]
    info = ""
    stagger = 1.0
    Re = 1.0
    Mach = 1.0
    solidity = 1.0

    open(filename) do f

        # skip header
        info = readline(f)
        stagger = parse(Float64, readline(f))
        Re = parse(Float64, readline(f))
        Mach = parse(Float64, readline(f))
        solidity = parse(Float64, readline(f))

        for line in eachline(f)
            parts = split(line)
            push!(inflow, parse(Float64, parts[1]))
            push!(cl, parse(Float64, parts[2]))
            push!(cd, parse(Float64, parts[3]))
        end
    end

    if !radians
        inflow *= pi / 180.0
    end

    return info, stagger, Re, Mach, solidity, inflow, cl, cd
end

"""
    writecascadefile(filename, info, Re, Mach, stagger, inflow, cl, cd, radians)

Cascade version of `writecascadefile` function from CCBlade. Writes solidity after Mach number
"""
function writecascadefile(
    filename, info, inflow, Re, stagger, solidity, Mach, cl, cd, radians
)
    open(filename, "w") do f
        @printf(f, "%s\n", info)
        @printf(f, "%.17g\n", stagger)
        @printf(f, "%.17g\n", Re)
        @printf(f, "%.17g\n", Mach)
        @printf(f, "%.17g\n", solidity)

        factor = 1.0
        if !radians
            factor = 180.0 / pi
        end

        for i in 1:length(inflow)
            @printf(f, "%.17g\t%.17g\t%.17g\n", inflow[i] * factor, cl[i], cd[i])
        end
    end

    return nothing
end

######################################################################
#                                                                    #
#                        CASCADE TYPES AND FUNCIONS                  #
#                                                                    #
######################################################################

##---------------------------------#
##        Stagger + Inflow         #
##---------------------------------#

#"""
#    StaggerInflowCAS(stagger, inflow, cl, cd, info, Mach)
#    StaggerInflowCAS(stagger, inflow, cl, cd, info)
#    StaggerInflowCAS(stagger, inflow, cl, cd)
#    read_StaggerInflowCAS(filenames::Vector{String}; radians=true)

#Data is fit with a recursive Akima spline.

#**Arguments**
#- `inflow::Vector{Float64}`: inflow angles
#- `stagger::Vector{Float64}`: stagger angle
#- `cl::Matrix{Float64}`: lift coefficients where cl[i, j] corresponds to inflow[i], stagger[j]
#- `cd::Matrix{Float64}`: drag coefficients where cd[i, j] corresponds to inflow[i], stagger[j]
#- `info::String`: a description of this airfoil data (just informational)
#- `Re::Float64`: Reynolds number data was taken at (just informational)
#- `Mach::Float64`: Mach number data was taken at (just informational)
#- `solidity::Float64`: solidity data was taken at (just informational)

#or

#filenames with one file per Reynolds number.

#**Arguments**
#- `filenames::Vector{String}`: name/path of files to read in, each at a different stagger angle in ascending order
#- `radians::Bool`: true if angle of attack in file is given in radians
#"""
#struct StaggerInflowCAS{TF,TS} <: DTCascade
#    stagger::Vector{TF}
#    inflow::Vector{TF}
#    cl::Matrix{TF}
#    cd::Matrix{TF}
#    info::TS # not used except for info in file
#    Re::TF # not used except for info in file
#    Mach::TF # not used except for info in file
#    solidity::TF # not used except for info in file
#end

#function StaggerInflowCAS(stagger, inflow, cl, cd)
#    return StaggerInflowCAS(
#        stagger, inflow, cl, cd, "DuctAPE written cascade", 0.0, 0.0, 0.0
#    )
#end

#function StaggerInflowCAS(filenames::Vector{String}; radians=true)

#    # read in first file
#    info, stagger1, Re_info, Mach_info, solidity_info, inflow, cl1, cd1 = parsecascadefile(
#        filenames[1], radians
#    )  # assumes common inflow across files, also common info etc.
#    ninflow = length(inflow)
#    nstagger = length(filenames)

#    cl = Array{Float64}(undef, ninflow, nstagger)
#    cd = Array{Float64}(undef, ninflow, nstagger)
#    stagger = Array{Float64}(undef, nstagger)

#    # iterate over files
#    for i in 1:nstagger
#        info, staggeri, _, _, _, _, cli, cdi = parsecascadefile(filenames[i], radians)
#        cl[:, i] = cli
#        cd[:, i] = cdi
#        stagger[i] = staggeri
#    end

#    return StaggerInflowCAS(
#        stagger, inflow, cl, cd, info, Re_info, Mach_info, solidity_info
#    )
#end

#function caseval(
#    cas::StaggerInflowCAS, stagger, inflow, Re=nothing, Mach=nothing, solidity=nothing
#)
#    cl = FLOWMath.interp2d(
#        FLOWMath.akima, cas.inflow, cas.stagger, cas.cl, [inflow], [stagger]
#    )[1]
#    cd = FLOWMath.interp2d(
#        FLOWMath.akima, cas.inflow, cas.stagger, cas.cd, [inflow], [stagger]
#    )[1]

#    return cl, cd
#end

#function writecascadefile(filenames, cas::StaggerInflowCAS; radians=true)
#    for i in 1:length(cas.stagger)
#        writecascadefile(
#            filenames[i],
#            cas.info,
#            cas.stagger[i],
#            cas.inflow,
#            cas.Re,
#            cas.Mach,
#            cas.solidity,
#            cas.cl[:, i],
#            cas.cd[:, i],
#            radians,
#        )
#    end
#    return nothing
#end

##---------------------------------#
##      Stagger + Inflow + Re      #
##---------------------------------#

#"""
#    StaggerInflowReCAS(stagger, inflow, Re, cl, cd, info)
#    StaggerInflowReCAS(stagger, inflow, Re, cl, cd)
#    StaggerInflowReCAS(filenames::Matrix{String}; radians=true)

#Data is fit with a recursive Akima spline.

#**Arguments**
#- `stagger::Vector{Float64}`: stagger angles
#- `alpha::Vector{Float64}`: inflow angles
#- `Re::Vector{Float64}`: Reynolds numbers
#- `cl::Array{Float64}`: lift coefficients where cl[i, j, k] corresponds to inflow[i], stagger[j], Re[k]
#- `cd::Array{Float64}`: drag coefficients where cd[i, j, k] corresponds to inflow[i], stagger[j], Re[k]
#- `info::String`: a description of this airfoil data (just informational)

#or files with one per Re/stagger combination

#**Arguments**
#- `filenames::Matrix{String}`: name/path of files to read in.  filenames[i, j] corresponds to stagger[i], Re[j] with stagger angle and Reynolds number in ascending order.
#- `radians::Bool`: true if angle of attack in file is given in radians
#"""
#struct StaggerInflowReCAS{TF,TS} <: DTCascade
#    stagger::Vector{TF}
#    inflow::Vector{TF}
#    Re::Vector{TF}
#    cl::Array{TF}
#    cd::Array{TF}
#    info::TS
#    Mach::TF
#    solidity::TF
#end

#function StaggerInflowReCAS(stagger, inflow, Re, cl, cd)
#    return StaggerInflowReCAS(
#        stagger, inflow, Re, cl, cd, "DuctAPE written cascade", 0.0, 0.0
#    )
#end

#function StaggerInflowReCAS(stagger, inflow, Re, cl, cd, info)
#    return StaggerInflowReCAS(stagger, inflow, Re, cl, cd, info, 0.0, 0.0)
#end

#function StaggerInflowReCAS(filenames::Matrix{String}; radians=true)

#    # - Get Dimensions - #
#    info, stagger1, Re1, Mach_info, solidity_info, inflow, cl1, cd1 = parsecascadefile(
#        filenames[1], radians
#    )  # assumes common inflow across files, also common info etc.

#    ninflow = length(inflow)
#    nstagger, nRe = size(filenames)

#    cl = Array{Float64}(undef, ninflow, nRe, nstagger)
#    cd = Array{Float64}(undef, ninflow, nRe, nstagger)
#    stagger = Array{Float64}(undef, nstagger)
#    Re = Array{Float64}(undef, nRe)

#    # iterate over files
#    for j in 1:nRe
#        for i in 1:nstagger
#            _, Rej, _, staggeri, _, clij, cdij = parsefile(filenames[i, j], radians)
#            cl[:, i, j] = clij
#            cd[:, i, j] = cdij
#            Re[j] = Rej
#            stagger[i] = staggeri  # NOTE: probably should add check to prevent user error here.
#        end
#    end

#    return StaggerInflowReCAS(stagger, inflow, Re, cl, cd, info)
#end

#function caseval(
#    cas::StaggerInflowReCAS, stagger, inflow, Re, Mach=nothing, solidity=nothing
#)
#    cl = FLOWMath.interp3d(
#        FLOWMath.akima, cas.stagger, cas.inflow, cas.Re, cas.cl, [stagger], [inflow], [Re]
#    )[1]
#    cd = FLOWMath.interp3d(
#        FLOWMath.akima, cas.stagger, cas.inflow, cas.Re, cas.cd, [stagger], [inflow], [Re]
#    )[1]

#    return cl, cd
#end

#function writecascadefile(filenames, cas::StaggerInflowReCAS; radians=true)
#    nre = length(cas.Re)
#    ns = length(cas.stagger)

#    for j in 1:nre
#        for i in 1:ns
#            writecascadefile(
#                filenames[i, j],
#                cas.info,
#                cas.stagger[i],
#                cas.inflow,
#                cas.Re[j],
#                cas.cl[:, i, j],
#                cas.cd[:, i, j],
#                cas.Mach,
#                cas.solidity,
#                radians,
#            )
#        end
#    end

#    return nothing
#end

##---------------------------------#
##   Stagger + Inflow + Re + Mach  #
##---------------------------------#
##TODO: need tests for these

#"""
#    StaggerInflowReMachCAS(stagger, inflow, Re, Mach, cl, cd, info, solidity)
#    StaggerInflowReMachCAS(stagger, inflow, Re, Mach, cl, cd)
#    StaggerInflowReMachCAS(filenames::Matrix{String}; radians=true)

#Data is fit with a recursive Akima spline.

#**Arguments**
#- `stagger::Vector{Float64}`: stagger angles
#- `inflow::Vector{Float64}`: inflow angles
#- `Re::Vector{Float64}`: Reynolds numbers
#- `Mach::Vector{Float64}`: Mach numbers
#- `cl::Array{Float64}`: lift coefficients where cl[i, j, k] corresponds to stagger[i], Re[j], Mach[k]
#- `cd::Array{Float64}`: drag coefficients where cd[i, j, k] corresponds to stagger[i], Re[j], Mach[k]
#- `info::String`: a description of this airfoil data (just informational)

#or files with one per Re/Mach combination

#**Arguments**
#- `filenames::Matrix{String}`: name/path of files to read in.  filenames[i, j, k] corresponds to stagger[i] Re[j] Mach[k] with each in ascending order.
#- `radians::Bool`: true if angle of attack in file is given in radians
#"""
#struct StaggerInflowReMachCAS{TF,TS} <: DTCascade
#    stagger::Vector{TF}
#    inflow::Vector{TF}
#    Re::Vector{TF}
#    Mach::Vector{TF}
#    cl::Array{TF}
#    cd::Array{TF}
#    info::TS
#    solidity::TF
#end

#function StaggerInflowReMachCAS(stagger, inflow, Re, Mach, cl, cd)
#    return StaggerInflowReMachCAS(
#        stagger, inflow, Re, Mach, cl, cd, "DuctAPE written cascade", 0.0
#    )
#end

#function StaggerInflowReMachCAS(stagger, inflow, Re, Mach, cl, cd, info)
#    return StaggerInflowReMachCAS(stagger, inflow, Re, Mach, cl, cd, info, 0.0)
#end

#function StaggerInflowReMachCAS(filenames::AbstractArray{String}; radians=true)
#    info, _, _, _, solidity, inflow, _, _ = parsecascadefile(filenames[1, 1, 1], radians)  # assumes common inflow and info across files
#    ninflow = length(inflow)
#    nstagger, nRe, nMach = size(filenames)

#    cl = Array{Float64}(undef, ninflow, nstagger, nRe, nMach)
#    cdrag = Array{Float64}(undef, ninflow, nstagger, nRe, nMach)
#    Re = Array{Float64}(undef, nRe)
#    Mach = Array{Float64}(undef, nMach)
#    stagger = Array{Float64}(undef, nstagger)

#    for k in 1:nMach
#        for j in 1:nRe
#            for i in 1:nstagger
#                _, staggeri, Rej, Machk, _, _, clijk, cdijk = parsecascadefile(
#                    filenames[i, j, k], radians
#                )
#                cl[:, i, j, k] = clijk
#                cdrag[:, i, j, k] = cdijk
#                stagger[i] = staggeri
#                Re[j] = Rej
#                Mach[k] = Machk
#            end
#        end
#    end

#    return StaggerInflowReMachCAS(stagger, inflow, Re, Mach, cl, cdrag, info, solidity)
#end

#function caseval(cas::StaggerInflowReMachCAS, stagger, inflow, Re, Mach, solidity=nothing)
#    cl = FLOWMath.interp4d(
#        FLOWMath.akima,
#        cas.stagger,
#        cas.inflow,
#        cas.Re,
#        cas.Mach,
#        cas.cl,
#        [stagger],
#        [inflow],
#        [Re],
#        [Mach],
#    )[1]
#    cd = FLOWMath.interp4d(
#        FLOWMath.akima,
#        cas.stagger,
#        cas.inflow,
#        cas.Re,
#        cas.Mach,
#        cas.cd,
#        [stagger],
#        [inflow],
#        [Re],
#        [Mach],
#    )[1]

#    return cl, cd
#end

#function writecascadefile(filenames, cas::StaggerInflowReMachCAS; radians=true)
#    for (k, Ma) in enumerate(cas.Mach)
#        for (j, Re) in enumerate(cas.Re)
#            for (i, stag) in enumerate(cas.stagger)
#                writecascadefile(
#                    filenames[i, j, k, ell],
#                    cas.info,
#                    stag,
#                    Re,
#                    Ma,
#                    cas.inflow,
#                    cas.cl[:, i, j, k],
#                    cas.cd[:, i, j, k],
#                    radians,
#                )
#            end
#        end
#    end

#    return nothing
#end

#-----------------------------------------#
# Inflow + Re + Stagger + Solidity + Mach #
#-----------------------------------------#
#TODO: need tests for these

"""
    InReStSoMaCAS(inflow, Re, stagger, solidity, Mach, cl, cd, info)
    InReStSoMaCAS(inflow, Re, stagger, solidity, Mach, cl, cd)
    InReStSoMaCAS(filenames::Matrix{String}; radians=true)

Data is fit recursively with Akima splines.

**Arguments**
- `inflow::Vector{Float64}`: inflow angles
- `Re::Vector{Float64}`: Reynolds numbers
- `stagger::Vector{Float64}`: stagger angles
- `solidity::Vector{Float64}`: local solidity
- `Mach::Vector{Float64}`: Mach numbers
- `cl::Array{Float64}`: lift coefficients where cl[i, j, k, ell] corresponds to stagger[i], Re[j], Mach[k], solidity[ell]
- `cd::Array{Float64}`: drag coefficients where cd[i, j, k, ell] corresponds to stagger[i], Re[j], Mach[k], solidity[ell]
- `info::String`: a description of this airfoil data (just informational)

or files with one per Re/Stagger/Solidty/Mach combination

**Arguments**
- `filenames::Matrix{String}`: name/path of files to read in.  filenames[i, j, k, ell] corresponds to Re[i] Stagger[j] Stagger[k] and Solidity[k] with each in ascending order.
- `radians::Bool`: true if angle of attack in file is given in radians
"""
struct InReStSoMaCAS{TF,TS} <: DTCascade
    inflow::Vector{TF}
    Re::Vector{TF}
    stagger::Vector{TF}
    solidity::Vector{TF}
    Mach::Vector{TF}
    cl::Array{TF}
    cd::Array{TF}
    info::TS
end

function InReStSoMaCAS(inflow, Re, stagger, solidity, Mach, cl, cd)
    return InReStSoMaCAS(
        inflow, Re, stagger, solidity, Mach, cl, cd, "DuctAPE written cascade"
    )
end

function InReStSoMaCAS(filenames::AbstractArray{String}; radians=true)
    info, _, _, _, _, inflow, _, _ = parsecascadefile(filenames[1, 1, 1, 1], radians)  # assumes common inflow and info across files
    ninflow = length(inflow)
    nRe, nstagger, nsolidity, nMach = size(filenames)

    cl = Array{Float64}(undef, ninflow, nRe, nstagger, nsolidity, nMach)
    cd = Array{Float64}(undef, ninflow, nRe, nstagger, nsolidity, nMach)
    Re = Array{Float64}(undef, nRe)
    stagger = Array{Float64}(undef, nstagger)
    solidity = Array{Float64}(undef, nsolidity)
    Mach = Array{Float64}(undef, nMach)

    for ell in 1:nMach
        for k in 1:nsolidity
            for j in 1:nstagger
                for i in 1:nRe
                    _, Rei, staggerj, solidityk, Machell, _, clijkl, cdijkl = parsecascadefile(
                        filenames[i, j, k, ell], radians
                    )
                    cl[:, i, j, k, ell] = clijkl
                    cd[:, i, j, k, ell] = cdijkl
                    Re[i] = Rei
                    stagger[j] = staggerj
                    solidity[k] = solidityk
                    Mach[ell] = Machell
                end
            end
        end
    end

    return InReStSoMaCAS(inflow, Re, stagger, solidity, Mach, cl, cd, info)
end

function caseval(cas::InReStSoMaCAS, inflow, Re, stagger, solidity, Mach)
    cl = interp5d(
        FLOWMath.akima,
        cas.inflow,
        cas.Re,
        cas.stagger,
        cas.solidity,
        cas.Mach,
        cas.cl,
        [inflow],
        [Re],
        [stagger],
        [solidity],
        [Mach],
    )[1]

    cd = interp5d(
        FLOWMath.akima,
        cas.inflow,
        cas.Re,
        cas.stagger,
        cas.solidity,
        cas.Mach,
        cas.cd,
        [inflow],
        [Re],
        [stagger],
        [solidity],
        [Mach],
    )[1]

    return cl, cd
end

function writecascadefile(filenames, cas::InReStSoMaCAS; radians=true)
    for (ell, Ma) in enumerate(cas.Ma)
        for (k, sol) in enumerate(cas.solidity)
            for (j, stag) in enumerate(cas.stagger)
                for (i, Re) in enumerate(cas.Re)
                    writecascadefile(
                        filenames[i, j, k, ell],
                        cas.info,
                        Re,
                        stag,
                        sol,
                        Ma,
                        cas.inflow,
                        cas.cl[:, i, j, k, ell],
                        cas.cd[:, i, j, k, ell],
                        radians,
                    )
                end
            end
        end
    end

    return nothing
end

# """
#     get_clcd(af, stagger, inflow, reynolds, mach, )

# Return lift and drag coefficients based on airfoil object type and flow conditions.

# **Arguments:**
#  - `cas::AFType` : Airfoil object either of a CCBlade airfoil type or custom type (custom if solidity information is included).
#  - `stagger::Vector{Float64}`: stagger angles
#  - `inflow::Vector{Float64}`: inflow angles
#  - `reynolds::Float` : Reynolds number
#  - `mach::Float` : Mach number

# **Returns:**
#  - `cl::Float` : section lift coefficient
#  - `cd::Float` : section drag coefficient
# """
#function get_clcd(af, stagger, inflow, reynolds, mach, solidity)

#    # if solidity undefined, use CCBlade functions
#    if stagger == nothing
#        return ccb.afeval(af, inflow, reynolds, mach)
#    else
#        #if solidity is defined, need to use custom functions
#        return caseval(af, stagger, inflow, reynolds, mach, solidity)
#    end
#end

"""
     interp5d(interp1d, x1data, x2data, x3data, x4data, fdata, x1pt, x2pt, x3pt, x4pt)

Same as FLOWMath.interp4d, ex1cept in five dimensions.
"""
function interp5d(
    interp1d, x1data, x2data, x3data, x4data, x5data, fdata, x1pt, x2pt, x3pt, x4pt, x5pt
)
    nd = length(x5data)
    nx1pt = length(x1pt)
    nx2pt = length(x2pt)
    nx3pt = length(x3pt)
    nx4pt = length(x4pt)
    nx5pt = length(x5pt)

    R = promote_type(eltype(x1pt), eltype(x2pt), eltype(x3pt), eltype(x4pt), eltype(x5pt))
    x5interp = Array{R}(undef, nd, nx1pt, nx2pt, nx3pt, nx4pt)
    output = Array{R}(undef, nx1pt, nx2pt, nx3pt, nx4pt, nx5pt)

    for i in 1:nd
        x5interp[i, :, :, :] .= interp4d(
            interp1d,
            x1data,
            x2data,
            x3data,
            x4data,
            fdata[:, :, :, :, i],
            x1pt,
            x2pt,
            x3pt,
            x4pt,
        )
    end
    for ell in 1:nx4pt
        for k in 1:nx3pt
            for j in 1:nx2pt
                for i in 1:nx1pt
                    output[i, j, k, ell, :] .= interp1d(
                        x5data, x5interp[:, i, j, k, ell], x5pt
                    )
                end
            end
        end
    end

    return output
end
