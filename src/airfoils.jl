#=
Types and functions related to airfoil data not included in CCBlade, i.e., airfoil data defined with an additional parameter related to solidity.

Formatting copied directly from CCBlade source code.
=#

using Printf: @printf

"""
overload `parsefile` function from CCBlade
assumes solidity parameter after mach number
"""
function parsefile(filename, radians, solidity)
    alpha = Float64[]
    cl = Float64[]
    cd = Float64[]
    info = ""
    Re = 1.0
    Mach = 1.0
    solidity = 1.0

    open(filename) do f

        # skip header
        info = readline(f)
        Re = parse(Float64, readline(f))
        Mach = parse(Float64, readline(f))
        if solidity
            solidity = parse(Float64, readline(f))
        end

        for line in eachline(f)
            parts = split(line)
            push!(alpha, parse(Float64, parts[1]))
            push!(cl, parse(Float64, parts[2]))
            push!(cd, parse(Float64, parts[3]))
        end
    end

    if !radians
        alpha *= pi / 180.0
    end

    if solidity
        return info, Re, Mach, solidity, alpha, cl, cd
    else
        return info, Re, Mach, alpha, cl, cd
    end
end

"""
overload `writefile` function from CCBlade to include solidity parameter in file header
"""
function writefile(filename, info, Re, Mach, solidity, alpha, cl, cd, radians)
    open(filename, "w") do f
        @printf(f, "%s\n", info)
        @printf(f, "%.17g\n", Re)
        @printf(f, "%.17g\n", Mach)
        @printf(f, "%.17g\n", solidity)

        factor = 1.0
        if !radians
            factor = 180.0 / pi
        end

        for i in 1:length(alpha)
            @printf(f, "%.17g\t%.17g\t%.17g\n", alpha[i] * factor, cl[i], cd[i])
        end
    end

    return nothing
end

#=
Add CCBlade-like airfoil types and functions for solidity parameter
=#
"""
    AlphaSolidityAF(alpha, solidity, cl, cd, info, Mach)
    AlphaSolidityAF(alpha, solidity, cl, cd, info)
    AlphaSolidityAF(alpha, solidity, cl, cd)
    read_AlphaSolidityAF(filenames::Vector{String}; radians=true)

Airfoil data that varies with angle of attack and Reynolds number.
Data is fit with a recursive Akima spline.

**Arguments**
- `alpha::Vector{Float64}`: angles of attack
- `solidity::Vector{Float64}`: solidity parameter
- `cl::Matrix{Float64}`: lift coefficients where cl[i, j] corresponds to alpha[i], solidity[j]
- `cd::Matrix{Float64}`: drag coefficients where cd[i, j] corresponds to alpha[i], solidity[j]
- `info::String`: a description of this airfoil data (just informational)
- `Mach::Float64`: Mach number data was taken at (just informational)

or

filenames with one file per Reynolds number.

**Arguments**
- `filenames::Vector{String}`: name/path of files to read in, each at a different Reynolds number in ascending order
- `radians::Bool`: true if angle of attack in file is given in radians
"""
struct AlphaSolidityAF{TF,TS} <: AFType
    alpha::Vector{TF}
    solidity::Vector{TF}
    cl::Matrix{TF}
    cd::Matrix{TF}
    info::TS # not used except for info in file
    Mach::TF # not used except for info in file
end

function AlphaSolidityAF(alpha, solidity, cl, cd, info)
    return AlphaSolidityAF(alpha, solidity, cl, cd, info, 0.0)
end
function AlphaSolidityAF(alpha, solidity, cl, cd)
    return AlphaSolidityAF(alpha, solidity, cl, cd, "CCBlade generated airfoil", 0.0)
end

function AlphaSolidityAF(filenames::Vector{String}; radians=true)
    info, solidity1, Mach, alpha, cl1, cd1 = parsefile(filenames[1], radians)  # assumes common alpha across files, also common info and common Mach
    nalpha = length(alpha)
    ncond = length(filenames)

    cl = Array{Float64}(undef, nalpha, ncond)
    cd = Array{Float64}(undef, nalpha, ncond)
    solidity = Array{Float64}(undef, ncond)
    cl[:, 1] = cl1
    cd[:, 1] = cd1
    solidity[1] = solidity1

    # iterate over remaining files
    for i in 2:ncond
        _, _, _, solidityi, _, cli, cdi = parsefile(filenames[i], radians)
        cl[:, i] = cli
        cd[:, i] = cdi
        solidity[i] = solidityi
    end

    return AlphaSolidityAF(alpha, solidity, cl, cd, info, Mach)
end

function afeval(af::AlphaSolidityAF, alpha, Re, Mach, solidity)
    cl = FLOWMath.interp2d(
        FLOWMath.akima, af.alpha, af.solidity, af.cl, [alpha], [solidity]
    )[1]
    cd = FLOWMath.interp2d(
        FLOWMath.akima, af.alpha, af.solidity, af.cd, [alpha], [solidity]
    )[1]

    return cl, cd
end

function write_af(filenames, af::AlphaSolidityAF; radians=true)
    for i in 1:length(af.solidity)
        writefile(
            filenames[i],
            af.info,
            af.Re,
            af.Mach,
            af.solidity[i],
            af.alpha,
            af.cl[:, i],
            af.cd[:, i],
            radians,
        )
    end
    return nothing
end

"""
    AlphaResolidityAF(alpha, Re, solidity, cl, cd, info)
    AlphaResolidityAF(alpha, Re, solidity, cl, cd)
    AlphaResolidityAF(filenames::Matrix{String}; radians=true)

Airfoil data that varies with angle of attack, Reynolds number, and solidity number.
Data is fit with a recursive Akima spline.

**Arguments**
- `alpha::Vector{Float64}`: angles of attack
- `Re::Vector{Float64}`: Reynolds numbers
- `solidity::Vector{Float64}`: solidity numbers
- `cl::Array{Float64}`: lift coefficients where cl[i, j, k] corresponds to alpha[i], Re[j], solidity[k]
- `cd::Array{Float64}`: drag coefficients where cd[i, j, k] corresponds to alpha[i], Re[j], solidity[k]
- `info::String`: a description of this airfoil data (just informational)

or files with one per Re/solidity combination

**Arguments**
- `filenames::Matrix{String}`: name/path of files to read in.  filenames[i, j] corresponds to Re[i] solidity[j] with Reynolds number and solidity number in ascending order.
- `radians::Bool`: true if angle of attack in file is given in radians
"""
struct AlphaResolidityAF{TF,TS} <: AFType
    alpha::Vector{TF}
    Re::Vector{TF}
    solidity::Vector{TF}
    cl::Array{TF}
    cd::Array{TF}
    info::TS
end

function AlphaResolidityAF(alpha, Re, solidity, cl, cd)
    return AlphaResolidityAF(alpha, Re, solidity, cl, cd, "CCBlade generated airfoil")
end

function AlphaResolidityAF(filenames::Matrix{String}; radians=true)
    info, _, _, _, alpha, _, _ = parsefile(filenames[1, 1], radians)  # assumes common alpha and info across files
    nalpha = length(alpha)
    nRe, nsolidity = size(filenames)

    cl = Array{Float64}(undef, nalpha, nRe, nsolidity)
    cd = Array{Float64}(undef, nalpha, nRe, nsolidity)
    Re = Array{Float64}(undef, nRe)
    solidity = Array{Float64}(undef, nsolidity)

    for j in 1:nsolidity
        for i in 1:nRe
            _, Rei, _, solidityj, _, clij, cdij = parsefile(filenames[i, j], radians)
            cl[:, i, j] = clij
            cd[:, i, j] = cdij
            Re[i] = Rei
            solidity[j] = solidityj  # NOTE: probably should add check to prevent user error here.
        end
    end

    return AlphaResolidityAF(alpha, Re, solidity, cl, cd, info)
end

function afeval(af::AlphaResolidityAF, alpha, Re, solidity)
    cl = FLOWMath.interp3d(
        FLOWMath.akima, af.alpha, af.Re, af.solidity, af.cl, [alpha], [Re], [solidity]
    )[1]
    cd = FLOWMath.interp3d(
        FLOWMath.akima, af.alpha, af.Re, af.solidity, af.cd, [alpha], [Re], [solidity]
    )[1]

    return cl, cd
end

function write_af(filenames, af::AlphaResolidityAF; radians=true)
    nre = length(af.Re)
    nm = length(af.solidity)

    for i in 1:nre
        for j in 1:nm
            writefile(
                filenames[i, j],
                af.info,
                af.Re[i],
                af.Mach,
                af.solidity[j],
                af.alpha,
                af.cl[:, i, j],
                af.cd[:, i, j],
                radians,
            )
        end
    end

    return nothing
end

"""
    AlphaReMachsolidityAF(alpha, Re, Mach, solidity, cl, cd, info)
    AlphaReMachsolidityAF(alpha, Re, Mach, solidity, cl, cd)
    AlphaReMachsolidityAF(filenames::Matrix{String}; radians=true)

Airfoil data that varies with angle of attack, Reynolds number, and Mach number.
Data is fit with a recursive Akima spline.

**Arguments**
- `alpha::Vector{Float64}`: angles of attack
- `Re::Vector{Float64}`: Reynolds numbers
- `Mach::Vector{Float64}`: Mach numbers
- `solidity::Vector{Float64}` : Solidity parameter
- `cl::Array{Float64}`: lift coefficients where cl[i, j, k] corresponds to alpha[i], Re[j], Mach[k]
- `cd::Array{Float64}`: drag coefficients where cd[i, j, k] corresponds to alpha[i], Re[j], Mach[k]
- `info::String`: a description of this airfoil data (just informational)

or files with one per Re/Mach combination

**Arguments**
- `filenames::Matrix{String}`: name/path of files to read in.  filenames[i, j] corresponds to Re[i] Mach[j] with Reynolds number and Mach number in ascending order.
- `radians::Bool`: true if angle of attack in file is given in radians
"""
struct AlphaReMachsolidityAF{TF,TS} <: AFType
    alpha::Vector{TF}
    Re::Vector{TF}
    Mach::Vector{TF}
    solidity::Vector{TF}
    cl::Array{TF}
    cd::Array{TF}
    info::TS
end

function AlphaReMachsolidityAF(alpha, Re, Mach, solidity, cl, cd)
    return AlphaReMachsolidityAF(
        alpha, Re, Mach, solidity, cl, cd, "CCBlade generated airfoil"
    )
end

function AlphaReMachsolidityAF(filenames::Matrix{String}; radians=true)
    info, _, _, _, alpha, _, _ = parsefile(filenames[1, 1], radians)  # assumes common alpha and info across files
    nalpha = length(alpha)
    nRe, nMach, nsolidity = size(filenames)

    cl = Array{Float64}(undef, nalpha, nRe, nMach)
    cd = Array{Float64}(undef, nalpha, nRe, nMach)
    Re = Array{Float64}(undef, nRe)
    Mach = Array{Float64}(undef, nMach)
    solidity = Array{Float64}(undef, nsolidity)

    for k in 1:nsolidity
        for j in 1:nMach
            for i in 1:nRe
                _, Rei, Machj, solidityk, _, clij, cdij = parsefile(
                    filenames[i, j], radians
                )
                cl[:, i, j] = clij
                cd[:, i, j] = cdij
                Re[i] = Rei
                Mach[j] = Machj  # NOTE: probably should add check to prevent user error here.
                solidity[k] = solidityk
            end
        end
    end

    return AlphaReMachsolidityAF(alpha, Re, Mach, solidity, cl, cd, info)
end

function afeval(af::AlphaReMachsolidityAF, alpha, Re, Mach, solidity)
    cl = FLOWMath.interp4d(
        FLOWMath.akima,
        af.alpha,
        af.Re,
        af.Mach,
        af.solidity,
        af.cl,
        [alpha],
        [Re],
        [Mach],
        [solidity],
    )[1]
    cd = FLOWMath.interp4d(
        FLOWMath.akima,
        af.alpha,
        af.Re,
        af.Mach,
        af.solidity,
        af.cd,
        [alpha],
        [Re],
        [Mach],
        [solidity],
    )[1]

    return cl, cd
end

function write_af(filenames, af::AlphaReMachsolidityAF; radians=true)
    nre = length(af.Re)
    nm = length(af.Mach)
    nc = length(af.solidity)

    for k in 1:nc
        for i in 1:nre
            for j in 1:nm
                writefile(
                    filenames[i, j],
                    af.info,
                    af.Re[i],
                    af.Mach[j],
                    af.solidity[k],
                    af.alpha,
                    af.cl[:, i, j],
                    af.cd[:, i, j],
                    radians,
                )
            end
        end
    end

    return nothing
end

"""
    get_clcd(af, alpha, reynolds, mach, solidity)

Return lift and drag coefficients based on airfoil object type and flow conditions.

**Arguments:**
 - `af::AFType` : Airfoil object either of a CCBlade airfoil type or custom type (custom if solidity information is included).
 - `alpha::Float` : angle of attack in degrees
 - `reynolds::Float` : Reynolds number
 - `mach::Float` : Mach number
 - `solidity::Float` : Solidity factor (local chord length over distance between radial blade stations on adjacent blades). Set to `nothing` if not being used.

**Returns:**
 - `cl::Float` : section lift coefficient
 - `cd::Float` : section drag coefficient
"""
function get_clcd(af, alpha, reynolds, mach, solidity)

    #convert alpha to radians
    alpha *= pi / 180.0

    # if solidity undefined, use CCBlade functions
    if solidity == nothing
        return ccb.afeval(af, alpha, reynolds, mach)
    else
        #if solidity is defined, need to use custom functions
        return afeval(af, alpha, reynolds, mach, solidity)
    end
end
