#=
Types and functions related to airfoil data not included in CCBlade, i.e., airfoil data defined with an additional parameter related to cascades.
=#

"""
overload `parsefile` function from CCBlade
assumes cascade parameter after mach number
"""
function parsefile(filename, radians, cascade)
    alpha = Float64[]
    cl = Float64[]
    cd = Float64[]
    info = ""
    Re = 1.0
    Mach = 1.0
    CAS = 1.0

    open(filename) do f

        # skip header
        info = readline(f)
        Re = parse(Float64, readline(f))
        Mach = parse(Float64, readline(f))
        if cascade
            CAS = parse(Float64, readline(f))
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

    if cascade
        return info, Re, Mach, CAS, alpha, cl, cd
    else
        return info, Re, Mach, alpha, cl, cd
    end
end

"""
overload `writefile` function from CCBlade to include cascade parameter in file header
"""
function writefile(filename, info, Re, Mach, CAS, alpha, cl, cd, radians)
    open(filename, "w") do f
        @printf(f, "%s\n", info)
        @printf(f, "%.17g\n", Re)
        @printf(f, "%.17g\n", Mach)
        @printf(f, "%.17g\n", CAS)

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
Add CCBlade-like airfoil types and functions for cascade parameter
=#
"""
    AlphaCASAF(alpha, CAS, cl, cd, info, Mach)
    AlphaCASAF(alpha, CAS, cl, cd, info)
    AlphaCASAF(alpha, CAS, cl, cd)
    read_AlphaCASAF(filenames::Vector{String}; radians=true)

Airfoil data that varies with angle of attack and Reynolds number.
Data is fit with a recursive Akima spline.

**Arguments**
- `alpha::Vector{Float64}`: angles of attack
- `CAS::Vector{Float64}`: Cascade parameter
- `cl::Matrix{Float64}`: lift coefficients where cl[i, j] corresponds to alpha[i], CAS[j]
- `cd::Matrix{Float64}`: drag coefficients where cd[i, j] corresponds to alpha[i], CAS[j]
- `info::String`: a description of this airfoil data (just informational)
- `Mach::Float64`: Mach number data was taken at (just informational)

or

filenames with one file per Reynolds number.

**Arguments**
- `filenames::Vector{String}`: name/path of files to read in, each at a different Reynolds number in ascending order
- `radians::Bool`: true if angle of attack in file is given in radians
"""
struct AlphaCASAF{TF,TS} <: AFType
    alpha::Vector{TF}
    CAS::Vector{TF}
    cl::Matrix{TF}
    cd::Matrix{TF}
    info::TS # not used except for info in file
    Mach::TF # not used except for info in file
end

AlphaCASAF(alpha, CAS, cl, cd, info) = AlphaCASAF(alpha, CAS, cl, cd, info, 0.0)
function AlphaCASAF(alpha, CAS, cl, cd)
    return AlphaCASAF(alpha, CAS, cl, cd, "CCBlade generated airfoil", 0.0)
end

function AlphaCASAF(filenames::Vector{String}; radians=true)
    info, CAS1, Mach, alpha, cl1, cd1 = parsefile(filenames[1], radians)  # assumes common alpha across files, also common info and common Mach
    nalpha = length(alpha)
    ncond = length(filenames)

    cl = Array{Float64}(undef, nalpha, ncond)
    cd = Array{Float64}(undef, nalpha, ncond)
    CAS = Array{Float64}(undef, ncond)
    cl[:, 1] = cl1
    cd[:, 1] = cd1
    CAS[1] = CAS1

    # iterate over remaining files
    for i in 2:ncond
        _, _, _, CASi, _, cli, cdi = parsefile(filenames[i], radians)
        cl[:, i] = cli
        cd[:, i] = cdi
        CAS[i] = CASi
    end

    return AlphaCASAF(alpha, CAS, cl, cd, info, Mach)
end

function afeval(af::AlphaCASAF, alpha, CAS, Mach)
    cl = FLOWMath.interp2d(FLOWMath.akima, af.alpha, af.CAS, af.cl, [alpha], [CAS])[1]
    cd = FLOWMath.interp2d(FLOWMath.akima, af.alpha, af.CAS, af.cd, [alpha], [CAS])[1]

    return cl, cd
end

function write_af(filenames, af::AlphaCASAF; radians=true)
    for i in 1:length(af.CAS)
        writefile(
            filenames[i],
            af.info,
            af.Re,
            af.Mach,
            af.CAS[i],
            af.alpha,
            af.cl[:, i],
            af.cd[:, i],
            radians,
        )
    end
    return nothing
end

"""
    AlphaReCASAF(alpha, Re, CAS, cl, cd, info)
    AlphaReCASAF(alpha, Re, CAS, cl, cd)
    AlphaReCASAF(filenames::Matrix{String}; radians=true)

Airfoil data that varies with angle of attack, Reynolds number, and CAS number.
Data is fit with a recursive Akima spline.

**Arguments**
- `alpha::Vector{Float64}`: angles of attack
- `Re::Vector{Float64}`: Reynolds numbers
- `CAS::Vector{Float64}`: CAS numbers
- `cl::Array{Float64}`: lift coefficients where cl[i, j, k] corresponds to alpha[i], Re[j], CAS[k]
- `cd::Array{Float64}`: drag coefficients where cd[i, j, k] corresponds to alpha[i], Re[j], CAS[k]
- `info::String`: a description of this airfoil data (just informational)

or files with one per Re/CAS combination

**Arguments**
- `filenames::Matrix{String}`: name/path of files to read in.  filenames[i, j] corresponds to Re[i] CAS[j] with Reynolds number and CAS number in ascending order.
- `radians::Bool`: true if angle of attack in file is given in radians
"""
struct AlphaReCASAF{TF,TS} <: AFType
    alpha::Vector{TF}
    Re::Vector{TF}
    CAS::Vector{TF}
    cl::Array{TF}
    cd::Array{TF}
    info::TS
end

function AlphaReCASAF(alpha, Re, CAS, cl, cd)
    return AlphaReCASAF(alpha, Re, CAS, cl, cd, "CCBlade generated airfoil")
end

function AlphaReCASAF(filenames::Matrix{String}; radians=true)
    info, _, _, _, alpha, _, _ = parsefile(filenames[1, 1], radians)  # assumes common alpha and info across files
    nalpha = length(alpha)
    nRe, nCAS = size(filenames)

    cl = Array{Float64}(undef, nalpha, nRe, nCAS)
    cd = Array{Float64}(undef, nalpha, nRe, nCAS)
    Re = Array{Float64}(undef, nRe)
    CAS = Array{Float64}(undef, nCAS)

    for j in 1:nCAS
        for i in 1:nRe
            _, Rei, _, CASj, _, clij, cdij = parsefile(filenames[i, j], radians)
            cl[:, i, j] = clij
            cd[:, i, j] = cdij
            Re[i] = Rei
            CAS[j] = CASj  # NOTE: probably should add check to prevent user error here.
        end
    end

    return AlphaReCASAF(alpha, Re, CAS, cl, cd, info)
end

function afeval(af::AlphaReCASAF, alpha, Re, CAS)
    cl = FLOWMath.interp3d(
        FLOWMath.akima, af.alpha, af.Re, af.CAS, af.cl, [alpha], [Re], [CAS]
    )[1]
    cd = FLOWMath.interp3d(
        FLOWMath.akima, af.alpha, af.Re, af.CAS, af.cd, [alpha], [Re], [CAS]
    )[1]

    return cl, cd
end

function write_af(filenames, af::AlphaReCASAF; radians=true)
    nre = length(af.Re)
    nm = length(af.CAS)

    for i in 1:nre
        for j in 1:nm
            writefile(
                filenames[i, j],
                af.info,
                af.Re[i],
                af.Mach,
                af.CAS[j],
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
    AlphaReMachCASAF(alpha, Re, Mach, CAS, cl, cd, info)
    AlphaReMachCASAF(alpha, Re, Mach, CAS, cl, cd)
    AlphaReMachCASAF(filenames::Matrix{String}; radians=true)

Airfoil data that varies with angle of attack, Reynolds number, and Mach number.
Data is fit with a recursive Akima spline.

**Arguments**
- `alpha::Vector{Float64}`: angles of attack
- `Re::Vector{Float64}`: Reynolds numbers
- `Mach::Vector{Float64}`: Mach numbers
- `CAS::Vector{Float64}` : cascade parameters
- `cl::Array{Float64}`: lift coefficients where cl[i, j, k] corresponds to alpha[i], Re[j], Mach[k]
- `cd::Array{Float64}`: drag coefficients where cd[i, j, k] corresponds to alpha[i], Re[j], Mach[k]
- `info::String`: a description of this airfoil data (just informational)

or files with one per Re/Mach combination

**Arguments**
- `filenames::Matrix{String}`: name/path of files to read in.  filenames[i, j] corresponds to Re[i] Mach[j] with Reynolds number and Mach number in ascending order.
- `radians::Bool`: true if angle of attack in file is given in radians
"""
struct AlphaReMachCASAF{TF,TS} <: AFType
    alpha::Vector{TF}
    Re::Vector{TF}
    Mach::Vector{TF}
    CAS::Vector{TF}
    cl::Array{TF}
    cd::Array{TF}
    info::TS
end

function AlphaReMachCASAF(alpha, Re, Mach, CAS, cl, cd)
    return AlphaReMachCASAF(alpha, Re, Mach, CAS, cl, cd, "CCBlade generated airfoil")
end

function AlphaReMachCASAF(filenames::Matrix{String}; radians=true)
    info, _, _, _, alpha, _, _ = parsefile(filenames[1, 1], radians)  # assumes common alpha and info across files
    nalpha = length(alpha)
    nRe, nMach, nCAS = size(filenames)

    cl = Array{Float64}(undef, nalpha, nRe, nMach)
    cd = Array{Float64}(undef, nalpha, nRe, nMach)
    Re = Array{Float64}(undef, nRe)
    Mach = Array{Float64}(undef, nMach)
    CAS = Array{Float64}(undef, nCAS)

    for k in 1:nCAS
        for j in 1:nMach
            for i in 1:nRe
                _, Rei, Machj, CASk, _, clij, cdij = parsefile(filenames[i, j], radians)
                cl[:, i, j] = clij
                cd[:, i, j] = cdij
                Re[i] = Rei
                Mach[j] = Machj  # NOTE: probably should add check to prevent user error here.
                CAS[k] = CASk
            end
        end
    end

    return AlphaReMachCASAF(alpha, Re, Mach, CAS, cl, cd, info)
end

function afeval(af::AlphaReMachCASAF, alpha, Re, Mach, CAS)
    cl = FLOWMath.interp4d(
        FLOWMath.akima,
        af.alpha,
        af.Re,
        af.Mach,
        af.CAS,
        af.cl,
        [alpha],
        [Re],
        [Mach],
        [CAS],
    )[1]
    cd = FLOWMath.interp4d(
        FLOWMath.akima,
        af.alpha,
        af.Re,
        af.Mach,
        af.CAS,
        af.cd,
        [alpha],
        [Re],
        [Mach],
        [CAS],
    )[1]

    return cl, cd
end

function write_af(filenames, af::AlphaReMachCASAF; radians=true)
    nre = length(af.Re)
    nm = length(af.Mach)
    nc = length(af.CAS)

    for k in 1:nc
        for i in 1:nre
            for j in 1:nm
                writefile(
                    filenames[i, j],
                    af.info,
                    af.Re[i],
                    af.Mach[j],
                    af.CAS[k],
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
