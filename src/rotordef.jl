#=
Rotor and associated component definition functions
TODO: consider just adding the cascade parameter stuff to what ccblade already does and using CCBlade functionality since it's basically the same.
=#

# for reference, ccblade does:
#
# """
#     AlphaReMachAF(alpha, Re, Mach, cl, cd, info)
#     AlphaReMachAF(alpha, Re, Mach, cl, cd)
#     AlphaReMachAF(filenames::Matrix{String}; radians=true)

# Airfoil data that varies with angle of attack, Reynolds number, and Mach number.
# Data is fit with a recursive Akima spline.

# **Arguments**
# - `alpha::Vector{Float64}`: angles of attack
# - `Re::Vector{Float64}`: Reynolds numbers
# - `Mach::Vector{Float64}`: Mach numbers
# - `cl::Array{Float64}`: lift coefficients where cl[i, j, k] corresponds to alpha[i], Re[j], Mach[k]
# - `cd::Array{Float64}`: drag coefficients where cd[i, j, k] corresponds to alpha[i], Re[j], Mach[k]
# - `info::String`: a description of this airfoil data (just informational)

# or

# files with one per Re/Mach combination

# **Arguments**
# - `filenames::Matrix{String}`: name/path of files to read in.  filenames[i, j] corresponds to Re[i] Mach[j] with Reynolds number and Mach number in ascending order.
# - `radians::Bool`: true if angle of attack in file is given in radians
# """
# struct AlphaReMachAF{TF,TS} <: AFType
#     alpha::Vector{TF}
#     Re::Vector{TF}
#     Mach::Vector{TF}
#     cl::Array{TF}
#     cd::Array{TF}
#     info::TS
# end

# function AlphaReMachAF(alpha, Re, Mach, cl, cd)
#     return AlphaReMachAF(alpha, Re, Mach, cl, cd, "CCBlade generated airfoil")
# end

# function AlphaReMachAF(filenames::Matrix{String}; radians=true)
#     info, _, _, alpha, _, _ = parsefile(filenames[1, 1], radians)  # assumes common alpha and info across files
#     nalpha = length(alpha)
#     nRe, nMach = size(filenames)

#     cl = Array{Float64}(undef, nalpha, nRe, nMach)
#     cd = Array{Float64}(undef, nalpha, nRe, nMach)
#     Re = Array{Float64}(undef, nRe)
#     Mach = Array{Float64}(undef, nMach)

#     for j in 1:nMach
#         for i in 1:nRe
#             _, Rei, Machj, _, clij, cdij = parsefile(filenames[i, j], radians)
#             cl[:, i, j] = clij
#             cd[:, i, j] = cdij
#             Re[i] = Rei
#             Mach[j] = Machj  # NOTE: probably should add check to prevent user error here.
#         end
#     end

#     return AlphaReMachAF(alpha, Re, Mach, cl, cd, info)
# end

# function afeval(af::AlphaReMachAF, alpha, Re, Mach)
#     cl = FLOWMath.interp3d(
#         FLOWMath.akima, af.alpha, af.Re, af.Mach, af.cl, [alpha], [Re], [Mach]
#     )[1]
#     cd = FLOWMath.interp3d(
#         FLOWMath.akima, af.alpha, af.Re, af.Mach, af.cd, [alpha], [Re], [Mach]
#     )[1]

#     return cl, cd
# end

# function write_af(filenames, af::AlphaReMachAF; radians=true)
#     nre = length(af.Re)
#     nm = length(af.Mach)

#     for i in 1:nre
#         for j in 1:nm
#             writefile(
#                 filenames[i, j],
#                 af.info,
#                 af.Re[i],
#                 af.Mach[j],
#                 af.alpha,
#                 af.cl[:, i, j],
#                 af.cd[:, i, j],
#                 radians,
#             )
#         end
#     end

#     return nothing
# end
