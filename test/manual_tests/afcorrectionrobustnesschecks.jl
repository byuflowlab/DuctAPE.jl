#=
script for generating plots showing how airfoil polar corrections affect the nominal polars
=#

######################################################################
#                                                                    #
#                                SETUP                               #
#                                                                    #
######################################################################

#---------------------------------#
#       Load Packages, etc.       #
#---------------------------------#

# - Get Project Directory - #
project_dir = dirname(dirname(@__FILE__))
if project_dir == ""
    project_dir = "."
end

# - load DuctAPE - #
using DuctAPE
const dt = DuctAPE

# - load plotting defaults - #
include(project_dir * "/visualize/plots_default.jl")
savepath = project_dir * "/test/figures/"

#---------------------------------#
#    Load Nominal Airfoil Data    #
#---------------------------------#
include(project_dir * "/test/data/naca_4412_raw.jl")
aoa = alpha
clo = cl
cdo = cd

######################################################################
#                                                                    #
#                          Stall Cutoffs                             #
#                                                                    #
######################################################################

aoaext, clext, cdext = dt.stalllimiters(
    aoa, cl, cd; cutoff_slope=0.1, N=20, blend_hardness=50
)

plco = plot(aoaext, clext; label="", color=2)
plot!(plco, aoa, clo; label="", color=1)
# - save figure - #
savefig(plco, savepath * "liftstallcutoffcheck.pdf")

pdco = plot(aoaext, cdext; label="", color=2)
plot!(pdco, aoa, cdo; label="", color=1)
# - save figure - #
savefig(pdco, savepath * "dragstallcutoffcheck.pdf")

aoa = aoaext
clo = clext
cdo = cdext

#---------------------------------#
#    Prandtl-Glauert Correction   #
#---------------------------------#
# - initialize plot - #
# have to create axes manually to make them thinner

# - plot nominal curve - #
pg = plot(aoa, clo; label="", color=1)

# - calcualte P-G correction
for ma in range(0.0, 1.0, 10)
    clpg = dt.prandtlglauert(clo, ma; verbose=true)

    # - plot corrected curve - #
    plot!(pg, aoa, clpg; label="", color=myred, ylim=(minimum(clpg), maximum(clpg) * 1.2))
end

plot!(pg, aoa, clo; label="", color=1, xlim=(-20 * pi / 180, 20 * pi / 180))
# - save figure - #
savefig(pg, savepath * "prandtlglauertcheck.pdf")

N = 1000
clss = zeros(N)
clnom = zeros(N)
pgs = plot(; xlabel="Mach number", ylabel=L"c_\ell")
for (i, ma) in enumerate(range(0.0, 1.5, N))
    clss[i] = dt.prandtlglauert(1.0, ma; verbose=true)
    clnom[i] = 1.0 / sqrt(1.0 - min(ma, 0.99)^2)
end
plot!(pgs, range(0.0, 1.5, N), clnom; label="nominal with cutoff")
plot!(pgs, range(0.0, 1.5, N), clss; linestyle=:dash, linewidth=2, label="smoothed")
# - save figure - #
savefig(pgs, savepath * "pgsmoothnesscheck.pdf")

