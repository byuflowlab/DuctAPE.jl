#=
script for generating plots showing how airfoil polar corrections affect the nominal polars
=#

######################################################################
#                                                                    #
#                                SETUP                               #
#                                                                    #
######################################################################

# - Get Project Directory - #
project_dir = dirname(@__FILE__)
if project_dir == ""
    project_dir = "."
end

# create save path
dispath = project_dir * "/"

# - load DuctAPE - #
using DuctAPE
const dt = DuctAPE

# - load plotting defaults - #
# include(project_dir * "/margin_plots_default.jl")
include(project_dir * "/smoothing_plots_default.jl")

#---------------------------------#
#    Load Nominal Airfoil Data    #
#---------------------------------#
include(project_dir * "/naca4412_smoothed_polar.jl")
aoa = ald[:, 1] # in degrees
clo = ald[:, 2]
cdo = ald[:, 3]

######################################################################
#                                                                    #
#                          Stall Cutoffs                             #
#                                                                    #
######################################################################

using Xfoil
using FLOWFoil
x, z = FLOWFoil.naca4(4, 4, 12)
aoarange = range(-25, 25, 51)
clxfoil, cdxfoil, _, _, conv = Xfoil.alpha_sweep(reverse(x), reverse(z), aoarange, 2e6)
good = findall(x -> x == 1, conv)
aoaext, clext, cdext = dt.stalllimiters(
    aoarange[good] * pi / 180,
    clxfoil[good],
    cdxfoil[good];
    cutoff_slope=0.1,
    N=20,
    blend_hardness=50,
)

plco = plot(
    aoarange[good] * pi / 180,
    clxfoil[good];
    color=1,
    label="",
    ylabel=L"c_\ell",
    xlabel="Angle of Attack (radians)",
    xlim=(-pi / 4, pi / 4),
    xticks=(round.([-pi / 4 + 0.1, 0.0, pi / 4 - 0.1], digits=1)),
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)
annotate!(plco, 0.4, 0, text("Nominal", :bottom, myblue, 8))

plot!(plco, aoaext, clext; label="", color=2, linestyle=:dash, linewidth=2)
annotate!(plco, -0.7, 1.9, text("Stall Limited", :left, myred, 8))
# - save figure - #
savefig(plco, dispath * "liftstall-cutoff.tikz")

pdco = plot(
    aoarange[good] * pi / 180,
    cdxfoil[good];
    color=1,
    label="",
    ylabel=L"c_d",
    xlabel="Angle of Attack (radians)",
    xlim=(-pi / 4, pi / 4),
    xticks=(round.([-pi / 4 + 0.1, 0.0, pi / 4 - 0.1], digits=1)),
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)
annotate!(pdco, -0.7, 0.3, text("Nominal", :left, myblue, 8))

plot!(pdco, aoaext, cdext; label="", color=2, linestyle=:dash, linewidth=2)
annotate!(pdco, -0.7, 0.265, text("Stall Limited", :left, myred, 8))
# - save figure - #
savefig(pdco, dispath * "dragstall-cutoff.tikz")

aoa = aoaext
clo = clext
cdo = cdext

######################################################################
#                                                                    #
#                          LIFT CORRECTIONS                          #
#                                                                    #
######################################################################

#---------------------------------#
#  Stagger & Solidity Correction  #
#---------------------------------#

##### ----- Vs Solidity Smoothness ----- #####

N = 100
clss = zeros(N, 2)
clssr = zeros(N, 2)
for (is, solidity) in enumerate(range(0.0, 3.0, N))
    for (i, stagger) in enumerate([0.0, pi / 4])
        clss[is, i] = dt.solidityandstaggerfactorsmooth(solidity, stagger;)
        clssr[is, i] = dt.solidityandstaggerfactor(solidity, stagger;)
    end
end

pss = plot(
    range(0.0, 3.0, N),
    clss[:, 1];
    label="",
    xlabel="Solidity",
    ylabel=L"c_{\ell_\mathrm{ss}}",
)
annotate!(pss, 1.5, 0.9, text("Nominal", :left, myblue, 8))
# plot!(pss, range(0.0, 3.0, N), clss[:, 2]; label="stagger active")
plot!(pss, range(0.0, 3.0, N), clssr[:, 1]; label="", linestyle=:dash, color=2, linewidth=2)
annotate!(pss, 1.5, 0.8, text("Smoothed", :left, myred, 8))
# plot!(pss, range(0.0, 3.0, N), clssr[:, 2]; label="unsmoothed", linestyle=:dash, color=2)
# - save figure - #
savefig(pss, dispath * "solidity-smoothed.tikz")

##### ----- Vs Stagger Smoothness ----- #####

N = 100
clss = zeros(N, 1)
clssr = zeros(N, 1)
pss = plot(; xlabel="Stagger (degrees)", ylabel=L"c_{\ell_\mathrm{ss}}")
for (i, solidity) in enumerate([2.0])
    for (is, stagger) in enumerate(range(1, 100, N) * pi / 180)
        clss[is, i] = dt.solidityandstaggerfactorsmooth(solidity, stagger;)
        clssr[is, i] = dt.solidityandstaggerfactor(solidity, stagger;)
    end
    plot!(pss, range(1, 100, N), clssr[:, i]; color=1, label="")
    plot!(
        pss, range(1, 100, N), clss[:, i]; color=2, linestyle=:dash, linewidth=2, label=""
    )
end
annotate!(pss, 5, 1.0, text("Nominal", :left, myblue, 8))
annotate!(pss, 5, 0.9, text("Smoothed", :left, myred, 8))
# - save figure - #
savefig(pss, dispath * "stagger-smoothed.tikz")

#---------------------------------#
#    Prandtl-Glauert Correction   #
#---------------------------------#

N = 100
clss = zeros(N)
clnom = zeros(N)
pgs = plot(; xlabel="Mach number", ylabel=L"c_{\ell_\text{pg}}")
for (i, ma) in enumerate(range(0.0, 1.5, N))
    clss[i] = dt.prandtlglauert(1.0, ma; verbose=true)
    clnom[i] = 1.0 / sqrt(1.0 - min(ma, 0.99)^2)
end
plot!(pgs, range(0.0, 1.5, N), clnom; label="")
annotate!(pgs, 0.5, 6, text("Nominal", :bottom, myblue, 6))
plot!(pgs, range(0.0, 1.5, N), clss; linestyle=:dash, linewidth=2, label="")
annotate!(pgs, 0.5, 6, text("Smoothed", :top, myred, 6))
# - save figure - #
savefig(pgs, dispath * "pg-smoothed-margin.tikz")

#---------------------------------#
#     transonic limiter correction    #
#---------------------------------#
# - calcualte transonic limiter correction
mcrit = 0.5
clcdmin = clo[findmin(cdo)[2]]
clmax, maxid = findmax(clo)
clmin, minid = findmin(clo)
dclda = (clmax - clmin) / (aoa[maxid] * pi / 180 - aoa[minid] * pi / 180)

N = 500
clnom = zeros(N)
clsmooth = zeros(N)
for (i, mach) in enumerate(range(0.0, 1.0, N))
    clnom[i] = dt.transonicliftlimiter(
        1.0, mach, clcdmin, clmax, clmin, dclda; mcrit=mcrit, verbose=true
    )
    clsmooth[i] = dt.transonicliftlimitersmooth(
        1.0, mach, clcdmin, clmax, clmin, dclda; mcrit=mcrit, verbose=true
    )
end

pclsm = plot(
    range(0, 1, N),
    clnom;
    label="",
    color=1,
    xlabel="Mach Number",
    ylabel=L"c_{\ell_\mathrm{lim}}",
)
plot!(pclsm, range(0, 1, N), clsmooth; label="", color=2, linestyle=:dash, linewidth=2)
annotate!(pclsm, 0.05, 0.95, text("Nominal", :left, myblue, 8))
annotate!(pclsm, 0.05, 0.9, text("Smoothed", :left, myred, 8))

# - save figure - #
savefig(pclsm, dispath * "cltranslim-smoothed.tikz")

cdnom = zeros(N)
cdsmooth = zeros(N)
for (i, mach) in enumerate(range(0.0, 1.0, N))
    cl = dt.transonicliftlimitersmooth(
        1.0, mach, clcdmin, clmax, clmin, dclda; mcrit=mcrit, verbose=true
    )
    cdnom[i] = dt.transonicdragaddition(1.0, cl, clcdmin, mach; mcrit=mcrit, verbose=true)
    cdsmooth[i] = dt.transonicdragadditionsmooth(
        1.0, cl, clcdmin, mach; mcrit=mcrit, verbose=true
    )
end

pcdsm = plot(
    range(0, 1, N),
    cdnom;
    label="",
    color=1,
    xlabel="Mach Number",
    ylabel=L"c_{d_\mathrm{lim}}",
)
plot!(pcdsm, range(0, 1, N), cdsmooth; label="", color=2, linestyle=:dash, linewidth=2)
annotate!(pcdsm, 0.05, 2.5, text("Nominal", :left, myblue, 8))
annotate!(pcdsm, 0.05, 2.3, text("Smoothed", :left, myred, 8))

# - save figure - #
savefig(pcdsm, dispath * "cdtranslim-smoothed.tikz")

