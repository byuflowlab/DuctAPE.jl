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

# - load DuctTAPE - #
using DuctTAPE
const dt = DuctTAPE

# - load plotting defaults - #
include(project_dir * "/visualize/plots_default.jl")
savepath = project_dir * "/test/figures/"
dispath = project_dir

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

######################################################################
#                                                                    #
#                          LIFT CORRECTIONS                          #
#                                                                    #
######################################################################

#---------------------------------#
#  Stagger & Solidity Correction  #
#---------------------------------#

##### ----- Vs AOA Smoothnes ----- #####
# - plot nominal curve - #
# pss = plot(aoa, clo; label="", color=1)

# # - calcualte solidity/stagger correction
# for solidity in range(0.0, 3.0, 10)
#     for stagger in range(-20.0, 100.0, 13)
#         clout, _ = dt.solidityandstagger(clo, solidity, stagger * pi / 180)
# - plot corrected curve - #
# plot!(pss, aoa, clout; label="", color=2)
# end
# end
# plot!(pss, aoa, clo; label="", color=1, xlim=(-20 * pi / 180, 20 * pi / 180))
# # - save figure - #
# savefig(pss, savepath * "soliditystaggercheck.pdf")

##### ----- Vs Solidity Smoothness ----- #####

N = 100
clss = zeros(N, 2)
clssr = zeros(N, 2)
for (is, solidity) in enumerate(range(0.0, 3.0, N))
    for (i, stagger) in enumerate([0.0, pi / 4])
        clss[is, i] = dt.solidityandstaggerfactorsmooth(solidity, stagger)
        clssr[is, i] = dt.solidityandstaggerfactor(solidity, stagger)
    end
end
pss = plot(
    range(0.0, 3.0, N),
    clss[:, 1];
    label="stagger inactive",
    xlabel="solidity",
    ylabel=L"c_\ell",
)
plot!(pss, range(0.0, 3.0, N), clss[:, 2]; label="stagger active")
plot!(pss, range(0.0, 3.0, N), clssr[:, 1]; label="unsmoothed", linestyle=:dash, color=1)
plot!(pss, range(0.0, 3.0, N), clssr[:, 2]; label="unsmoothed", linestyle=:dash, color=2)
# - save figure - #
savefig(pss, savepath * "soliditysmoothnesscheck.pdf")

##### ----- Vs Stagger Smoothness ----- #####

N = 100
clss = zeros(N, 4)
clssr = zeros(N, 4)
pss = plot(; xlabel="stagger angle", ylabel=L"c_\ell")
for (i, solidity) in enumerate([0.01, 1.0, 2.0, 3.0])
    for (is, stagger) in enumerate(range(1, 100, N) * pi / 180)
        clss[is, i] = dt.solidityandstaggerfactorsmooth(solidity, stagger)
        clssr[is, i] = dt.solidityandstaggerfactor(solidity, stagger)
    end
    plot!(
        pss,
        range(1, 100, N),
        clssr[:, i];
        linestyle=:dash,
        color=i,
        label="before smoothing",
    )
    plot!(pss, range(1, 100, N), clss[:, i]; color=i, label="solidity=$(solidity)")
end
# - save figure - #
savefig(pss, savepath * "staggersmoothnesscheck1.pdf")

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

#---------------------------------#
#     transonic limiter correction    #
#---------------------------------#
# - initialize plot - #
# have to create axes manually to make them thinner

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
savefig(pclsm, savepath * "cltranssmoothness.pdf")

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
savefig(pcdsm, savepath * "cdtranssmoothness.pdf")

# - plot nominal curve - #
pml = plot(aoa, clo; label="", color=1)
pdlim = plot(aoa, cdo; label="", color=1)

for mach in range(0.0, 1.0, 11)
    cllim =
        dt.transonicliftlimitersmooth.(
            clo, mach, clcdmin, clmax, clmin, dclda; mcrit=mcrit, verbose=true
        )
    plot!(pml, aoa, cllim; label="", color=myred)
    cdlim = dt.transonicdragaddition(cdo, cllim, clcdmin, mach; mcrit=mcrit, verbose=true)
    plot!(pdlim, aoa, cdlim; label="", color=myred)
end
# plot!(pml, aoa, clo; label="", color=1, xlim=(-20 * pi / 180, 20 * pi / 180))
# plot!(pdlim, aoa, cdo; label="", color=1, xlim=(-20 * pi / 180, 20 * pi / 180))

# - save figure - #
savefig(pml, savepath * "clminmaxlimitcheck.pdf")
savefig(pdlim, savepath * "transdragcheck.pdf")
