project_dir = dirname(dirname(@__FILE__))
if project_dir == ""
    project_dir = "."
end

# TODO: change datapath to test directory or validation data directory
datapath = project_dir * "/validation/data/"
savepath = project_dir * "/validation/data/"
# savepath = project_dir * "/dfdc_airfoil_approximations/"
afname = "dfdc_naca4412"

using LsqFit
using FLOWMath
const fm = FLOWMath
using Plots
using DuctAPE.C4Blade
const c4b = C4Blade

function linfit(aoas, cl0, dclda)
    return dclda .* aoas .+ cl0
end

function quadfit(cls, cdmin, clcdmin, dcddcl2)
    return dcddcl2 * (cls .- clcdmin) .^ 2 .+ cdmin
end

# - Read in Polar - #
info, Re_ref, mach, aoa, cl, cd = c4b.parsefile(datapath * "naca4412.dat", true)
clsp = fm.Akima(aoa, cl)
cdsp = fm.Akima(aoa, cd)

# - Plot Polar - #
plot(aoa * 180 / pi, cl)

# - Chop Polar Down for Fits and replot - #
aoachop = range(-5, 10, 100) * pi / 180
clchop = clsp(aoachop)
cdchop = cdsp(aoachop)
pcl = plot(aoachop, clchop)
pcd = plot(cdchop, clchop)

# - Get cdmin and cl at cdmin as well as zero lift angle of attack- #
cdmin, cdminid = findmin(cdchop)
clcdmin = clchop[cdminid]

_, zeroclid = findmin(abs.(clchop))
alpha0 = fm.linear(
    [clchop[zeroclid - 1], clchop[zeroclid]],
    [aoachop[zeroclid - 1], aoachop[zeroclid]],
    0.0,
)

_, zeroalphaid = findmin(abs.(aoachop))
cl0 = fm.linear(
    [aoachop[zeroalphaid - 1], aoachop[zeroalphaid]],
    [clchop[zeroalphaid - 1], clchop[zeroalphaid]],
    0.0,
)

# - Fit drag curve - #
p0 = [0.05]
quadfitwrap(cls, p) = quadfit(clchop, cdmin, clcdmin, p[1])
cdfit = LsqFit.curve_fit(quadfitwrap, clchop, cdchop, p0)

dcddcl2 = cdfit.param[1]

# - Plot fitted curve - #
plot!(pcd, quadfit(clchop, cdmin, clcdmin, cdfit.param[1]), clchop)

# - Fit Lift Curve - #
p0 = [2.0 * pi]
linfitwrap(aoas, p) = linfit(aoas, cl0, p[1])
clfit = LsqFit.curve_fit(linfitwrap, aoachop, clchop, p0)

dclda = clfit.param[1]

# - Plot fitted curve - #
plot!(pcl, aoachop, linfit(aoachop, cl0, clfit.param[1]))

# - Select clmin, clmax - #
clmin, clminid = findmin(clchop)
clmaxid = findlast(x -> x < 0.1, clfit.resid)
clmax = clchop[clmaxid]

# - Manually set various paramters - #
Re_exp = -0.5
mcrit = 0.7
cmcon = 0.0
dclda_stall = 0.1
dcl_stall = 0.1

# - `alpha0::Float` : zero lift angle of attack
# - `clmax::Float` : maximum cl
# - `clmin::Float` : minimum cl
# - `dclda::Float` : lift curve slope (1/radians)
# - `dclda_stall::Float` :  lift curve slope post-stall (1/radians)
# - `dcl_stall::Float` : cl increment from initial to total stall.
# - `cdmin::Float` : minimum cd
# - `cldmin::Float` : cl at cdmin
# - `dcdcl2::Float` : quadratic curve factor for cd curve (d(cd)/d(cl^2))
# - `cmcon::Float` : pitching moment constant
# - `Re_ref::Float` : reference Reynolds number at which cd values apply
# - `Re_exp::Float` : Reynolds number exponent scaling (cd = cd*(Re/Re_ref)^Re_exp)
# - `mcrit::Float` : critical Mach number

f = open(savepath * afname * ".jl", "w")

write(
    f,
    "alpha0      = $alpha0
    clmax       = $clmax
    clmin       = $clmin
    dclda       = $dclda
    dclda_stall = $dclda_stall
    dcl_stall   = $dcl_stall
    cdmin       = $cdmin
    clcdmin     = $clcdmin
    dcddcl2     = $dcddcl2
    cmcon       = $cmcon
    Re_ref      = $Re_ref
    Re_exp      = $Re_exp
    mcrit       = $mcrit",
)

close(f)
