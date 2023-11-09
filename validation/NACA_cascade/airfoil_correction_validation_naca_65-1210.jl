#=
Comparison against NACA cascade experiments with Airfoil polars and corrections applied
TODO: move this file into the test or docs folder along with other plotted items
=#

#---------------------------------#
#              SETUP              #
#---------------------------------#

# - Get Project Directory - #
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

# create save path
savepath = project_dir * "/validation/NACA_cascade/figures/"
dispath =
    project_dir * "/../../Writing/dissertation/src/ductsolvercontents/ductsolverfigures/"

using DuctTAPE
const dt = DuctTAPE
using FLOWMath
using Xfoil

# - load plotting defaults - #
include(project_dir * "/visualize/plots_default.jl")

# - load experimental data - #
include(project_dir * "/test/data/NACA_cascade_data_65-1210.jl")

#---------------------------------#
#           Load Geometry         #
#---------------------------------#

# - Get Airfoil Coordinates - #
x, z = dt.naca65c(1.2)

# save coordinates for tikz
f = open(dispath * "naca-651210scaled.dat", "w")
for (x, z) in zip(eachrow(x), eachrow(z))
    write(f, "$(x[1]) $(z[1])\n")
end
close(f)

# plot for double checking
plot(x, z; aspectratio=1, label="", axis=false, ticks=false)
savefig(savepath * "naca-651210scaled.pdf")

# Reverse Data for Xfoil
reverse!(x)
reverse!(z)

#---------------------------------#
#       Get Airfoil Polars        #
#---------------------------------#

anglesofattack = range(-15, 50; step=1.0)

# - Get Xfoil data at Re=245000 - #
cl245, cd245, _, _, conv245 = Xfoil.alpha_sweep(
    x,
    z,
    anglesofattack,
    245000;# reinit=true, percussive_maintenance=true
)
gid245 = findall(c -> c == 1, conv245)
a245g = anglesofattack[gid245]
cl245g = cl245[gid245]
cd245g = cd245[gid245]

# spline data
cl245sp = FLOWMath.Akima(a245g, cl245g)
cd245sp = FLOWMath.Akima(a245g, cd245g)

# plot to check things
ppl245 = plot(
    a245g,
    cl245g;
    label="Xfoil Data",
    xlabel="Angle of Attack",
    ylabel=L"c_\ell",
    title="Re = 245,000",
    xlim=(minimum(a245g), maximum(a245g)),
)
ppd245 = plot(
    a245g,
    cd245g;
    label="Xfoil Data",
    xlabel="Angle of Attack",
    ylabel=L"c_d",
    title="Re=245,000",
    xlim=(minimum(a245g), maximum(a245g)),
)

# - Get Xfoil data at Re=200000 - #
cl200, cd200, _, _, conv200 = Xfoil.alpha_sweep(
    x,
    z,
    anglesofattack,
    200000; #reinit=true, percussive_maintenance=true
)
gid200 = findall(c -> c == 1, conv200)
a200g = anglesofattack[gid200]
cl200g = cl200[gid200]
cd200g = cd200[gid200]

# spline data
cl200sp = FLOWMath.Akima(a200g, cl200g)
cd200sp = FLOWMath.Akima(a200g, cd200g)

# plot to check things
ppl200 = plot(
    a200g,
    cl200g;
    label="Xfoil Data",
    xlabel="Angle of Attack",
    ylabel=L"c_\ell",
    title="Re=200,000",
    xlim=(minimum(a200g), maximum(a200g)),
)
ppd200 = plot(
    a200g,
    cd200g;
    label="Xfoil Data",
    xlabel="Angle of Attack",
    ylabel=L"c_d",
    title="Re=200,000",
    xlim=(minimum(a200g), maximum(a200g)),
)

#---------------------------------#
#   Plot Experimental Lift Data   #
#---------------------------------#
##### ----- Beta = 30 ----- #####
pbl30 = plot(; xlabel="Angle of Attack (degrees)", ylabel=L"c_\ell", title=L"\beta=30")
plot!(
    pbl30,
    cl_b30s100[:, 1],
    cl_b30s100[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.0",
)
plot!(
    pbl30,
    cl_b30s125[:, 1],
    cl_b30s125[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.25",
)
plot!(
    pbl30,
    cl_b30s150[:, 1],
    cl_b30s150[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.5",
)

##### ----- Beta = 45 ----- #####
pbl45 = plot(; xlabel="Angle of Attack (degrees)", ylabel=L"c_\ell", title=L"\beta=45")
plot!(
    pbl45,
    cl_b45s050[:, 1],
    cl_b45s050[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=0.5",
)
plot!(
    pbl45,
    cl_b45s075[:, 1],
    cl_b45s075[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=0.75",
)
plot!(
    pbl45,
    cl_b45s100[:, 1],
    cl_b45s100[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.0",
)
plot!(
    pbl45,
    cl_b45s125[:, 1],
    cl_b45s125[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.25",
)
plot!(
    pbl45,
    cl_b45s150[:, 1],
    cl_b45s150[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.5",
)

##### ----- Beta = 60 ----- #####
pbl60 = plot(; xlabel="Angle of Attack (degrees)", ylabel=L"c_\ell", title=L"\beta=60")
plot!(
    pbl60,
    cl_b60s050[:, 1],
    cl_b60s050[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=0.5",
)
plot!(
    pbl60,
    cl_b60s075[:, 1],
    cl_b60s075[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=0.75",
)
plot!(
    pbl60,
    cl_b60s100[:, 1],
    cl_b60s100[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.0",
)
plot!(
    pbl60,
    cl_b60s125[:, 1],
    cl_b60s125[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.25",
)
plot!(
    pbl60,
    cl_b60s150[:, 1],
    cl_b60s150[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.5",
)

##### ----- Beta = 70 ----- #####
pbl70 = plot(; xlabel="Angle of Attack (degrees)", ylabel=L"c_\ell", title=L"\beta=70")
plot!(
    pbl70,
    cl_b70s100[:, 1],
    cl_b70s100[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.0",
)
plot!(
    pbl70,
    cl_b70s125[:, 1],
    cl_b70s125[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.25",
)
plot!(
    pbl70,
    cl_b70s150[:, 1],
    cl_b70s150[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.5",
)

#---------------------------------#
#   Plot Experimental Drag Data   #
#---------------------------------#
##### ----- Beta = 30 ----- #####
pbd30 = plot(; xlabel="Angle of Attack (degrees)", ylabel=L"c_d", title=L"\beta=30")
plot!(
    pbd30,
    cd_b30s100[:, 1],
    cd_b30s100[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.0",
)
plot!(
    pbd30,
    cd_b30s125[:, 1],
    cd_b30s125[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.25",
)
plot!(
    pbd30,
    cd_b30s150[:, 1],
    cd_b30s150[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.5",
)

##### ----- Beta = 45 ----- #####
pbd45 = plot(; xlabel="Angle of Attack (degrees)", ylabel=L"c_d", title=L"\beta=45")
plot!(
    pbd45,
    cd_b45s050[:, 1],
    cd_b45s050[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=0.5",
)
plot!(
    pbd45,
    cd_b45s075[:, 1],
    cd_b45s075[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=0.75",
)
plot!(
    pbd45,
    cd_b45s100[:, 1],
    cd_b45s100[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.0",
)
plot!(
    pbd45,
    cd_b45s125[:, 1],
    cd_b45s125[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.25",
)
plot!(
    pbd45,
    cd_b45s150[:, 1],
    cd_b45s150[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.5",
)

##### ----- Beta = 60 ----- #####
pbd60 = plot(; xlabel="Angle of Attack (degrees)", ylabel=L"c_d", title=L"\beta=60")
plot!(
    pbd60,
    cd_b60s050[:, 1],
    cd_b60s050[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=0.5",
)
plot!(
    pbd60,
    cd_b60s075[:, 1],
    cd_b60s075[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=0.75",
)
plot!(
    pbd60,
    cd_b60s100[:, 1],
    cd_b60s100[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.0",
)
plot!(
    pbd60,
    cd_b60s125[:, 1],
    cd_b60s125[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.25",
)
plot!(
    pbd60,
    cd_b60s150[:, 1],
    cd_b60s150[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.5",
)

##### ----- Beta = 70 ----- #####
pbd70 = plot(; xlabel="Angle of Attack (degrees)", ylabel=L"c_d", title=L"\beta=70")
plot!(
    pbd70,
    cd_b70s100[:, 1],
    cd_b70s100[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.0",
)
plot!(
    pbd70,
    cd_b70s125[:, 1],
    cd_b70s125[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.25",
)
plot!(
    pbd70,
    cd_b70s150[:, 1],
    cd_b70s150[:, 2];
    seriestype=:scatter,
    markershape=:utriangle,
    label=L"\sigma=1.5",
)

######################################################################
#                                                                    #
#                         Apply Corrections                          #
#                                                                    #
######################################################################

# - Estimate Mach number - #
mach = 95.0 / 1125.33 #mach based on flow speed in feet per second

# - Estimate Polar Features for Re = 245,000 - #
clmin245, clminid245 = findmin(cl245g)
clmax245, clmaxid245 = findmax(cl245g[1:25])
clcdmin245 = cl245g[findmin(cd245g)[2]]
dclda245 = (clmax245 - clmin245) / ((a245g[clmaxid245] - a245g[clminid245]) * pi / 180.0)

# - stall limit data for Re = 245,000 - #
a245ext, cl245ext, cd245ext = dt.stalllimiters(
    a245g * pi / 180.0,
    cl245g,
    cd245g;
    cutoff_slope=0.1,
    N=20,
    blend_hardness=50,
    clminid=clminid245,
    clmaxid=clmaxid245,
)
plot!(
    ppl245,
    a245ext * 180 / pi,
    cl245ext;
    label="stall limited data",
    linestyle=:dash,
    linewidth=2,
)
plot!(
    ppd245,
    a245ext * 180 / pi,
    cd245ext;
    label="stall limited data",
    linestyle=:dash,
    linewidth=2,
)

# - Estimate Polar Features for Re = 200,000 - #
clmin200, clminid200 = findmin(cl200g)
clmax200, clmaxid200 = findmax(cl200g[1:23])
clcdmin200 = cl200g[findmin(cd200g)[2]]
dclda200 = (clmax200 - clmin200) / ((a200g[clmaxid200] - a200g[clminid245]) * pi / 180.0)

# - stall limit data for Re = 200,000 - #
a200ext, cl200ext, cd200ext = dt.stalllimiters(
    a200g * pi / 180.0,
    cl200g,
    cd200g;
    cutoff_slope=0.1,
    N=20,
    blend_hardness=50,
    clminid=clminid200,
    clmaxid=clmaxid200,
)
plot!(
    ppl200,
    a200ext * 180 / pi,
    cl200ext;
    label="stall limited data",
    linestyle=:dash,
    linewidth=2,
)
plot!(
    ppd200,
    a200ext * 180 / pi,
    cd200ext;
    label="stall limited data",
    linestyle=:dash,
    linewidth=2,
)

savefig(ppl245, savepath * "naca_651210_cl_re245000.pdf")
savefig(ppd245, savepath * "naca_651210_cd_re245000.pdf")
savefig(ppl200, savepath * "naca_651210_cl_re200000.pdf")
savefig(ppd200, savepath * "naca_651210_cd_re200000.pdf")
#---------------------------------#
#         Beta = 30 Cases         #
#---------------------------------#

##### ----- Solidity = 1.0 ----- #####
cl_b30_s100_r245 = FLOWMath.akima(a245ext, cl245ext, cl_b30s100[:, 1] * pi / 180)
cd_b30_s100_r245 = FLOWMath.akima(a245ext, cd245ext, cd_b30s100[:, 1] * pi / 180)

for (a, cl, cd) in
    zip(cl_b30s100[:, 1], eachrow(cl_b30_s100_r245), eachrow(cd_b30_s100_r245))
    stagger = (30.0 - a) * pi / 180.0
    solidity = 1.0
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin245, clmax245, clmin245, dclda245
    )
end
plot!(pbl30, cl_b30s100[:, 1], cl_b30_s100_r245; label="1.0 af corr", color=1)
plot!(pbd30, cd_b30s100[:, 1], cd_b30_s100_r245; label="1.0 af corr", color=1)

##### ----- Solidity = 1.25 ----- #####
cl_b30_s125_r245 = FLOWMath.akima(a245ext, cl245ext, cl_b30s125[:, 1] * pi / 180)
cd_b30_s125_r245 = FLOWMath.akima(a245ext, cd245ext, cd_b30s125[:, 1] * pi / 180)

for (a, cl, cd) in
    zip(cl_b30s125[:, 1], eachrow(cl_b30_s125_r245), eachrow(cd_b30_s125_r245))
    stagger = (30.0 - a) * pi / 180.0
    solidity = 1.25
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin245, clmax245, clmin245, dclda245
    )
end
plot!(pbl30, cl_b30s125[:, 1], cl_b30_s125_r245; label="1.25 af corr", color=2)
plot!(pbd30, cd_b30s125[:, 1], cd_b30_s125_r245; label="1.25 af corr", color=2)

##### ----- Solidity = 1.5 ----- #####
cl_b30_s150_r245 = FLOWMath.akima(a245ext, cl245ext, cl_b30s150[:, 1] * pi / 180)
cd_b30_s150_r245 = FLOWMath.akima(a245ext, cd245ext, cd_b30s150[:, 1] * pi / 180)
for (a, cl, cd) in
    zip(cl_b30s150[:, 1], eachrow(cl_b30_s150_r245), eachrow(cd_b30_s150_r245))
    stagger = (30.0 - a) * pi / 180.0
    solidity = 1.5
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin245, clmax245, clmin245, dclda245
    )
end
plot!(pbl30, cl_b30s150[:, 1], cl_b30_s150_r245; label="1.5 af corr", color=3)
plot!(pbd30, cd_b30s150[:, 1], cd_b30_s150_r245; label="1.5 af corr", color=3)

savefig(pbl30, savepath * "naca_651210_cl_b30_comp.pdf")
savefig(pbd30, savepath * "naca_651210_cd_b30_comp.pdf")

#---------------------------------#
#         Beta = 45 Cases         #
#---------------------------------#

##### ----- Solidity = 0.5 ----- #####
cl_b45_s050_r200 = FLOWMath.akima(a200ext, cl200ext, cl_b45s050[:, 1] * pi / 180)
cd_b45_s050_r200 = FLOWMath.akima(a200ext, cd200ext, cd_b45s050[:, 1] * pi / 180)

for (a, cl, cd) in
    zip(cl_b45s050[:, 1], eachrow(cl_b45_s050_r200), eachrow(cd_b45_s050_r200))
    stagger = (45.0 - a) * pi / 180.0
    solidity =0.5
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin200, clmax200, clmin200, dclda200
    )
end
plot!(pbl45, cl_b45s050[:, 1], cl_b45_s050_r200; label="0.5 af corr", color=1)
plot!(pbd45, cd_b45s050[:, 1], cd_b45_s050_r200; label="0.5 af corr", color=1)

##### ----- Solidity = 0.75 ----- #####
cl_b45_s075_r200 = FLOWMath.akima(a200ext, cl200ext, cl_b45s075[:, 1] * pi / 180)
cd_b45_s075_r200 = FLOWMath.akima(a200ext, cd200ext, cd_b45s075[:, 1] * pi / 180)

for (a, cl, cd) in
    zip(cl_b45s075[:, 1], eachrow(cl_b45_s075_r200), eachrow(cd_b45_s075_r200))
    stagger = (45.0 - a) * pi / 180.0
    solidity =0.75
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin200, clmax200, clmin200, dclda200
    )
end
plot!(pbl45, cl_b45s075[:, 1], cl_b45_s075_r200; label="0.75 af corr", color=2)
plot!(pbd45, cd_b45s075[:, 1], cd_b45_s075_r200; label="0.75 af corr", color=2)

##### ----- Solidity = 1.0 ----- #####
cl_b45_s100_r245 = FLOWMath.akima(a245ext, cl245ext, cl_b45s100[:, 1] * pi / 180)
cd_b45_s100_r245 = FLOWMath.akima(a245ext, cd245ext, cd_b45s100[:, 1] * pi / 180)

for (a, cl, cd) in
    zip(cl_b45s100[:, 1], eachrow(cl_b45_s100_r245), eachrow(cd_b45_s100_r245))
    stagger = (45.0 - a) * pi / 180.0
    solidity = 1.0
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin245, clmax245, clmin245, dclda245
    )
end
plot!(pbl45, cl_b45s100[:, 1], cl_b45_s100_r245; label="1.0 af corr", color=3)
plot!(pbd45, cd_b45s100[:, 1], cd_b45_s100_r245; label="1.0 af corr", color=3)

##### ----- Solidity = 1.25 ----- #####
cl_b45_s125_r245 = FLOWMath.akima(a245ext, cl245ext, cl_b45s125[:, 1] * pi / 180)
cd_b45_s125_r245 = FLOWMath.akima(a245ext, cd245ext, cd_b45s125[:, 1] * pi / 180)

for (a, cl, cd) in
    zip(cl_b45s125[:, 1], eachrow(cl_b45_s125_r245), eachrow(cd_b45_s125_r245))
    stagger = (45.0 - a) * pi / 180.0
    solidity = 1.25
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin245, clmax245, clmin245, dclda245
    )
end
plot!(
    pbl45,
    cl_b45s125[:, 1],
    cl_b45_s125_r245;
    label="1.25 af corr",
    color=4,
    linestyle=:dash,
)
plot!(
    pbd45,
    cd_b45s125[:, 1],
    cd_b45_s125_r245;
    label="1.25 af corr",
    color=4,
    linestyle=:dash,
)

##### ----- Solidity = 1.5 ----- #####
cl_b45_s150_r245 = FLOWMath.akima(a245ext, cl245ext, cl_b45s150[:, 1] * pi / 180)
cd_b45_s150_r245 = FLOWMath.akima(a245ext, cd245ext, cd_b45s150[:, 1] * pi / 180)
for (a, cl, cd) in
    zip(cl_b45s150[:, 1], eachrow(cl_b45_s150_r245), eachrow(cd_b45_s150_r245))
    stagger = (45.0 - a) * pi / 180.0
    solidity = 1.5
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin245, clmax245, clmin245, dclda245
    )
end
plot!(
    pbl45, cl_b45s150[:, 1], cl_b45_s150_r245; label="1.5 af corr", color=5, linestyle=:dash
)
plot!(
    pbd45, cd_b45s150[:, 1], cd_b45_s150_r245; label="1.5 af corr", color=5, linestyle=:dash
)

savefig(pbl45, savepath * "naca_651210_cl_b45_comp.pdf")
savefig(pbd45, savepath * "naca_651210_cd_b45_comp.pdf")

#---------------------------------#
#         Beta = 60 Cases         #
#---------------------------------#

##### ----- Solidity = 0.5 ----- #####
cl_b60_s050_r200 = FLOWMath.akima(a200ext, cl200ext, cl_b60s050[:, 1] * pi / 180)
cd_b60_s050_r200 = FLOWMath.akima(a200ext, cd200ext, cd_b60s050[:, 1] * pi / 180)

for (a, cl, cd) in
    zip(cl_b60s050[:, 1], eachrow(cl_b60_s050_r200), eachrow(cd_b60_s050_r200))
    stagger = (60.0 - a) * pi / 180.0
    solidity = 0.5
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin200, clmax200, clmin200, dclda200
    )
end
plot!(pbl60, cl_b60s050[:, 1], cl_b60_s050_r200; label="0.5 af corr", color=1)
plot!(pbd60, cd_b60s050[:, 1], cd_b60_s050_r200; label="0.5 af corr", color=1)

##### ----- Solidity = 0.75 ----- #####
cl_b60_s075_r200 = FLOWMath.akima(a200ext, cl200ext, cl_b60s075[:, 1] * pi / 180)
cd_b60_s075_r200 = FLOWMath.akima(a200ext, cd200ext, cd_b60s075[:, 1] * pi / 180)

for (a, cl, cd) in
    zip(cl_b60s075[:, 1], eachrow(cl_b60_s075_r200), eachrow(cd_b60_s075_r200))
    stagger = (60.0 - a) * pi / 180.0
    solidity =0.75
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin200, clmax200, clmin200, dclda200
    )
end
plot!(pbl60, cl_b60s075[:, 1], cl_b60_s075_r200; label="0.75 af corr", color=2)
plot!(pbd60, cd_b60s075[:, 1], cd_b60_s075_r200; label="0.75 af corr", color=2)

##### ----- Solidity = 1.0 ----- #####
cl_b60_s100_r245 = FLOWMath.akima(a245ext, cl245ext, cl_b60s100[:, 1] * pi / 180)
cd_b60_s100_r245 = FLOWMath.akima(a245ext, cd245ext, cd_b60s100[:, 1] * pi / 180)

for (a, cl, cd) in
    zip(cl_b60s100[:, 1], eachrow(cl_b60_s100_r245), eachrow(cd_b60_s100_r245))
    stagger = (60.0 - a) * pi / 180.0
    solidity = 1.0
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin245, clmax245, clmin245, dclda245
    )
end
plot!(pbl60, cl_b60s100[:, 1], cl_b60_s100_r245; label="1.0 af corr", color=3)
plot!(pbd60, cd_b60s100[:, 1], cd_b60_s100_r245; label="1.0 af corr", color=3)

##### ----- Solidity = 1.25 ----- #####
cl_b60_s125_r245 = FLOWMath.akima(a245ext, cl245ext, cl_b60s125[:, 1] * pi / 180)
cd_b60_s125_r245 = FLOWMath.akima(a245ext, cd245ext, cd_b60s125[:, 1] * pi / 180)

for (a, cl, cd) in
    zip(cl_b60s125[:, 1], eachrow(cl_b60_s125_r245), eachrow(cd_b60_s125_r245))
    stagger = (60.0 - a) * pi / 180.0
    solidity = 1.25
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin245, clmax245, clmin245, dclda245
    )
end
plot!(
    pbl60,
    cl_b60s125[:, 1],
    cl_b60_s125_r245;
    label="1.25 af corr",
    color=4,
    linestyle=:dash,
)
plot!(
    pbd60,
    cd_b60s125[:, 1],
    cd_b60_s125_r245;
    label="1.25 af corr",
    color=4,
    linestyle=:dash,
)

##### ----- Solidity = 1.5 ----- #####
cl_b60_s150_r245 = FLOWMath.akima(a245ext, cl245ext, cl_b60s150[:, 1] * pi / 180)
cd_b60_s150_r245 = FLOWMath.akima(a245ext, cd245ext, cd_b60s150[:, 1] * pi / 180)
for (a, cl, cd) in
    zip(cl_b60s150[:, 1], eachrow(cl_b60_s150_r245), eachrow(cd_b60_s150_r245))
    stagger = (60.0 - a) * pi / 180.0
    solidity = 1.5
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin245, clmax245, clmin245, dclda245
    )
end
plot!(
    pbl60, cl_b60s150[:, 1], cl_b60_s150_r245; label="1.5 af corr", linestyle=:dash, color=5
)
plot!(
    pbd60, cd_b60s150[:, 1], cd_b60_s150_r245; label="1.5 af corr", linestyle=:dash, color=5
)

savefig(pbl60, savepath * "naca_651210_cl_b60_comp.pdf")
savefig(pbd60, savepath * "naca_651210_cd_b60_comp.pdf")

#---------------------------------#
#         Beta = 70 Cases         #
#---------------------------------#

##### ----- Solidity = 1.0 ----- #####
cl_b70_s100_r245 = FLOWMath.akima(a245ext, cl245ext, cl_b70s100[:, 1] * pi / 180)
cd_b70_s100_r245 = FLOWMath.akima(a245ext, cd245ext, cd_b70s100[:, 1] * pi / 180)

for (a, cl, cd) in
    zip(cl_b70s100[:, 1], eachrow(cl_b70_s100_r245), eachrow(cd_b70_s100_r245))
    stagger = (70.0 - a) * pi / 180.0
    solidity = 1.0
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin245, clmax245, clmin245, dclda245
    )
end
plot!(pbl70, cl_b70s100[:, 1], cl_b70_s100_r245; label="1.0 af corr", color=1)
plot!(pbd70, cd_b70s100[:, 1], cd_b70_s100_r245; label="1.0 af corr", color=1)

##### ----- Solidity = 1.25 ----- #####
cl_b70_s125_r245 = FLOWMath.akima(a245ext, cl245ext, cl_b70s125[:, 1] * pi / 180)
cd_b70_s125_r245 = FLOWMath.akima(a245ext, cd245ext, cd_b70s125[:, 1] * pi / 180)

for (a, cl, cd) in
    zip(cl_b70s125[:, 1], eachrow(cl_b70_s125_r245), eachrow(cd_b70_s125_r245))
    stagger = (70.0 - a) * pi / 180.0
    solidity = 1.25
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin245, clmax245, clmin245, dclda245
    )
end
plot!(pbl70, cl_b70s125[:, 1], cl_b70_s125_r245; label="1.25 af corr", color=2)
plot!(pbd70, cd_b70s125[:, 1], cd_b70_s125_r245; label="1.25 af corr", color=2)

##### ----- Solidity = 1.5 ----- #####
cl_b70_s150_r245 = FLOWMath.akima(a245ext, cl245ext, cl_b70s150[:, 1] * pi / 180)
cd_b70_s150_r245 = FLOWMath.akima(a245ext, cd245ext, cd_b70s150[:, 1] * pi / 180)
for (a, cl, cd) in
    zip(cl_b70s150[:, 1], eachrow(cl_b70_s150_r245), eachrow(cd_b70_s150_r245))
    stagger = (70.0 - a) * pi / 180.0
    solidity = 1.5
    dt.corrected_clcd!(
        cl, cd, mach, solidity, stagger, clcdmin245, clmax245, clmin245, dclda245
    )
end
plot!(pbl70, cl_b70s150[:, 1], cl_b70_s150_r245; label="1.5 af corr", color=3)
plot!(pbd70, cd_b70s150[:, 1], cd_b70_s150_r245; label="1.5 af corr", color=3)

savefig(pbl70, savepath * "naca_651210_cl_b70_comp.pdf")
savefig(pbd70, savepath * "naca_651210_cd_b70_comp.pdf")
