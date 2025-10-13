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
project_dir = dirname(@__FILE__)
if project_dir == ""
    project_dir = "."
end

# create save path
dispath = project_dir

# - load DuctTAPE - #
using DuctTAPE
const dt = DuctTAPE

# - load plotting defaults - #
include(project_dir * "/margin_plots_default.jl")

#---------------------------------#
#    Load Nominal Airfoil Data    #
#---------------------------------#
include(project_dir * "/naca4412_smoothed_polar.jl")
aoa = ald[:, 1] # in degrees
clo = ald[:, 2]
cdo = ald[:, 3]

######################################################################
#                                                                    #
#                          LIFT CORRECTIONS                          #
#                                                                    #
######################################################################

#---------------------------------#
#  Stagger & Solidity Correction  #
#---------------------------------#
# - initialize plot - #
# have to create axes manually to make them thinner
pss = plot([minimum(aoa); maximum(aoa)], [0.0; 0.0]; linewidth=0.25, color=:black, label="")
plot!(pss, [0.0; 0.0], [minimum(clo); maximum(clo)]; linewidth=0.25, color=:black, label="")

# - plot nominal curve - #
plot!(pss, aoa, clo; label="", color=1)
annotate!(
    pss,
    aoa[end]-4,
    maximum(clo),
    text("Nominal", :right, myblue, 6),
)

# - calcualte P-G correction
solidity = 1.0
stagger = pi / 4.0
clss = dt.solidityandstagger(clo, solidity,stagger)

# - plot corrected curve - #
plot!(pss, aoa, clss; label="", color=2)
annotate!(
    pss,
    aoa[end]+2,
    maximum(clss)/2,
    text("Corrected", :right, myred, 6),
)

# - save figure - #
savefig(pss, "soliditystagger-correction-margin.pdf")
savefig(pss, "soliditystagger-correction-margin.tikz")

#---------------------------------#
#    Prandtl-Glauert Correction   #
#---------------------------------#
# - initialize plot - #
# have to create axes manually to make them thinner
pg = plot([minimum(aoa); maximum(aoa)], [0.0; 0.0]; linewidth=0.25, color=:black, label="")
plot!(pg, [0.0; 0.0], [minimum(clo); maximum(clo)]; linewidth=0.25, color=:black, label="")

# - plot nominal curve - #
plot!(pg, aoa, clo; label="", color=1)
annotate!(
    pg,
    aoa[findfirst(a -> a > maximum(clo) / 2 + 0.1, clo)],
    maximum(clo) / 2,
    text("Nominal", :left, myblue, 6),
)

# - calcualte P-G correction
ma = 0.5
clpg = dt.prandtlglauert(clo, ma)

# - plot corrected curve - #
plot!(pg, aoa, clpg; label="", color=myred, ylim=(minimum(clpg), maximum(clpg) * 1.2))
annotate!(
    pg,
    aoa[findfirst(a -> a > maximum(clpg) - 0.2, clpg)],
    maximum(clo) + 0.2,
    text("Corrected", :right, myred, 6),
)

# - save figure - #
savefig(pg, "prandtlglauert-correction-margin.tikz")

#---------------------------------#
# transonic cl limiter correction #
#---------------------------------#
# - initialize plot - #
# have to create axes manually to make them thinner
pml = plot([minimum(aoa); maximum(aoa)], [0.0; 0.0]; linewidth=0.25, color=:black, label="")
plot!(pml, [0.0; 0.0], [minimum(clo); maximum(clo)]; linewidth=0.25, color=:black, label="")

# - plot nominal curve - #
plot!(pml, aoa, clo; label="", color=1)
annotate!(pml, aoa[findmax(clo)[2]] - 3, maximum(clo), text("Nominal", :right, myblue, 6))

# - calcualte transonic limiter correction
ma = 0.8
clcdmin = clo[findmin(cdo)[2]]
clmax, maxid = findmax(clo)
clmin, minid = findmin(clo)
dclda = (clmax - clmin) / (aoa[maxid] * pi / 180 - aoa[minid] * pi / 180)
cllim = dt.transonicliftlimiter.(clo, ma, clcdmin, clmax, clmin, dclda)

# - plot corrected curve - #
plot!(pml, aoa, cllim; label="", color=myred)
annotate!(
    pml,
    aoa[findfirst(a -> a > maximum(clo) / 2, clo)],
    maximum(clo) / 2 - 0.1,
    text("Limited", :left, myred, 6),
)

# - save figure - #
savefig(pml, "clminmaxlimit-correction-margin.tikz")

######################################################################
#                                                                    #
#                          DRAG CORRECTIONS                          #
#                                                                    #
######################################################################

#---------------------------------#
#    Reynolds Number Correction   #
#---------------------------------#

# - initialize plot - #
# have to create axes manually to make them thinner
pdre = plot(
    [minimum(aoa); maximum(aoa)], [0.0; 0.0]; linewidth=0.25, color=:black, label=""
)
plot!(pdre, [0.0; 0.0], [0.0; maximum(cdo)]; linewidth=0.25, color=:black, label="")

# - plot nominal curve - #
plot!(pdre, aoa, cdo; label="", color=1)
annotate!(
    pdre, aoa[findmax(cdo)[2]], minimum(cdo) + 0.05, text("Nominal", :right, myblue, 6)
)

# - calcualte Re correction
reref = 2e6
re = 5e6
reexp = 0.5
cdre = dt.redrag(cdo, re, reref; reexp=reexp)

# - plot corrected curve - #
plot!(pdre, aoa, cdre; label="", color=myred)
annotate!(pdre, aoa[findmax(cdo)[2]], minimum(cdre), text("Adjusted", :right, myred, 6))

# - save figure - #
savefig(pdre, "redrag-correction-margin.tikz")

#---------------------------------#
#    Transonic Drag Correction    #
#---------------------------------#

# - initialize plot - #
# have to create axes manually to make them thinner
pdlim = plot(
    [minimum(aoa); maximum(aoa)], [0.0; 0.0]; linewidth=0.25, color=:black, label=""
)
plot!(pdlim, [0.0; 0.0], [0.0; maximum(cdo)]; linewidth=0.25, color=:black, label="")

# - plot nominal curve - #
plot!(pdlim, aoa, cdo; label="", color=1)
annotate!(
    pdlim, aoa[findmax(cdo)[2]], minimum(cdo) + 0.05, text("Nominal", :right, myblue, 6)
)

# - calcualte transonic correction
mach = 0.8
mcrit = 0.7
cdlim = dt.transonicdragaddition(cdo, cllim, clcdmin, mach; mcrit=mcrit)

plot!(pdlim, [0.0; 0.0], [0.0; maximum(cdlim)]; linewidth=0.25, color=:black, label="")
# - plot corrected curve - #
plot!(pdlim, aoa, cdlim; label="", color=myred)
annotate!(
    pdlim, aoa[findmax(cdo)[2]] + 2, maximum(cdlim) / 2, text("Augmented", :right, myred, 6)
)

# - save figure - #
savefig(pdlim, "transdrag-correction-margin.tikz")

