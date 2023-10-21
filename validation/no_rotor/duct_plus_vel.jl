#---------------------------------#
#              SETUP              #
#---------------------------------#

# - Get Project Directory - #
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

# create save path
savepath = project_dir * "/validation/no_rotor/figs/"

# - load DuctTAPE - #
using DuctTAPE
const dt = DuctTAPE

# - load plotting defaults - #
include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/plots_default.jl")

# - load geometry - #
# read data file
include(project_dir * "/test/data/naca_662-015_smooth.jl")
# put coordinates together
coordinates = reverse(duct_coordinates; dims=1)

npan = 161
repanel = dt.repanel_airfoil(coordinates; N=npan, normalize=false)

#---------------------------------#
#             Paneling            #
#---------------------------------#
##### ----- Generate Panels ----- #####
# panels = dt.generate_panels([repanel];body=true)
panels = dt.generate_panels([repanel])

xn = panels.node[:, 1]
xcp = panels.controlpoint[:, 1]

#---------------------------------#
#       Operating Conditions      #
#---------------------------------#

# - Add Artificial Velocity - #
weibull(x, k=5, lam=1 / 2) = 1 - exp(-(x / lam)^k)
vartificial = weibull.(panels.controlpoint[1:80, 1])

# Define freestream on panels
Vinf = 1.0 #magnitude doesn't matter yet.
Vs = [1.0 0.0] # axisymmetric, so no radial component
Vsmat = repeat(Vs, size(panels.controlpoint, 1)) # need velocity on each panel
Vsmat[1:80, 1] .+= vartificial
Vsmat .*= Vinf

plot(panels.controlpoint[:, 1], Vsmat[:, 1]; label="", xlabel="x", ylabel="Vinf")
savefig(savepath * "vinf-artificial.pdf")

#---------------------------------#
#        Induced Velocities       #
#---------------------------------#

# - Initial System Matrices - #
AICn, AICt = dt.vortex_aic_boundary_on_boundary(
    panels.controlpoint,
    panels.normal,
    panels.tangent,
    panels.node,
    panels.nodemap,
    panels.influence_length,
)

kids = [
    size(AICn)[1]+1 1
    size(AICn)[1]+1 size(AICn)[2]
]

LHS = zeros(size(AICn)[1] + 1, size(AICn)[2])

dt.add_kutta!(LHS, AICn, kids)

RHS = dt.freestream_influence_vector(panels.normal, Vsmat)
push!(RHS, 0.0)

#---------------------------------#
#             Solving             #
#---------------------------------#
gamb = LHS \ RHS

pg = plot(xn, gamb; xlabel="x", ylabel="node strengths", label="")
savefig(pg, savepath * "duct-gammas.pdf")

#---------------------------------#
#         Post-Processing         #
#---------------------------------#

## --- Velocity Contributions --- ###

# get tangent
Vtan = [dt.dot(v, t) for (v, t) in zip(eachrow(Vsmat), eachrow(panels.tangent))]

# add in body induced tangent velocity
Vtan .+= AICt * gamb

# add in jump term
jump = (gamb[1:(end - 1)] + gamb[2:end]) / 2
Vtan .-= jump / 2.0

### --- Steady Surface Pressure --- ###
cp = 1.0 .- (Vtan / Vinf) .^ 2

f = open(savepath * "ductplusvel-xvcp-$(npan-1)-panels.jl", "w")
write(f, "ductxcp = [\n")
for (x, c) in zip(panels.controlpoint[:, 1], cp)
    write(f, "$x $c\n")
end
write(f, "]")
close(f)

f = open(savepath * "ductplusvel-xvvs-$(npan-1)-panels.jl", "w")
write(f, "ductxvvs = [\n")
for (x, c) in zip(panels.controlpoint[:, 1], Vtan)
    write(f, "$x $c\n")
end
write(f, "]")
close(f)

#---------------------------------#
#             PLOTTING            #
#---------------------------------#
cut = 0

pcp = plot(;
    xlabel="x",
    ylabel=L"c_p",
    yflip=true,
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)
plot!(pcp, xcp[cut+1:end-cut], cp[cut+1:end-cut]; label="", color=myred)
include(savepath * "duct-xvcp-$(npan-1)-panels.jl")
plot!(pcp, ductxcp[:, 1], ductxcp[:, 2]; label="", color=1)
savefig(pcp, savepath * "ductcp-plus-artificialvel.pdf")

pvs = plot(;
    xlabel="x",
    ylabel=L"\frac{V_s}{V_\infty}",
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)
plot!(pvs, xcp[cut+1:end-cut], Vtan[cut+1:end-cut] ./ Vinf; color=myred, label="")
include(savepath * "duct-xvvs-$(npan-1)-panels.jl")
plot!(pvs, ductxvvs[:, 1], ductxvvs[:, 2]; color=myblue, label="")
savefig(pvs, savepath * "ductvel-plus-artificialvel.pdf")

