#---------------------------------#
#              SETUP              #
#---------------------------------#

# - Get Project Directory - #
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

# create save path
savepath = project_dir * "/validation/no_rotor/"

# - load DuctTAPE - #
using DuctTAPE
const dt = DuctTAPE

# - load plotting defaults - #
include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/plots_default_new.jl")

# - load experimental data - #
include(project_dir * "/test/data/naca_662-015.jl")

# - load geometry - #
# read data file
include(project_dir * "/test/data/naca_662-015_smooth.jl")
# put coordinates together
coordinates = reverse(duct_coordinates; dims=1)
npan = 5000
coordinates = dt.repanel_airfoil(coordinates; N=npan, normalize=false)

#---------------------------------#
#             Paneling            #
#---------------------------------#
##### ----- Generate Panels ----- #####
panels = dt.generate_panels([coordinates]; body=true)

##### ----- Visualize to Check ----- #####
visualize_paneling(;
    body_panels=panels,
    coordinates=[coordinates],
    controlpoints=true,
    nodes=true,
    normals=true,
    savepath=savepath,
    filename=["duct-geometry.pdf"],
)

xn = panels.node[:, 1]
xcp = panels.controlpoint[:, 1]

#---------------------------------#
#       Operating Conditions      #
#---------------------------------#

# Define freestream on panels
Vinf = 1.0 #magnitude doesn't matter yet.
Vs = Vinf * [1.0 0.0] # axisymmetric, so no radial component
Vsmat = repeat(Vs, size(panels.controlpoint, 1)) # need velocity on each panel

#---------------------------------#
#        Induced Velocities       #
#---------------------------------#

# - Initial System Matrices - #
LHS = dt.vortex_influence_matrix(
    panels.controlpoint, panels.normal, panels.node, panels.influence_length
)
RHS = dt.freestream_influence_vector(panels.normal, Vsmat)

#---------------------------------#
#             Solving             #
#---------------------------------#
gamb = LHS \ RHS

cut = max(15, round(Int, npan / 50))
pg = plot(
    xn[(cut + 1):(end - cut)],
    gamb[(cut + 1):(end - cut)];
    xlabel="x",
    ylabel="panel strengths",
    label="",
)
savefig(pg, savepath * "duct-gammas.pdf")

#---------------------------------#
#         Post-Processing         #
#---------------------------------#

### --- Velocity Contributions --- ###
Vb = similar(Vsmat) .= 0.0
dt.vfromvortices!(Vb, panels.controlpoint, panels.node, panels.influence_length, gamb)
Vb
Vb .+= Vsmat

# get vz at nodes
# why does this not work when multiplying by influence lengths, but half works when you don't?
vznode =
    abs.([
        dt.smoke_ring_vz(r, l) * g for
        (r, l, g) in zip(panels.node[:, 2], panels.influence_length, gamb)
    ])

# average nodes to get values at cps
vzcp = (vznode[2:end] .+ vznode[1:(end - 1)]) / 2
# get vz at TE cp
push!(vzcp, (vznode[1] + vznode[end]) / 2)
vzcp = [vzcp zeros(length(vzcp))]

# add half of vzcp at each cp
# Vb .+= vzcp / 2

# get tangent
Vtan = [dt.dot(v, t) for (v, t) in zip(eachrow(Vb), eachrow(panels.tangent))]

### --- Steady Surface Pressure --- ###
cp = 1.0 .- (Vtan / Vinf) .^ 2

#---------------------------------#
#             PLOTTING            #
#---------------------------------#
pp = plot(; xlabel="x", ylabel=L"c_p", yflip=true)
plot!(
    pp,
    pressurexupper,
    pressureupper;
    seriestype=:scatter,
    color=myblue[1],
    markershape=:utriangle,
    label="exp outer",
)
plot!(
    pp,
    pressurexlower,
    pressurelower;
    seriestype=:scatter,
    color=myblue[1],
    markershape=:dtriangle,
    label="exp inner",
)
# cut = 15
plot!(pp, xcp[(cut + 1):(end - cut)], cp[(cut + 1):(end - cut)]; label="DuctAPE")

savefig(savepath * "duct-pressure-comp.pdf")
