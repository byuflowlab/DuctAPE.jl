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

#---------------------------------#
#             Paneling            #
#---------------------------------#
##### ----- Generate Panels ----- #####
panels = dt.generate_panels(coordinates)

##### ----- Generate Control Points and Normals ----- #####
using LinearAlgebra
cpnormal = zeros(size(coordinates, 1) - 2, 2)
for i in 1:(panels.npanels - 1)
    p1 = panels.controlpoint[i, :]
    p2 = panels.controlpoint[i + 1, :]
    that = (p1 .- p2) / dt.norm(p1 .- p2)
    that = [that; 0.0]
    cpnormal[i, :] = LinearAlgebra.cross(that,[0.0;0.0;1.0])[1:2]
end
cptuple = (;
    controlpoint=coordinates[2:(end - 1), :], normal=cpnormal, npanels=size(cpnormal, 1)
)

##### ----- Visualize to Check ----- #####
visualize_paneling(;
    # body_panels=panels,
    body_panels=cptuple,
    coordinates=nothing,
    controlpoints=true,
    nodes=false,
    normals=true,
    savepath=savepath,
    filename="duct-geometry.pdf",
)

xs = panels.controlpoint[:, 1]

#---------------------------------#
#       Operating Conditions      #
#---------------------------------#

# Define freestream on panels
Vinf = 1.0 #magnitude doesn't matter yet.
Vs = Vinf * [1.0 0.0] # axisymmetric, so no radial component
Vsmat = repeat(Vs, panels.npanels) # need velocity on each panel

#---------------------------------#
#        Induced Velocities       #
#---------------------------------#

# - Initial System Matrices - #
LHS = dt.vortex_panel_influence_matrix(panels, cptuple)
RHS = dt.freestream_influence_vector(cptuple.normal, Vsmat)

# - Adding Kutta Condition - #
LHS = [LHS; zeros(size(LHS, 2))']
LHS[end, 1] = 1.0
LHS[end, end] = 1.0
push!(RHS, 0.0)

#---------------------------------#
#             Solving             #
#---------------------------------#
gamb = LHS \ RHS

pg = plot(xs, gamb; xlabel="x", ylabel="panel strengths", label="")
savefig(pg, savepath * "duct-gammas.pdf")

#---------------------------------#
#         Post-Processing         #
#---------------------------------#

### --- Velocity Contributions --- ###
Vb = similar(Vsmat) .= 0.0
dt.vfromvortexpanels!(Vb, panels.controlpoint, panels.controlpoint, panels.len, gamb)
Vb .+= Vsmat
Vtan = [dt.dot(v, t) for (v, t) in zip(eachrow(Vb), eachrow(panels.tangent))]
dt.norm.(eachrow(Vb))

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
cut = 0
plot!(pp, xs[(cut + 1):(end - cut)], cp[(cut + 1):(end - cut)]; label="DuctAPE")

savefig(savepath * "duct-pressure-comp.pdf")
