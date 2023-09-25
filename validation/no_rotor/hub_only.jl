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
include(project_dir * "/visualize/plots_default_new.jl")
include(project_dir * "/visualize/visualize_geometry.jl")

# - load geometry - #
# read data file
include(project_dir * "/test/data/bodyofrevolutioncoords.jl")
# hub final r-coordinate needs to be set to zero so that it's not negative
r_hub[end] = 0.0
# put coordinates together
coordinates = [x_hub[1:end-1] r_hub[1:end-1]]

#---------------------------------#
#             Paneling            #
#---------------------------------#
##### ----- Generate Panels ----- #####
panels = dt.generate_panels(coordinates)

##### ----- Visualize to Check ----- #####
visualize_paneling(;
    body_panels=panels,
    coordinates=[coordinates],
    controlpoints=true,
    nodes=true,
    TEnodes=true,
    normals=true,
    savepath=savepath,
    filename=["hub-geometry.pdf"],
)

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
LHS = dt.vortex_panel_influence_matrix(panels, panels)
# put this in the vortex panel influence matrix stuff
# LHSLU = dt.lu!(LHS, dt.NoPivot(); check=false) # we shouldn't need a pivot since we have self-induced velocities, and we don't want to throw an error if the factorization didn't work
# lufail = !dt.issuccess(LHSLU) # we can pass this as part of the optimization fail flag

# note that this is not the body strengths, but rather the system RHS which will be overwritten to be the body strengths upon solving in place.
# gamb = dt.freestream_influence_vector(panels.normal, Vsmat)
RHS = dt.freestream_influence_vector(panels.normal, Vsmat)

#---------------------------------#
#             Solving             #
#---------------------------------#
# - Use LinearAlgebra's ldiv! to solve linear system using factorized matrix
# note that gamb is overwritten to be the RHS vector and then overwritten to be the body strengths
# dt.ldiv!(LHSLU, gamb)
gamb = LHS\RHS

#---------------------------------#
#         Post-Processing         #
#---------------------------------#

### --- Velocity Contributions --- ###
# - Body-induced Surface Velocity - #
# TODO: write new surface velocity calculation function.  may need to name this something else as it likely already exists. also probably want to make sure it does what you think it should
Vb = similar(Vsmat) .= 0.0
dt.vfromvortexpanels!(
    Vb, panels.controlpoint, panels.controlpoint, panels.len, gamb
)
Vb .+= Vsmat
Vtan = [dt.dot(v,t) for (v,t) in zip(eachrow(Vb), eachrow(panels.tangent))]
dt.norm.(eachrow(Vb))

# get x-coordinates for plotting
xs = panels.controlpoint[:, 1]
#---------------------------------#
#             PLOTTING            #
#---------------------------------#
pg = plot(; xlabel="x", ylabel=L"\mathrm{Panel~strengths~}(\gamma^B)")
plot!(pg,xs,gamb,label="")
savefig(savepath * "hub-gammas.pdf")

pp = plot(; xlabel="x", ylabel=L"\mathrm{Normalized~surface~velocity~}(V_s/V_\infty)")
plot!(
    pp,
    Vs_over_Vinf_x,
    Vs_over_Vinf_vs;
    seriestype=:scatter,
    color=myblue[1],
    markershape=:utriangle,
    label="experimental",
)
plot!(pp, xs, Vtan ./ Vinf; label="DuctTAPE")

savefig(savepath * "hub-vel-comp.pdf")
