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
coordinates = [x_hub r_hub]

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
# TODO: write new LHS function that includes LU decomposition and returns LU decomposed object
LHSLU = dt.vortex_panel_influence_matrix(panels.nodes, panels)
# put this in the vortex panel influence matrix stuff
# LHSLU = lu!(LHS, NoPivot(); check=false) # we shouldn't need a pivot since we have self-induced velocities, and we don't want to throw an error if the factorization didn't work
# lufail = !issuccess(LHSLU) # we can pass this as part of the optimization fail flag

# note that this is not the body strengths, but rather the system RHS which will be overwritten to be the body strengths upon solving in place.
gamb = dt.freestream_influence_vector(panels.normal, Vsmat)

#---------------------------------#
#             Solving             #
#---------------------------------#
# - Use LinearAlgebra's ldiv! to solve linear system using factorized matrix
# note that gamb is overwritten to be the RHS vector and then overwritten to be the body strengths
la.ldiv!(LHSLU, gamb)

#---------------------------------#
#         Post-Processing         #
#---------------------------------#

### --- Velocity Contributions --- ###
# - Body-induced Surface Velocity - #
# TODO: write new surface velocity calucation function.  may need to name this something else as it likely already exists. also probably want to make sure it does what you think it should
Vb = dt.vfromvortexpanels(panels.controlpoint, panels.controlpoint, gamb)

#---------------------------------#
#             PLOTTING            #
#---------------------------------#
pp = plot(; xlabel="x", ylabel=L"V_s/V_\infty")
plot!(
    pp,
    Vs_over_Vinf_x,
    Vs_over_Vinf_vs;
    seriestype=:scatter,
    color=myblue[1],
    markershape=:utriangle,
    label="experimental",
)
xs = panels.controlpoint[:, 1]
plot!(pp, xs, dt.norm.(eachrow(Vtot)) ./ Vinf; label="DuctTAPE")

savefig(savepath * "hub-vel-comp.pdf")
