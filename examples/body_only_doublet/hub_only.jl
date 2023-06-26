#---------------------------------#
#              SETUP              #
#---------------------------------#

# - Get Project Directory - #
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

# create save path
savepath = project_dir*"/examples/body_only_doublet/"

# - load DuctTAPE - #
using DuctTAPE
const dt = DuctTAPE

# - load plotting defaults - #
include(project_dir*"/visualize/visualize_geometry.jl")
include(project_dir*"/visualize/plots_default_new.jl")

# - load geometry - #
# read data file
include(project_dir*"/test/data/bodyofrevolutioncoords.jl")
# hub final r-coordinate needs to be set to zero so that it's not negative
r_hub[end]=0.0
# put coordinates together
hub_coordinates = [x_hub r_hub]

#---------------------------------#
#             Paneling            #
#---------------------------------#
##### ----- Generate Panels ----- #####
#= use new Neumann paneling functions in panel.jl
generate_panels(coordinates::Vector{Matrix{TF}}) where {TF}
note that input must be a vector of matrices, so even if only modeling one body, still
needs to be a vector.
Also note, multiple dispatch enabled to allow for single body, simply calls the function
placing the coordinate matrix in a vector as an input.
=#
hub_panels = dt.generate_panels(hub_coordinates)

##### ----- Visualize to Check ----- #####
visualize_paneling(hub_panels; coordinates=hub_coordinates, control_points=true, nodes=true, normals=true, savepath=savepath, filename="hub-geometry.pdf")

#---------------------------------#
#        Induced Velocities       #
#---------------------------------#

# TODO: need to coordinate with Eduardo on who is adding what where for this stuff
