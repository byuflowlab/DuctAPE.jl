project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

savedir = project_dir*"/examples/new_panels/"

using DuctTAPE
const dt = DuctTAPE
include(project_dir*"/visualize/plots_default_new.jl")

# define coordinates
x1 = [1.0; 0.5; 0.0; 0.5; 1.0]
r1 = [2.0; 1.5; 2.0; 2.5; 2.0]

x2 = [0.0; 0.5; 1.0]
r2 = [0.0; 0.5; 0.0]

c1 = [x1 r1]
c2 = [x2 r2]

coordinates = [c1,c2]

# generate panels
panels = dt.generate_panels(coordinates)

# plot generated panels
plot(aspectratio=1, xlabel="x", ylabel="r")

# plot inputs
plot!(x1,r1, color=:black, linestyle=:dot, label="input coordinates")
plot!(x2,r2, color=:black, linestyle=:dot, label="")

#plot control points
plot!(panels.controlpoint[:,1], panels.controlpoint[:,2], color=red[2], seriestype=:scatter, markershape=:rect, label="Control Points")

#plot nodes
for i in 1:length(panels.len)
    lab= i==1 ? "Nodes" : ""
    plot!(panels.nodes[i,:,1],panels.nodes[i,:,2], label=lab, color=blue[2], seriestype=:scatter)
end

#plot normal
for i in 1:length(panels.len)
    lab= i==1 ? "Normals" : ""
plot!([0.0;panels.normal[i,1]].+panels.controlpoint[i,1],[0.0;panels.normal[i,2]].+panels.controlpoint[i,2], label=lab, color=blue[3])
end

savefig(savedir*"/checkpaneling.pdf")
