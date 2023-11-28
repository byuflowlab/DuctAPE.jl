project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

savepath = project_dir*"/test/visual_tests/"

using DuctAPE
const dt = DuctAPE

include(project_dir * "/visualize/plots_default_new.jl")
include(project_dir * "/visualize/visualize_geometry.jl")

# panels and bodies
# define coordinates
x1 = [2.0; 0.5; 0.0; 0.5; 2.0]
r1 = [2.0; 1.5; 2.0; 2.5; 2.0]

x2 = [0.0; 0.5; 2.0]
r2 = [0.0; 0.5; 0.0]

c1 = [x1 r1]
c2 = [x2 r2]

coordinates = [c1, c2]

panels = dt.generate_panels(coordinates; body=true)

visualize_paneling(;
    body_panels=panels,
    # rotor_panels=nothing,
    # wake_panels=nothing,
    coordinates=coordinates,
    controlpoints=true,
    nodes=true,
    # wakeinterfaceid=[],
    # prescribedpanels=nothing,
    normals=true,
    normal_scaling=0.1,
    savepath=savepath,
    filename=["constant_vortex_paneling_check.pdf"],
    legendloc=:best,
    zoom=false,
    limits=nothing,
    nodemarkersize=2,
    cpmarkersize=2,
)

