project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

using DuctTAPE
const dt = DuctTAPE

datapath = project_dir * "/examples/dfdc_partial/"
savepath = datapath * "outputs/"

include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/visualize_flowfield.jl")
include(project_dir * "/visualize/plots_default_new.jl")
# pyplot()

# read in parameters file
include(datapath * "lewis_with_rotor.casenasa.jl")

# run analyze_propulsor function
out, converged_states, inputs, initial_states, convergeflag = @time dt.analyze_propulsor(
    duct_coordinates,
    hub_coordinates,
    paneling_constants,
    rotor_parameters,
    freestream,
    reference_parameters;
    debug=false,
    verbose=true,
    maximum_linesearch_step_size=1e6,
    iteration_limit=50,
)

## -- Plotting -- ##

@time "Visualize 2D" begin
    println("Visualizing Paneling")
    # check normals
    visualize_paneling(;
        body_panels=inputs.body_doublet_panels,
        coordinates=[duct_coordinates, hub_coordinates],
        controlpoints=true,
        nodes=true,
        normals=true,
        normal_scaling=0.1,
        savepath=savepath,
        filename=["lc-bodygeometry.pdf"; "lc-bodygeomtry.png"],
        legendloc=:right,
    )

    # everything
    visualize_paneling(;
        body_panels=inputs.body_doublet_panels,
        rotor_panels=inputs.rotor_source_panels,
        wake_panels=inputs.wake_vortex_panels,
        coordinates=[duct_coordinates, hub_coordinates],
        controlpoints=true,
        nodes=false,
        normals=false,
        normal_scaling=0.1,
        savepath=savepath,
        filename=["lc-fullgeometry.pdf"; "lc-fullgeometry.png"],
        legendloc=:outerright,
    )
end
