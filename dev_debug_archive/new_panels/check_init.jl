project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

datapath = project_dir * "/dev_debug_archive/new_panels/"
savepath = datapath

include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/plots_default_new.jl")

include(datapath * "run_body_only.jl")

using DuctAPE
const dt = DuctAPE

using FLOWMath
const fm = FLOWMath

function init_inputs(; datapath="", filename="test.case.jl")
    println("Reading Parameters")
    include(datapath * filename)

    println("Initializing Inputs")
    # initialize various inputs used in analysis
    inputs = dt.precomputed_inputs(
        duct_coordinates,
        hub_coordinates,
        paneling_constants,
        rotor_parameters,
        freestream,
        reference_parameters;
        debug=false,
    )

    return inputs
end

function init_states(inputs)
    println("Initializing States")
    return dt.initialize_states(inputs)
end

function check_mu_init(inputs)
    states = init_states(inputs)

    mub, gamw, Gamr, sigr = dt.extract_state_variables(states, inputs)

    panels = inputs.body_doublet_panels

    mu = run_body_only(panels, inputs.Vinf, inputs.prescribedpanels)

    println("mu_init = ", mub)
    println("mu_body_only = ", mu)
    println("mu_init - mu_body_only: ")
    display(mub .- mu)

    return nothing
end

function plot_inputs(inputs)
    println("\tplotting body")
    # plot body panels
    visualize_paneling(;
        body_panels=inputs.body_doublet_panels,
        controlpoints=true,
        nodes=true,
        normals=true,
        normal_scaling=0.02,
        prescribedpanels=inputs.prescribedpanels,
        savepath=savepath,
        filename="test-body-geometry.pdf",
    )

    println("\tplotting rotors")
    # plot rotor panels
    visualize_paneling(;
        rotor_panels=inputs.rotor_source_panels,
        controlpoints=true,
        nodes=true,
        normals=false,
        savepath=savepath,
        filename="test-rotor-geometry.pdf",
    )

    println("\tplotting wake")
    # plot wake panels
    visualize_paneling(;
        wake_panels=inputs.wake_vortex_panels,
        controlpoints=true,
        nodes=true,
        wakeinterfaceid=reduce(
            vcat, [inputs.ductwakeinterfaceid; inputs.hubwakeinterfaceid]
        ),
        normals=false,
        savepath=savepath,
        filename="test-wake-geometry.pdf",
    )

    println("\tplotting everything")
    # plot everythings panels
    visualize_paneling(;
        body_panels=inputs.body_doublet_panels,
        rotor_panels=inputs.rotor_source_panels,
        wake_panels=inputs.wake_vortex_panels,
        controlpoints=true,
        nodes=true,
        normals=false,
        wakeinterfaceid=reduce(
            vcat, [inputs.ductwakeinterfaceid, inputs.hubwakeinterfaceid]
        ),
        prescribedpanels=inputs.prescribedpanels,
        savepath=savepath,
        filename="test-all-geometry.pdf",
    )

    return nothing
end

# inputs = init_inputs(; datapath=datapath,filename="test.case..jl")
inputs = init_inputs(;
    datapath="", filename=project_dir * "/test/data/basic_two_rotor_for_test.jl"
)
println("Plotting")
plot_inputs(inputs)
check_mu_init(inputs)
# states = init_states(inputs)
# out = dt.post_process(states, inputs)
;
