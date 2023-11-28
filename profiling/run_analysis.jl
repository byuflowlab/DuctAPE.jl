# - Set up - #
project_dir = dirname(dirname(@__FILE__))
if project_dir == ""
    project_dir = "."
end

datapath = project_dir * "/profiling/"
savepath = project_dir * "/profiling/outputs/"

using DuctAPE
const dt = DuctAPE

# - Include Inputs file - #
include(datapath * "inputs_file.jl")

# - Run Analysis - #
#=
out contains all the post-processed stuff.
converged states are the converged circulation and panel strengths in one big vector (use extract_states to get them separated)
inputs is the tuple created in the precomputation stage and contains all the geometry, etc.
initial_states contains the circulation and panel strengths that were initialized and given to the solver to start with.
convergeflag is just a boolean to tell if the solver converged or not.  If you set verbose to false, this is helpful, if verbose is true, you'll see if it converged.
=#
out, converged_states, inputs, initial_states, convergeflag = dt.analyze_propulsor(
    duct_coordinates,
    hub_coordinates,
    paneling_constants,
    rotor_parameters,
    freestream,
    reference_parameters;
    debug=false, #whatever this did, it's probably broken, so don't set this to true unless you're using it.
    verbose=true, #this basically just sets if NLsolve prints stuff
    maximum_linesearch_step_size=1e6, #this may be worth fiddling with, I haven't touched it since I first decided this was a good value.
    iteration_limit=50, #you should converge in less than 5 iterations for reasonable cases. I've seen some of the CCBlade comparison cases go up to 80ish though.
)

