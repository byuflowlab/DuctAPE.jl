using Test

println("Compiling Package")
using DuctAPE
const dt = DuctAPE

using FLOWMath
const fm = FLOWMath

using StaticArrays
const sa = StaticArrays

using PreallocationTools

using ForwardDiff
const frd = ForwardDiff
# using ReverseDiff
# const rd = ReverseDiff
using FiniteDiff
const fnd = FiniteDiff

include("test_utils.jl")

println("Running Tests...")

# - pre-process related tests - #
include("afcorrections.jl")
include("panel_generation_tests.jl")
include("induced_velocities.jl")
include("influence_coefficients.jl")
include("linear_system_assembly.jl")
include("pre_processing_tests.jl")

# - solve related tests - #
include("iteration_step_tests.jl")
include("state_estimation.jl")
include("relaxation_tests.jl")
include("wake_aero_tests.jl")

# - post process related tests - #
include("post_processing_tests.jl")
include("thermodynamics_tests.jl")
include("boundary_layer_tests.jl")

##########################################################
#EVERYTHING BELOW THIS POINT NEEDS TO BE UPDATED
##########################################################

# - Active Development - #
# include("solve_checks.jl")

# - Need to update, add, fix, etc. - #
# include("body_aero_tests.jl")
# # TODO: need to add source and body induced velocities to induced velocity on rotor test
# include("rotor_aero_tests.jl")
# TODO: need to update function names and argument order for functions as well a tests for them
# include("test_cascade_functions.jl")
# # include("coupled_geometry.jl")
# include("derivative_checks.jl")
