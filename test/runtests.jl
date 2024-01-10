using Test

println("Compiling Package")
using DuctAPE
const dt = DuctAPE
# using CCBlade
# const ccb = CCBlade
using FLOWMath
const fm = FLOWMath
using ForwardDiff
const frd = ForwardDiff
using FiniteDiff
const fnd = FiniteDiff
using StaticArrays
const sa = StaticArrays

include("test_utils.jl")

println("Running Tests...")

# - Active Development - #
# include("derivative_checks.jl")

# - Should be Working - #
include("afcorrections.jl")
include("panel_generation_tests.jl")
include("induced_velocities.jl")
include("influence_coefficients.jl")
include("linear_system_assembly.jl")
include("pre_processing_tests.jl")
include("relaxation_tests.jl")

# - add back in after body is validated and need to get rotor/wake working again.
# include("blade_element_aero_lookups.jl")
# include("post_processing_tests.jl")

# - Not necessary for now - #
# include("test_cascade_functions.jl")

# - Need to update, add, fix, etc. - #
# include("body_aero_tests.jl")
# # TODO: need to add source and body induced velocities to induced velocity on rotor test
# include("rotor_aero_tests.jl")
# # TODO: need to update most of this test after having added sources and bodies in.
# include("dimension_tests.jl")

# # BROKEN:
# include("constant_initialization.jl")
# include("wake_aero_tests.jl")
# # include("aero_coefficient_tests.jl")
# # include("coupled_geometry.jl")
# # include("rotor_post_processing_tests.jl")
