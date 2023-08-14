using Test

println("Compiling Package")
using DuctTAPE
const dt = DuctTAPE
using ForwardDiff
const frd = ForwardDiff
using FiniteDiff
const fnd = FiniteDiff
using FLOWMath
const fm = FLOWMath
# using ImplicitAD

println("Running Tests...")
include("compare_objects.jl")

# - Active Development - #
# include("derivative_checks.jl")
include("velocity_probe_tests.jl")

# # - Should be Working - #
# include("test_cascade_functions.jl")
# include("blade_element_aero_lookups.jl")
# include("new_panel_tests.jl")
# include("basic_singularity.jl")
# include("post_processing_tests.jl")




# - Need to update, add, fix, etc. - #
# include("body_aero_tests.jl")
# include("body_geometry_tests.jl")
# # TODO: need to add source and body induced velocities to induced velocity on rotor test
# include("rotor_aero_tests.jl")
# # TODO: need to update most of this test after having added sources and bodies in.
# include("dimension_tests.jl")

# # BROKEN:
# include("constant_initialization.jl")
# include("wake_aero_tests.jl")

# # include("aero_coefficient_tests.jl")
# # include("basic_geometry.jl")
# # include("geometry_tests.jl")
# # include("rotor_post_processing_tests.jl")
