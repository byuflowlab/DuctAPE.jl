using Test

println("Compiling Package")
using DuctTAPE
const dt = DuctTAPE
using CCBlade
const ccb = CCBlade
using FLOWFoil
const ff = FLOWFoil
using ForwardDiff
const frd = ForwardDiff
using FiniteDiff
const fnd = FiniteDiff
using FLOWMath
const fm = FLOWMath
using ImplicitAD

println("Running Tests...")
include("compare_objects.jl")

# Active Development
include("post_processing_tests.jl")

# # WORKING:
include("body_aero_tests.jl")

# # In Development
include("body_geometry_tests.jl")
# # TODO: need to add source and body induced velocities to induced velocity on rotor test
# include("rotor_aero_tests.jl")
# # TODO: need to update most of this test after having added sources and bodies in.
# include("dimension_tests.jl")
# # TODO: need tofigure out what to do for testing derivative
# include("derivative_checks.jl")

# # BROKEN:
# include("basic_singularity.jl")
# include("constant_initialization.jl")
# include("wake_aero_tests.jl")

# # include("aero_coefficient_tests.jl")
# # include("basic_geometry.jl")
# # include("geometry_tests.jl")
# # include("rotor_post_processing_tests.jl")
