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

#TODO: fix all these after updates
# include("body_aero_tests.jl")
# include("geometry_tests.jl")
# include("aero_coefficient_tests.jl")
# include("rotor_aero_tests.jl")
# include("wake_aero_tests.jl")
# include("dimension_tests.jl")

include("basic_singularity.jl")
