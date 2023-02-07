using Test

println("Compiling Package")
using DuctTAPE
const dt = DuctTAPE
using CCBlade
const ccb = CCBlade
using ForwardDiff
const frd = ForwardDiff
using FiniteDiff
const fnd = FiniteDiff

include("body_aero_tests.jl")
include("geometry_tests.jl")
include("aero_coefficient_tests.jl")
include("rotor_aero_tests.jl")

include("wake_aero_tests.jl")
