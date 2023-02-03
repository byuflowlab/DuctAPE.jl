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

include("rotor_aero_tests.jl")
include("wake_aero_tests.jl")
include("body_aero_tests.jl")

