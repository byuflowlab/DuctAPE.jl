using FiniteDiff
using ForwardDiff
using ReverseDiff
using FLOWFoil
const ff = FLOWFoil

######################################################################
#                                                                    #
#                         DuctTAPE Example                           #
#                                                                    #
######################################################################

println("DuctTAPE Current Status:")
include("../examples/lilium_ish.jl")

#---------------------------------#
#        FiniteDiff Jacobian      #
#---------------------------------#
findiff_j = FiniteDiff.finite_difference_jacobian(wrapper, x0)

#---------------------------------#
#       ForwardDiff Jacobian      #
#---------------------------------#
fordiff_j = ForwardDiff.jacobian(wrapper, x0)

#---------------------------------#
#             Compare             #
#---------------------------------#
println("FiniteDiff ≈≈ ForwardDiff: ", isapprox(findiff_j, fordiff_j; atol=1e-3))
println("Maximum Diff: ", maximum(findiff_j .- fordiff_j))
