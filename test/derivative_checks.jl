using FiniteDiff
using ForwardDiff
using ReverseDiff
using FLOWFoil
const ff = FLOWFoil

######################################################################
#                                                                    #
#                          FLOWFoil Checks                           #
#                                                                    #
######################################################################

#---------------------------------#
#             Geometry            #
#---------------------------------#

include("data/bodyofrevolutioncoords.jl")
xhub = x
rhub = r

include("data/naca_662-015.jl")
xduct = x
rduct = r

#---------------------------------#
#            Function             #
#---------------------------------#

function run_flow_foil(rductin)
    coordinates = ([xduct rductin], [xhub rhub])

    method = ff.AxisymmetricProblem(Vortex(Constant()), Neumann(), [false, true])
    flow_angle = [0.0]
    reynolds = [-1.0]
    mach = [-1.0]

    # Generate Problem Object
    problem = ff.define_problem(method, coordinates, flow_angle, reynolds, mach)

    # Generate Panel Geometry
    panels = ff.generate_panels(method, coordinates)

    # Generate Influence Mesh
    mesh = ff.generate_mesh(method, panels)

    # Assemble Linear System
    system = ff.generate_inviscid_system(method, panels, mesh)

    # Solve Linear System
    solution = ff.solve(system)

    return solution.x
end

#---------------------------------#
#        FiniteDiff Jacobian      #
#---------------------------------#
findiff_j = FiniteDiff.finite_difference_jacobian(run_flow_foil, rduct)

#---------------------------------#
#       ForwardDiff Jacobian      #
#---------------------------------#
fordiff_j = ForwardDiff.jacobian(run_flow_foil, rduct)

#---------------------------------#
#             Compare             #
#---------------------------------#
println("FLOWFoil Only:")
println("FiniteDiff ≈≈ ForwardDiff: ", isapprox(findiff_j, fordiff_j; atol=1e-3))

######################################################################
#                                                                    #
#                         DuctTAPE Example                           #
#                                                                    #
######################################################################

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
println("DuctTAPE Current Status:")
println("FiniteDiff ≈≈ ForwardDiff: ", isapprox(findiff_j, fordiff_j; atol=1e-3))
