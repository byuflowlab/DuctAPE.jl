println("\nBODY INFLUENCE COEFFICIENT TESTS")
#---------------------------------#
# Body System Influence Matrices  #
#---------------------------------#
@testset "Body System Influence Matrices" begin

    #TODO: what kind of test would be useful here?
end

#---------------------------------#
#   Off-Body Influence Matrices   #
#---------------------------------#
##TODO: below here are the off-body cases

#@testset "Vortex Induced Velocities" begin
#    # setup
#    node = [1.0; 1.0]
#    controlpoint = [0.0; 2.0]
#    influence_length = 1.0
#    xi, rho, m, r_influence = dt.calculate_xrm(controlpoint, node)

#    # individual components
#    vz = dt.vortex_ring_vz(xi, rho, m, r_influence, influence_length)
#    vr = dt.vortex_ring_vr(xi, rho, m, r_influence)

#    # components together
#    vel1 = dt.vortex_induced_velocity(controlpoint, node, influence_length)
#    @test vel1 == [vz; vr]

#    # components in place
#    vel2 = similar(vel1) .= 0.0
#    dt.vortex_induced_velocity!(vel2, controlpoint, node, influence_length)

#    @test vel1 == vel2
#end

#@testset "Vortex Coefficient Functions" begin
#    controlpoints = [1.0 1.0; 2.0 1.0]
#    nodes = [1.0 1.0; 1.0 2.0]
#    influence_lengths = 2.0 * ones(2)
#    strengths = 2.0 * ones(2)

#    #TODO: set up control points(s) and node(s), and get individual values to test this first one, then the others are all relative to that.
#    v = zeros(2, 2, 2)

#    for (ic, cp) in enumerate(eachrow(controlpoints))
#        for (in, np) in enumerate(eachrow(nodes))
#            dt.vortex_induced_velocity!(
#                view(v, ic, in, :), cp, np, influence_lengths[in], strengths[in]
#            )
#        end
#    end

#    AICcomp1 = dt.influencefromvortices(controlpoints, nodes, influence_lengths, strengths)
#    @test AICcomp1[:, :, 1] == v[:, :, 1]
#    @test AICcomp1[:, :, 2] == v[:, :, 2]

#    AICcomp2 = similar(AICcomp1) .= 0.0
#    dt.influencefromvortices!(AICcomp2, controlpoints, nodes, influence_lengths, strengths)
#    @test AICcomp1 == AICcomp2

#    V1 = dt.vfromvortices(controlpoints, nodes, influence_lengths, strengths)
#    @test reduce(vcat, V1) == reduce(vcat, sum(AICcomp2; dims=2))

#    V2 = similar(V1) .= 0.0
#    dt.vfromvortices!(V2, controlpoints, nodes, influence_lengths, strengths)
#    @test V2 == V1
#end
