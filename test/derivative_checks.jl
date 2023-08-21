include("derivative_wrappers.jl")

@testset "Automaic Derivatives" begin
    Vinf = 20.0
    chords = [
        0.089142
        0.079785
        0.0713
        0.063979
        0.057777
        0.052541
        0.048103
        0.044316
        0.041061
        0.038243
    ]
    ductr = [
        0.835
        0.834496
        0.832676
        0.828558
        0.820991
        0.810449
        0.797709
        0.784153
        0.771912
        0.764112
        0.7605
        0.760275
        0.762716
        0.767502
        0.774294
        0.7827
        0.792384
        0.802639
        0.812809
        0.822477
        0.835
        0.847523
        0.857191
        0.867361
        0.877616
        0.8873
        0.895706
        0.902498
        0.907284
        0.909725
        0.9095
        0.905888
        0.898088
        0.885847
        0.872291
        0.859551
        0.849009
        0.841442
        0.837324
        0.835504
        0.835
    ]
    r = collect(0.16:0.06669722222222223:0.760275)

    inputs = [Vinf; chords]
    # outs = @time dt_full_wrapper(inputs)

    # println("\trunning pre/post AD tests")
    # # - Pre and Post computation test Test - #
    # precompinputs = [Vinf; chords]
    # # ForwardDiff Jacobian
    # fordiff_j = ForwardDiff.jacobian(dt_prepost_wrapper, precompinputs)
    # # FiniteDiff Jacobian
    # findiff_j = FiniteDiff.finite_difference_jacobian(dt_prepost_wrapper, precompinputs)
    # # compare, note scaling is an issue here, there are some really large and really small numbers
    # @test isapprox(findiff_j, fordiff_j; atol=1.2e-2)

    # - Full Solver Test - #
    # ForwardDiff Jacobian
    println("\tCalculating ForwardDiff Jacobian")
    # cnfg = ForwardDiff.JacobianConfig(dt_full_wrapper, inputs, ForwardDiff.Chunk{41}())
    # fordiff_j = zeros(4, length(inputs))
    # fordiff_j = @time ForwardDiff.jacobian(dt_full_wrapper, inputs, cnfg)
    fordiff_j = @time ForwardDiff.jacobian(dt_full_wrapper, inputs)

    fo = open("fordiff_j.jl", "w")
    write(fo, "fordiff_j=")
    write(fo, "$(fordiff_j)")
    close(fo)

    # FiniteDiff Jacobian
    println("\tCalculating FiniteDiff Jacobian")
    # findiff_j = similar(fordiff_j) .= 0.0
    # @time FiniteDiff.finite_difference_jacobian!(findiff_j, dt_full_wrapper, inputs)
    findiff_j = @time FiniteDiff.finite_difference_jacobian(dt_full_wrapper, inputs)

    fi = open("findiff_j.jl", "w")
    write(fi, "findiff_j=")
    write(fi, "$(findiff_j)")
    close(fi)

    @test all(isapprox.(findiff_j, fordiff_j; atol=2e-3))
end
