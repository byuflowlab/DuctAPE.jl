@testset "Airfoil Caching Tests" begin

    # - Test AlphaAF - #
    alphas = collect(range(-10.0, 10.0, 21))
    cl = collect(range(0.0, 20.0, 21))
    cd = collect(range(0.0, 0.2, 21))
    af = dt.c4b.AlphaAF(alphas, cl, cd)
    total_length = [0]
    cache = dt.allocate_airfoil_cache([[af]], total_length, 1, 1)
    vec = zeros(total_length[])
    afs = dt.withdraw_airfoil_cache(vec, cache)
    afs[1] = af

    # - Test AlphaReAF - #
    alphas = collect(range(-10.0, 10.0, 21))
    res = [1, 2]
    cl = reshape(collect(range(0.0, 40.0, 21 * 2)), (21, 2))
    cd = reshape(collect(range(0.0, 0.4, 21 * 2)), (21, 2))
    af = dt.c4b.AlphaReAF(alphas, res, cl, cd)
    total_length = [0]
    dims, cons = dt.allocate_airfoil_cache([[af]], total_length, 1, 1)

    # - Test AlphaMachAF - #
    alphas = collect(range(-10.0, 10.0, 21))
    machs = [1, 2]
    cl = reshape(collect(range(0.0, 40.0, 21 * 2)), (21, 2))
    cd = reshape(collect(range(0.0, 0.4, 21 * 2)), (21, 2))
    af = dt.c4b.AlphaMachAF(alphas, machs, cl, cd)
    total_length = [0]
    dims, cons = dt.allocate_airfoil_cache([[af]], total_length, 1, 1)

    # - Test AlphaReMachAF - #
    alphas = collect(range(-10.0, 10.0, 21))
    res = [1, 2]
    machs = [3, 4]
    cl = reshape(collect(range(0.0, 80.0, 21 * 4)), (21, 2, 2))
    cd = reshape(collect(range(0.0, 0.8, 21 * 4)), (21, 2, 2))
    af = dt.c4b.AlphaReMachAF(alphas, res, machs, cl, cd)
    total_length = [0]
    dims, cons = dt.allocate_airfoil_cache([[af]], total_length, 1, 1)
end
