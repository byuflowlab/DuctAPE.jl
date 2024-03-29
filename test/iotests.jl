println("\nI/O TESTS")

@testset "State Variable Vectorization and Reshaping" begin
    v1 = [1 4; 2 5; 3 6]
    v2 = [7 10; 8 11; 9 12]
    v3 = [i for i in 13:26]

    vars, dims = dt.vectorize_velocity_states(v1, v2, v3)

    @test vars == [i for i in 1:26]

    V1, V2, V3 = dt.extract_state_variables(dt.SolverOptions(), vars, dims)

    @test V1 == v1
    @test V2 == v2
    @test V3 == v3

    v1 = [1 4; 2 5; 3 6]
    v2 = [7 10; 8 11; 9 12]
    v3 = [i for i in 13:26]

    vars, dims = dt.vectorize_strength_states(v1, v2, v3)

    @test vars == [i for i in 1:26]

    V1, V2, V3 = dt.extract_state_variables(dt.CSORSolverOptions(), vars, dims)

    @test V1 == v1
    @test V2 == v2
    @test V3 == v3
end
