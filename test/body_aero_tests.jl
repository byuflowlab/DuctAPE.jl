@testset "Body Aerodynamic Tests" begin
    body_vortex_strengths = ones(3)
    A_body_to_body = [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    body_boundary_conditions = ones(4)
    wake_gammas = ones(2, 3)
    A_wake_to_body = [[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0] for i in 1:2]

    dt.calculate_body_vortex_strengths!(
        body_vortex_strengths,
        A_body_to_body,
        body_boundary_conditions,
        wake_gammas,
        A_wake_to_body,
    )#, Sigmas, A_rotor_to_body)

    @test body_vortex_strengths == -ones(3)
end
