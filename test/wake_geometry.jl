@testset "Fill in Wake Strengths" begin
    gamw = [1.0; 2.0; 3.0]
    rotor_indices = [1]
    num_wake_x_panels = 5

    wake_vortex_strengths = DuctAPE.fill_out_wake_strengths(gamw, rotor_indices, num_wake_x_panels)

    @test wake_vortex_strengths == repeat(gamw; inner=(1, num_wake_x_panels))

    gamw = [1.0 10.0; 2.0 20.0; 3.0 30.0]
    rotor_indices = [1; 5]
    num_wake_x_panels = 10

    wake_vortex_strengths = DuctAPE.fill_out_wake_strengths(gamw, rotor_indices, num_wake_x_panels)

    @test wake_vortex_strengths[:,1:4] == repeat(gamw[:,1]; inner=(1, 4))
    @test wake_vortex_strengths[:,5:end] == repeat(gamw[:,2]; inner=(1, 6))

end
