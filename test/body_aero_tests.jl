@testset "Body Aerodynamic Tests" begin
    gamb = ones(4)
    A_bb = [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    b_bf = ones(4)
    gamw = ones(2, 4)
    A_bw = [
        [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0] for i in 1:2
    ]
    sigr = ones(4,2)
    A_br = [
        [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0] for i in 1:2
    ]
    kidx = [1 4]

    dt.calculate_body_vortex_strengths!(gamb, A_bb, b_bf, kidx, A_bw, gamw, A_br, sigr)

    A = [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    A[1, :] .-= A[end, :]
    A[:, 1] .-= A[:, end]


    # 1.0 minus 2.0 for wakes minus 2 for rotors = -3
    b = -3.0*ones(4)
    b[1] -= b[end]

    x = A[1:3,1:3]\b[1:3]

    @test [x;-x[1]] == gamb

end
