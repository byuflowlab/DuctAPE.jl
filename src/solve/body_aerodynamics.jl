"""
    calculate_body_strengths_residual(gamb, A_bb, b_bf, vz_bw, vr_bw, gamw)

Calculate body vortex strengths
kid is kutta indices, where kid[1] is the row/column to keep and kid[2] is the row/column to subtract from kid[1] and then delete (1 is the first duct panel, and 2 is the Nth duct panel)
"""
function calculate_body_vortex_strengths!(
    gamb, A_bb, b_bf, gamw, A_bw, A_pw, sigr, A_br, A_pr, RHS; post=false
)

    # problem dimensions
    nbn = size(A_bw, 1) # number of body nodes
    nrotor = size(A_br, 3) # number of rotors

    # add freestream contributions to right hand side
    RHS .= b_bf

    if post
        RHSw = similar(RHS) .= 0
        RHSr = similar(RHS) .= 0
    end

    # add wake vortex sheet contributions to right hand side
    RHS[1:nbn] .-= A_bw * gamw
    # add wake influence on psuedo control points
    # TODO: may want to make this not hard coded at some point
    RHS[nbn + 1] -= A_pw[1, :]' * gamw
    RHS[nbn + 4] -= A_pw[2, :]' * gamw

    if post
        RHSw[1:nbn] .-= A_bw * gamw
        RHSw[nbn + 1] -= A_pw[1, :]' * gamw
        RHSw[nbn + 4] -= A_pw[2, :]' * gamw
    end

    # add rotor source sheet contributions to right hand side
    for jrotor in 1:nrotor
        # get induced velocity in the x-direction
        @views RHS[1:nbn] .-= A_br[:, :, jrotor] * sigr[:, jrotor]
        # add rotor influence on pseudo control points
        # TODO: may want to make this not hard coded at some point
        @views RHS[nbn + 1] -= A_pr[1, :, jrotor]' * sigr[:, jrotor]
        @views RHS[nbn + 4] -= A_pr[2, :, jrotor]' * sigr[:, jrotor]

        if post
            @views RHSr[1:nbn] .-= A_br[:, :, jrotor] * sigr[:, jrotor]
            @views RHSr[nbn + 1] -= A_pr[1, :, jrotor]' * sigr[:, jrotor]
            @views RHSr[nbn + 4] -= A_pr[2, :, jrotor]' * sigr[:, jrotor]
        end
    end

    if post # return gamb, wake RHS, and first rotor RHS
        return ldiv!(gamb, A_bb, RHS), RHSw, RHSr
    else
        return ldiv!(gamb, A_bb, RHS)
        # gamb .= A_bb \ RHS
        # return gamb
    end
end
