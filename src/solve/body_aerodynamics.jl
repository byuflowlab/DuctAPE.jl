"""
    calculate_body_strengths_residual(gamb, A_bb_LU, b_bf, vz_bw, vr_bw, gamw)

Calculate body vortex strengths
kid is kutta indices, where kid[1] is the row/column to keep and kid[2] is the row/column to subtract from kid[1] and then delete (1 is the first duct panel, and 2 is the Nth duct panel)
"""
function calculate_body_vortex_strengths!(
    gamb, A_bb_LU, b_bf, gamw, A_bw, A_pw, sigr, A_br, A_pr, A_bb, rhs; post=false
)

    # problem dimensions
    nbn = size(A_bw, 1) # number of body nodes
    nrotor = size(A_br, 3) # number of rotors

    # add freestream contributions to right hand side
    rhs .= b_bf

    if post
        RHSw = similar(RHS) .= 0
        RHSr = similar(RHS) .= 0
    end

    # add wake vortex sheet contributions to right hand side
    rhs[1:nbn] .-= A_bw * gamw

    # add wake influence on psuedo control points
    # TODO: may want to make this not hard coded at some point
    rhs[nbn + 1] -= @view(A_pw[1, :])' * gamw
    rhs[nbn + 4] -= @view(A_pw[2, :])' * gamw

    if post
        RHSw[1:nbn] .-= A_bw * gamw
        RHSw[nbn + 1] -= @view(A_pw[1, :])' * gamw
        RHSw[nbn + 4] -= @view(A_pw[2, :])' * gamw
    end

    # add rotor source sheet contributions to right hand side
    nws = size(sigr, 1)
    for jrotor in 1:nrotor
        rotorrange = (nws * (jrotor - 1) + 1):(nws * jrotor)
        # get induced velocity in the x-direction
        @views rhs[1:nbn] .-= A_br[:, rotorrange] * sigr[:, jrotor]
        # add rotor influence on pseudo control points
        # TODO: may want to make this not hard coded at some point
        @views rhs[nbn + 1] -= A_pr[1, rotorrange]' * sigr[:, jrotor]
        @views rhs[nbn + 4] -= A_pr[2, rotorrange]' * sigr[:, jrotor]

        if post
            @views RHSr[1:nbn] .-= A_br[:, rotorrange] * sigr[:, jrotor]
            @views RHSr[nbn + 1] -= A_pr[1, rotorrange]' * sigr[:, jrotor]
            @views RHSr[nbn + 4] -= A_pr[2, rotorrange]' * sigr[:, jrotor]
        end
    end

    # if post # return gamb, wake RHS, and first rotor RHS
    # return ImplicitAD.implicit_linear(A_bb, gamb; lsolve=ldiv!, Af=A_bb_LU), RHSw, RHSr
    # return ldiv!(gamb, A_bb_LU, RHS), RHSw, RHSr
    # else
        # printval("gamb_just_before_iad= ", gamb)
        # println()

        # use ImplicitAD overwrite gamb
        # return ImplicitAD.implicit_linear(A_bb, gamb; lsolve=ldiv!, Af=A_bb_LU)
        # gamb .= ImplicitAD.implicit_linear(A_bb, rhs; lsolve=ldiv!, Af=A_bb_LU)
        # return gamb

        # use ldiv! in place for gamb
        # return ldiv!(gamb, A_bb_LU, RHS)

        # just left divide and overwrite gamb
        gamb .= A_bb_LU \ rhs
        return gamb

    # end
end
