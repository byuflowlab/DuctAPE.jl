"""
    calculate_body_vortex_strengths!(
        gamb, A_bb_LU, b_bf, gamw, A_bw, A_pw, sigr, A_br, A_pr, A_bb, rhs; post=false
    ) -> gamb

Computes the bound vortex strengths (`gamb`) on the body panels by assembling and solving
a linear system that accounts for contributions from the freestream, wake vortex sheets,
and rotor source panels.

# Arguments
- `gamb`: Output vector where the solved body vortex strengths will be stored.
- `A_bb_LU`: LU factorization of the body-body influence matrix `A_bb` for efficient solving.
- `b_bf`: Right-hand-side vector contribution from the freestream.
- `gamw`: Wake vortex strengths.
- `A_bw`: Influence matrix of the wake vortex sheet on the body.
- `A_pw`: Influence matrix of the wake vortex sheet on pseudo control points.
- `sigr`: Source strengths on the rotor surfaces.
- `A_br`: Influence matrix of the rotor source sheets on the body.
- `A_pr`: Influence matrix of the rotor source sheets on pseudo control points.
- `A_bb`: The full body-body influence matrix (used for ImplicitAD system tracing).
- `rhs`: Preallocated vector for the assembled right-hand side of the system.
- `post` (keyword, default = `false`): If `true`, stores intermediate residuals for post-processing
  (via `RHSw` and `RHSr`, though they are not returned from this function).

# Returns
- `gamb`: The computed body vortex strengths satisfying the linear system.
"""
function calculate_body_vortex_strengths!(
    gamb, A_bb_LU, b_bf, gamw, A_bw, A_pw, sigr, A_br, A_pr, A_bb, rhs; post=false
)

    # problem dimensions
    nbn = size(A_bw, 1) # number of body nodes
    nrotor = size(A_br, 3) # number of rotors

    # add freestream contributions to right hand side
    rhs .= -b_bf

    if post
        RHSw = similar(RHS) .= 0
        RHSr = similar(RHS) .= 0
    end

    # add wake vortex sheet contributions to right hand side
    rhs[1:nbn] .-= A_bw * gamw

    # add wake influence on psuedo control points
    # TODO: may want to make this not hard coded at some point
    rhs[nbn + 1] -= view(A_pw, 1, :)' * gamw
    rhs[nbn + 4] -= view(A_pw, 2, :)' * gamw

    if post
        RHSw[1:nbn] .-= A_bw * gamw
        RHSw[nbn + 1] -= view(A_pw, 1, :)' * gamw
        RHSw[nbn + 4] -= view(A_pw, 2, :)' * gamw
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

    gamb .= ImplicitAD.implicit_linear(A_bb, rhs; lsolve=ldiv!, Af=A_bb_LU)
    return gamb

    # gamb .= A_bb_LU \ rhs
    # return gamb

end
