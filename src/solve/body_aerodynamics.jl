"""
    calculate_body_strengths_residual(gamb, A_bb, b_bf, vz_bw, vr_bw, gamw)

Calculate body vortex strengths
kid is kutta indices, where kid[1] is the row/column to keep and kid[2] is the row/column to subtract from kid[1] and then delete (1 is the first duct panel, and 2 is the Nth duct panel)
"""
function calculate_body_vortex_strengths!(
    gamb, A_bb, b_bf, gamw, A_bw, sigr, A_br, RHS; debug=false
)

    # problem dimensions
    nbn = size(A_bw, 1) # number of body nodes
    nrotor = size(A_br, 3) # number of rotors

    # add freestream contributions to right hand side
    # note: the negative was already included in the precomputation for the freestream.
    RHS .= b_bf

    # add wake vortex sheet contributions to right hand side
    # note: the subtraction was not included in the other coefficients in the precomputation, so we need to subract
    # get induced velocity in the x-direction
    RHS[1:nbn] .-= A_bw * gamw

    # add rotor source sheet contributions to right hand side
    # note: the subtraction is not included in the coefficients here, so we need to subract
    for jrotor in 1:nrotor
        # get induced velocity in the x-direction
        @views RHS[1:nbn] .-= A_br[:, :, jrotor] * sigr[:, jrotor]
    end

    return ldiv!(gamb, A_bb, RHS)
end
