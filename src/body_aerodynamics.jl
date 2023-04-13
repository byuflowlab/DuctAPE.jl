"""
    calculate_body_strengths_residual(gamb, A_bb, b_bf, vx_bw, vr_bw, wake_gamma)

Calculate body vortex strengths
kid is kutta indices, where kid[1] is the row/column to keep and kid[2] is the row/column to subtract from kid[1] and then delete (1 is the first duct panel, and 2 is the Nth duct panel)
"""
function calculate_body_vortex_strengths!(gamb, A_bb, b_bf, kidx, A_bw, wake_gamma, A_br, sigr)

    # problem dimensions
    nwake = length(A_bw) # number of wake sheets
    nrotor = length(A_br) # number of rotors

    # add freestream contributions to right hand side
    # note: the negative was already included in the precomputation for the freestream.
    b = similar(gamb) .= b_bf

    # add wake vortex sheet contributions to right hand side
    # note: the subtraction was not included in the other coefficients in the precomputation, so we need to subract
    for jwake in 1:nwake

        # get induced velocity in the x-direction
        b .-= A_bw[jwake] * view(wake_gamma, jwake, :)
    end

    # add rotor source sheet contributions to right hand side
    # note: the subtraction is not included in the coefficients here, so we need to subract
    for jrotor in 1:nrotor

        # get induced velocity in the x-direction
        b .-= A_br[jrotor] * view(sigr, : , jrotor)
    end

    #return the residual, rather than solving the linear system.
    #TODO: in the solve functions, need to make sure that we don't subtract things again.
    # return A * gamb - b

    # TODO: do these later after things are working if you need to get things faster.
    # TODO: decide if linear system needs to be solved or not
    # TODO: precompute factorization, test if implicit_linear is faster
    # return ldiv!(gamb, factorize(A_bb), b)

    view(gamb, :) .= solve_body_system(A_bb, b, kidx)

    return nothing
end

"""
Apply Kutta Condition and Solve Linear System
"""
function solve_body_system(A, b, kidx)

    # - Apply the Subtractive Kutta Condition - #
    for i in 1:length(kidx[:, 1])
        # subtract "Nth" row from "1st" row
        A[kidx[i, 1], :] .-= A[kidx[i, 2], :]
        # subtract "Nth" column from "1st" column
        A[:, kidx[i, 1]] .-= A[:, kidx[i, 2]]
        # subtract "Nth" row from "1st" row
        b[kidx[i, 1]] -= b[kidx[i, 2]]
    end

    x = ImplicitAD.implicit_linear(
        A[1:end .∉ [kidx[:, 2]], 1:end .∉ [kidx[:, 2]]], b[1:end .∉ [kidx[:, 2]]]
    )

    # recover the final body panel strength value
    for i in 1:length(kidx[:, 1])
        insert!(x, kidx[i, 2], -x[kidx[i, 1]])
    end

    return x
end
