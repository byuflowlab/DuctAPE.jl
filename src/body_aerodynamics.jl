"""
    calculate_body_strengths_residual(gamb, A_bb, b_bf, vx_bw, vr_bw, wake_gamma)

Calculate body vortex strengths
kid is kutta indices, where kid[1] is the row/column to keep and kid[2] is the row/column to subtract from kid[1] and then delete (1 is the first duct panel, and 2 is the Nth duct panel)
"""
function calculate_body_vortex_strengths!(
    gamb, A_bb, b_bf, kidx, A_bw, wake_gamma, A_br, sigr; debug=false
)

    # problem dimensions
    nwake = length(A_bw) # number of wake sheets
    nrotor = length(A_br) # number of rotors

    # add freestream contributions to right hand side
    # note: the negative was already included in the precomputation for the freestream.
    b = similar(gamb) .= b_bf

    if debug
        #initialize extra outputs
        bfree = copy(b_bf)
        bwake = similar(bfree) .= 0.0
        brotor = similar(bfree) .= 0.0
    end

    # add wake vortex sheet contributions to right hand side
    # note: the subtraction was not included in the other coefficients in the precomputation, so we need to subract
    for jwake in 1:nwake

        # get induced velocity in the x-direction
        b .-= A_bw[jwake] * view(wake_gamma, jwake, :)
        if debug
            bwake .-= A_bw[jwake] * view(wake_gamma, jwake, :)
        end
    end

    # add rotor source sheet contributions to right hand side
    # note: the subtraction is not included in the coefficients here, so we need to subract
    for jrotor in 1:nrotor

        # get induced velocity in the x-direction
        b .-= A_br[jrotor] * view(sigr, :, jrotor)
        if debug
            brotor .-= A_br[jrotor] * view(sigr, :, jrotor)
        end
    end

    #return the residual, rather than solving the linear system.
    #TODO: in the solve functions, need to make sure that we don't subtract things again.
    # return A * gamb - b

    # TODO: do these later after things are working if you need to get things faster.
    # TODO: decide if linear system needs to be solved or not
    # TODO: precompute factorization, test if implicit_linear is faster
    # return ldiv!(gamb, factorize(A_bb), b)

    view(gamb, :) .= solve_body_system(A_bb, b, kidx)
    # gamb[:] .= solve_body_system(A_bb, b, kidx)

    if debug
        return bfree, bwake, brotor
    else
        return nothing
    end
end

"""
Apply Kutta Condition and Solve Linear System
"""
function solve_body_system(A_bb, b, kidx)

    # - Apply the Subtractive Kutta Condition - #
    for i in 1:length(kidx[:, 1])
        # subtract "Nth" row from "1st" row
        A_bb[kidx[i, 1], :] .-= A_bb[kidx[i, 2], :]
        # subtract "Nth" column from "1st" column
        A_bb[:, kidx[i, 1]] .-= A_bb[:, kidx[i, 2]]
        # subtract "Nth" row from "1st" row
        b[kidx[i, 1]] -= b[kidx[i, 2]]
    end

    # Solve with subtractive Kutta Condition applied
    # x = ImplicitAD.implicit_linear(
    #     A_bb[1:end .∉ [kidx[:, 2]], 1:end .∉ [kidx[:, 2]]], b[1:end .∉ [kidx[:, 2]]]
    # )

    x = A_bb[1:end .∉ [kidx[:, 2]], 1:end .∉ [kidx[:, 2]]] \ b[1:end .∉ [kidx[:, 2]]]

    # - put things back after solving - #
    for i in 1:length(kidx[:, 1])
        # add "Nth" row back to "1st" row
        A_bb[kidx[i, 1], :] .+= A_bb[kidx[i, 2], :]
        # add "Nth" column back to "1st" column
        A_bb[:, kidx[i, 1]] .+= A_bb[:, kidx[i, 2]]
        # add "Nth" row back to "1st" row
        b[kidx[i, 1]] += b[kidx[i, 2]]
    end

    # recover the final body panel strength value
    for i in 1:length(kidx[:, 1])
        insert!(x, kidx[i, 2], -x[kidx[i, 1]])
    end

    return x
end
