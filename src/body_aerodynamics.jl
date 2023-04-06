"""
    calculate_body_strengths_residual(gamb, A_bb, b_bf, vx_bw, vr_bw, gamw)

Calculate body vortex strengths
kid is kutta indices, where kid[1] is the row/column to keep and kid[2] is the row/column to subtract from kid[1] and then delete (1 is the first duct panel, and 2 is the Nth duct panel)
"""
function calculate_body_strengths_residual(
    gamb, A_bb, b_bf, kid, vx_bw, vr_bw, gamw, vx_br, vr_br, sigr
)

    # TODO: need to have the A_bb and b_bf matrix and vector be the full sized one, without applying the kutta condition.  then we can apply the kutta condition here after updating the boundary conditions.

    # problem dimensions
    nbody = length(gamb) # number of body panels
    nwake = length(vx_bw) # number of wakes (divided into streamlines)
    nrotor = length(vx_br) # number of rotors

    # add freestream contributions to right hand side
    b = similar(gamb) .= b_bf

    # add wake vortex sheet contributions to right hand side
    for iwake in 1:nwake

        # get induced velocity in the x-direction
        Vx = vx_bw[iwake] * view(gamw, :, iwake)

        # get induced velocity in the r-direction
        Vr = vr_bw[iwake] * view(gamw, :, iwake)

        # add tangential component of velocity to right hand side
        for ibody in 1:nbody
            sa, ca = sincos(body_panels.panel_angle[ibody])
            Vs = Vx[ibody] * ca + Vr[ibody] * sa
            b[ibody] -= Vs
        end
    end

    # # add rotor source sheet contributions to right hand side
    # for irotor in 1:nrotor

    #     # get induced velocity in the x-direction
    #     Vx = vx_br[irotor] * view(sigr, :, irotor)

    #     # get induced velocity in the r-direction
    #     Vr = vr_br[irotor] * view(sigr, :, irotor)

    #     # add normal component of velocity to right hand side
    #     for ibody in 1:nbody
    #         sa, ca = sincos(body_panels[ibody].panel_angle)
    #         Vs = Vx[ibody] * ca + Vr[ibody] * sa
    #         b[ibody] -= Vs
    #     end
    # end

    # - Apply the Subtractive Kutta Condition - #

    #return the residual, rather than solving the linear system.
    #TODO: in the solve functions, need to make sure that we don't subtract things again.
    # return A * gamb - b

    # TODO: do these later after things are working if you need to get things faster.
    # TODO: decide if linear system needs to be solved or not
    # TODO: precompute factorization, test if implicit_linear is faster
    # return ldiv!(gamb, factorize(A_bb), b)

    return solve_body_system(A, b)
end

"""
generates body only linear system
returns A matrix and b vector
"""
function generate_body_linear_system(
    body_panels,
    body_mesh;
    method=ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true]),
)
    system = ff.generate_inviscid_system(method, body_panels, body_mesh)

    return system.A, system.b
end

"""
"""
function solve_body_system(A, b)
    return ImplicitAD.implicit_linear(A, b)
end
