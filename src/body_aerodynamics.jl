"""
    calculate_body_vortex_strengths!(Γb, A_bb, b_bf, vx_bw, vr_bw, Γw)

Calculate body vortex strengths
"""
function calculate_body_vortex_strengths!(
    Γb, A_bb, b_bf, vx_bw, vr_bw, Γw, vx_br, vr_br, Σr
)

#TODO: need to consider how the subtractive Kutta conditions affects things here.
#TODO: build full system then apply subtractive stuff? flowfoil already did the subtractive stuff for the body only case...
#TODO: perhaps need to move axisymmetric content from FLOWFoil over to DuctTAPE completely...
#TODO, or separate out function in flowfoil to make original matrix, then add a separate funciton to add back diagonal correction and another to apply kutta condition to system?

    # problem dimensions
    nbody = length(Γb) # number of body panels
    nwake = length(vx_bw) # number of wakes (divided into streamlines)
    nrotor = length(vx_br) # number of rotors

    # add freestream contributions to right hand side
    b = similar(Γb) .= b_bf

    # add wake vortex sheet contributions to right hand side
    for iwake in 1:nwake

        # get induced velocity in the x-direction
        Vx = vx_bw[iwake] * view(Γw, :, iwake)

        # get induced velocity in the r-direction
        Vr = vr_bw[iwake] * view(Γw, :, iwake)

        # add tangential component of velocity to right hand side
        for ibody in 1:nbody
            sa, ca = sincos(body_panels.panel_angle[ibody])
            Vt = Vx[ibody] * ca + Vr[ibody] * sa
            b[ibody] -= Vt
        end
    end

    # add rotor source sheet contributions to right hand side
    for irotor in 1:nrotor

        # get induced velocity in the x-direction
        Vx = vx_br[irotor] * view(Σr, :, irotor)

        # get induced velocity in the r-direction
        Vr = vr_br[irotor] * view(Σr, :, irotor)

        # add normal component of velocity to right hand side
        for ibody in 1:nbody
            sa, ca = sincos(body_panels[ibody].panel_angle)
            Vt = Vx[ibody] * ca + Vr[ibody] * sa
            b[ibody] -= Vt
        end
    end

    #return the residual, rather than solving the linear system.
    #TODO: in the solve functions, need to make sure that we don't subtract things again.
    return A*Γb - b

    # TODO: precompute factorization, test if implicit_linear is faster
    # return ldiv!(Γb, factorize(A_bb), b)
end
