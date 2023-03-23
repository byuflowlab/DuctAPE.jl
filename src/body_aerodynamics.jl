
"""
    calculate_body_vortex_strengths!(Γb, A_bb, b_fb, Ax_wb, Ar_wb, Γw)

Calculate body vortex strengths
"""
function calculate_body_vortex_strengths!(Γb, A_bb, b_fb, Ax_wb, Ar_wb, Γw, Ax_rb, Ar_rb, Σr)

    # problem dimensions
    nbody = length(Γb) # number of body panels
    nwake = length(Ax_wb) # number of wakes (divided into streamlines)
    nrotor = length(Ax_rb) # number of rotors

    # add freestream contributions to right hand side
    b = similar(Γb) .= b_fb

    # add wake vortex sheet contributions to right hand side
    for iwake in 1:nwake

        # get induced velocity in the x-direction
        Vx = Ax_wb[iwake] * view(Γw, :, iwake)

        # get induced velocity in the r-direction
        Vr = Ar_wb[iwake] * view(Γw, :, iwake)

        # add tangential component of velocity to right hand side
        for ibody in 1:nbody
            sa, ca = sincos(body_panels.panel_angle[ibody])
            Vt = Vx[ibody]*ca + Vr[ibody]*sa
            b[ibody] -= Vt
        end

    end

    # add rotor source sheet contributions to right hand side
    for irotor in 1:nrotor

        # get induced velocity in the x-direction
        Vx = Ax_rb[irotor] * view(Σr, :, irotor)

        # get induced velocity in the r-direction
        Vr = Ar_rb[irotor] * view(Σr, :, irotor)

        # add normal component of velocity to right hand side
        for ibody in 1:nbody
            sa, ca = sincos(body_panels[ibody].panel_angle)
            Vt = Vx[ibody]*ca + Vr[ibody]*sa
            b[ibody] -= Vt
        end

    end

    # TODO: precompute factorization, test if implicit_linear is faster

    return ldiv!(Γb, factorize(A_bb), b)
end
