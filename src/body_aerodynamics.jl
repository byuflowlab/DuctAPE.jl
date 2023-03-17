
"""
    calculate_body_vortex_strengths!(Γb, A_bb, b_fb, Ax_wb, Ar_wb, Γw)

Calculate body vortex strengths
"""
function calculate_body_vortex_strengths!(Γb, A_bb, b_fb, Ax_wb, Ar_wb, Γw, Ax_rb, Ar_rb, Σr)

    # problem dimensions
    nb = length(Γb) # number of body panels
    nw = length(Ax_wb) # number of wakes
    nr = length(Ax_rb) # number of rotors

    # add freestream contributions to right hand side (b_fb = -V∞ ̂n)
    b = similar(Γb) .= b_fb

    # add wake vortex sheet contributions to right hand side
    for iw in 1:nw

        # get induced velocity in the x-direction
        Vx = Ax_wb[iw] * view(Γw, :, iw)

        # get induced velocity in the r-direction
        Vr = Ar_wb[iw] * view(Γw, :, iw)

        # add normal component of velocity to right hand side
        for ib in 1:nb
            sa, ca = sincos(body_panels[ib].panel_angle)
            Vn = Vx[ib]*ca + Vr[ib]*sa
            b[ib] -= Vn
        end

    end

    # add rotor source sheet contributions to right hand side
    for ir in 1:nr

        # get induced velocity in the x-direction
        Vx = Ax_rb[ir] * view(Σr, :, ir)

        # get induced velocity in the r-direction
        Vr = Ar_rb[ir] * view(Σr, :, ir)

        # add normal component of velocity to right hand side
        for ib in 1:nb
            sa, ca = sincos(body_panels[ib].panel_angle)
            Vn = Vx[ib]*ca + Vr[ib]*sa
            b[ib] -= Vn
        end

    end

    # TODO: precompute factorization, test if implicit_linear is fast

    return ldiv!(Γb, factorize(A_bb), b)
end
