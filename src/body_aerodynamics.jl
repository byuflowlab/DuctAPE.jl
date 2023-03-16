
"""
    calculate_body_vortex_strengths!(Γ_b, A_bb, b0, Γ_w, Ax_wb, Ar_wb)

Calculate body vortex strengths
"""
function calculate_body_vortex_strengths!(Γ_b, A_bb, b0, Γ_w, Ax_wb, Ar_wb, Σ_r, Ax_rb, Ar_rb)

    # number of body panels
    nb = length(Γ_b)
    nw = length(Ax_wb)
    nr = length(Ax_rb)

    # initialize right hand side (b0 = -V∞ ̂n)
    b = similar(Γ_b) .= b0

    # add wake contributions to right hand side
    for iw in 1:nw

        # get induced velocity in the x-direction
        Vx = Ax_wb[iw] * view(Γ_w, iw,:)

        # get induced velocity in the r-direction
        Vr = Ar_wb[iw] * view(Γ_w, iw,:)

        # add normal component of velocity to right hand side
        for ib in 1:nb
            sa, ca = sincos(body_panels.panel_angle)
            Vn = Vx[ib]*ca + Vr[ib]*sa
            b[ib] -= Vn
        end

    end

    # add rotor contributions to right hand side
    for ir in 1:nr

        # get induced velocity in the x-direction
        Vx = Ax_rb[ir] .* view(Σ_r, i, :)

        # get induced velocity in the r-direction
        Vr = Ar_rb[ir] .* view(Σ_r, i, :)

        # add normal component of velocity to right hand side
        for ib in 1:nb
            sa, ca = sincos(body_panels.panel_angle)
            Vn = Vx[ib]*ca + Vr[ib]*sa
            b[ib] -= Vn
        end

    end

    # calculate and return solution
    return Γ_b .= A_bb\b
end
