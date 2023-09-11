#---------------------------------#
#           LHS Matrix            #
#---------------------------------#
"""
body on body coefficients
"""
function generate_AIC(v_bb, panels)
    N = length(panels.panel_length)

    AIC = zeros(eltype(panels.panel_center), N)

    for r in 1:N # for rows in matrix
        for c in 1:N # for columns in matrix
            AIC[r, c] = vdotn(v_bb[r, c, :], panels.panel_normal[r])
        end # for columns
    end # for rows

    return AIC
end

#---------------------------------#
#          RHS Vectors            #
#---------------------------------#
"""
generate aic vector for unit induced velocities on body
"""
function generate_aic_vec(v_ai, panels)
    naic = length(panels.panel_normal)
    aic = zeros(TF, naic)

    for ia in 1:naic
        aic[ia] = vdotn(v_ai, panels.panel_normal[ia])
    end

    return aic
end

#---------------------------------#
#      Auxiliary Functions        #
#---------------------------------#
"""
V dot n
"""
function vdotn(vvec, nhat)
    return vvec[1] * nhat[1] + vvec[2] * nhat[2]
end

