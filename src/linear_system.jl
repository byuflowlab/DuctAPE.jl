#########################################
#                                       #
#              LHS Matrix               #
#                                       #
#########################################
"""
body on body coefficients
"""
function generate_LHS()

    return aic
end

#########################################
#                                       #
#             RHS Vectors               #
#                                       #
#########################################

"""
Vinf induced boundary condition coefficients
"""
function generate_aicvinf(Vinfx, panels)

    naic = length(panels.panel_normal)
    aicvinf = zeros(TF, naic)

    for ia in 1:naic
        aicvinf[ia] = vdotn(Vinfx, panels.panel_normal[ia])
    end

    return aicvinf
end

"""
Wake induced boundary condition coefficients
"""
function generate_aicwake()
    return aicwake
end

"""
rotor induced boundary condition coefficients
"""
function generate_aicrotor()
    return aicrotor
end

#########################################
#                                       #
#         Auxiliary Functions           #
#                                       #
#########################################

"""
V dot n
"""
function vdotn(vvec, nhat)
    return vvec[1]*nhat[1] + vvec[2]*nhat[2]
end

