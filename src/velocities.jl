#---------------------------------#
#    Unit Induced Velocities      #
#---------------------------------#
"""
"""
function vortex_band_induced_velocity(xi, rho, m, r_influence, ds_influence)

    # get unit induced velocities of vortex ring
    vx = vortex_ring_vx(xi, rho, m, r_influence)
    vr = vortex_ring_vr(xi, rho, m, r_influence)

    # multiply by ds to get band unit induced velocity
    return [vx; vr].*ds_influence
end

"""
"""
function vortex_ring_vx(xi, rho, m, r_influence)

    #get the first denominator
    den1 = 2.0 * pi * r_influence * sqrt(xi^2 + (rho + 1.0)^2)

    #get numerator and denominator of second fraction
    num2 = 2 * (rho - 1)
    den2 = xi^2 + (rho - 1)^2

    #get values for elliptic integrals
    K, E = get_elliptics(m)

    # set self-induced case to zero for the off-body case
    if xi^2 + (rho - 1.0)^2 <= eps()
        return 0.0
    else
        # negative here is due to our convention that the vortex is postive clockwise
        return -1.0 / den1 * (K - (1.0 + num2 / den2) * E)
    end

end

"""
"""
function vortex_ring_vr(xi, rho, m, r_influence)

    #get numerator and denominator of first fraction
    num1 = xi / rho
    den1 = 2.0 * pi * r_influence * sqrt(xi^2 + (rho + 1.0)^2)

    num2 = 2 * rho
    den2 = xi^2 + (rho - 1)^2

    #get values for elliptic integrals
    K, E = get_elliptics(m)

    # set self-induced case to zero for the off-body case
    if xi^2 + (rho - 1.0)^2 <= eps()
        return 0.0
    else
        return num1 / den1 * (K - (1.0 + num2 / den2) * E)
    end
end

