######################################################################
#                                                                    #
#                     Linear Panel Integration                       #
#                                                                    #
######################################################################

#---------------------------------#
#             Nominal             #
#---------------------------------#

"""

**Arguments:**
- `t::Float` : sample point in range (0,1) selected by quadrature.
"""
function nominal_vortex_induced_velocity_sample(
    t, node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
)

    # Transform from (0,1) to actual position on panel
    z = linear_transform(nondimrange, [node1[1]; node2[1]], t)
    r = linear_transform(nondimrange, [node1[2]; node2[2]], t)

    # get relative geometry
    xi, rho, m, r0 = calculate_xrm(controlpoint, [z; r])

    # Get velocity components at sample points
    vz = vortex_ring_vz(xi, rho, m, r0, 1.0) #shouldn't need influence length
    vr = vortex_ring_vr(xi, rho, m, r0)

    #=
    assemble output components in the format:
        [x_j r_j; x_{j+1} r_{j+1}]
    and scale by influence panel length
        (due to transformation of integration range to/from range=(0,1))
     =#
    return [vz * [1.0 - t; t] vr * [1.0 - t, t]] * influence_length
end

"""
`V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: [vz_j vr_j; vz_{j+1} vr_{j+1}]
"""
function nominal_vortex_panel_integration(
    node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0], debug=false
)

    # Define function to integrate
    function fsample(t)
        return nominal_vortex_induced_velocity_sample(
            t, node1, node2, influence_length, controlpoint; nondimrange=nondimrange
        )
    end

    V, err = quadgk(fsample, 0.0, 1.0)

    if debug
        return V, err
    else
        return V
    end
end

#---------------------------------#
#             Singular            #
#---------------------------------#

"""
"""
function subtracted_singular_vortex_influence(node, controlpoint)
    rmag2 = (controlpoint[1] - node[1])^2 + (controlpoint[2] - node[2])^2
    num1z = controlpoint[1] - node[1]
    num1r = node[2] - controlpoint[2]
    den1 = 2.0 * pi * rmag2
    den2 = 64.0 * node[2]^2

    axial = num1r / den1 - log(rmag2 / den2) / (8.0 * pi * node[2])
    radial = num1z / den1

    return axial, radial
end

"""
"""
function analytically_integrated_vortex_influence(r, influence_length)
    #mine
    # axial = (influence_length / (2.0 * pi * r)) * (1.0 + log(8.0 * r / influence_length))
    #DFDC
    axial = (1 / (2.0 * pi * r)) * (1.0 + log(16.0 * r / influence_length))
    radial = 0.0
    return axial, radial
end

"""

**Arguments:**
- `t::Float` : sample point in range (0,1) selected by quadrature.
"""
function self_vortex_induced_velocity_sample(
    t, node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
)

    # Transform from (0,1) to actual position on panel
    z0 = linear_transform(nondimrange, [node1[1]; node2[1]], t)
    r0 = linear_transform(nondimrange, [node1[2]; node2[2]], t)

    # get relative geometry
    xi, rho, m, _ = calculate_xrm(controlpoint, [z0; r0])
    # println("xi: ", xi)
    # println("rho: ", rho)
    # println("m: ", m)
    # println("r0: ", r0)

    # Get velocity components at sample points
    vz = vortex_ring_vz(xi, rho, m, r0, 1.0) #shouldn't need influence length
    vr = vortex_ring_vr(xi, rho, m, r0)
    # println("vz: ", vz)
    # println("vr: ", vr)

    # Get singular piece to subtract
    vzs, vrs = subtracted_singular_vortex_influence([z0; r0], controlpoint)
    # println("vzs: ", vzs)
    # println("vrs: ", vrs)

    # Get analytic piece to add back in
    vza, vra = analytically_integrated_vortex_influence(controlpoint[2], influence_length)
    # println("vza: ", vza)
    # println("vra: ", vra)

    #=
    assemble output components in the format:
        [x_j r_j; x_{j+1} r_{j+1}]
    and scale by influence panel length
        (due to transformation of integration range to/from range=(0,1))
     =#
    return [(vz * [1.0 - t; t] .- vzs / 2.0) (vr * [1.0 - t, t] .- vrs / 2.0)]
end

"""
`V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: [vz_j vr_j; vz_{j+1} vr_{j+1}]
"""
function self_vortex_panel_integration(
    node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0], debug=false
)

    # Define function to integrate
    function fsample(t)
        return self_vortex_induced_velocity_sample(
            t, node1, node2, influence_length, controlpoint; nondimrange=nondimrange
        )
    end

    V, err = quadgk(fsample, 0.0, 0.5, 1.0; atol=1e-6)

    vza, vra = analytically_integrated_vortex_influence(controlpoint[2], influence_length)

    V .*= influence_length
    V[:, 1] .+= vza / (2.0 * influence_length)

    if debug
        return V, err
    else
        return V
    end
end