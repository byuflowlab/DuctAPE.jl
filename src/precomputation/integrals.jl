######################################################################
#                                                                    #
#                 Linear Vortex Panel Integration                    #
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
    # return StaticArrays.SMatrix{2,2}(
    #     [vz * [1.0 - t; t] vr * [1.0 - t; t]] * influence_length
    # )
    return StaticArrays.SVector{4}(
        [vz * (1.0 - t); vz * t; vr * (1.0 - t); vr * t] * influence_length
    )
    #  return [vz * [1.0 - t; t] vr * [1.0 - t; t]] * influence_length
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

    V, err = quadgk(fsample, 0.0, 1.0; maxevals=1e2)

    if debug
        return reshape(V, (2, 2)), err
        # return V, err
    else
        return reshape(V, (2, 2))
        # return V
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
    den2 = 64.0 * controlpoint[2]^2

    if isapprox(rmag2, 0.0) || isapprox(controlpoint[2], 0.0)
        axial = 0.0
    else
        axial = num1r / den1 - log(rmag2 / den2) / (8.0 * pi * controlpoint[2])
    end

    if isapprox(controlpoint[2], 0.0)
        radial = 0.0
    else
        radial = num1z / den1
    end

    return axial, radial
end

"""
"""
function analytically_integrated_vortex_influence(r, influence_length)
    axial = (influence_length / (4.0 * pi * r)) * (1.0 + log(16.0 * r / influence_length))
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

    # Get velocity components at sample points
    vz = vortex_ring_vz(xi, rho, m, r0, 1.0) #shouldn't need influence length
    vr = vortex_ring_vr(xi, rho, m, r0)

    # Get singular piece to subtract
    vzs, vrs = subtracted_singular_vortex_influence([z0; r0], controlpoint)

    #=
    assemble output components in the format:
        [x_j r_j; x_{j+1} r_{j+1}]
    and scale by influence panel length
        (due to transformation of integration range to/from range=(0,1))
     =#
    return [
        vz * (1.0 - t) - vzs / 2.0
        vz * t .- vzs / 2.0
        vr * (1.0 - t) .- vrs / 2.0
        vr * t .- vrs / 2.0
    ]
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

    V, err = quadgk(fsample, 0.0, 0.5, 1.0; maxevals=1e2)

    vza, _ = analytically_integrated_vortex_influence(controlpoint[2], influence_length)

    V .*= influence_length
    V[1:2] .+= vza / 2.0

    if debug
        return reshape(V, (2, 2)), err
    else
        return reshape(V, (2, 2))
    end
end

######################################################################
#                                                                    #
#                 Linear Source Panel Integration                    #
#                                                                    #
######################################################################

#---------------------------------#
#             Nominal             #
#---------------------------------#

"""

**Arguments:**
- `t::Float` : sample point in range (0,1) selected by quadrature.
"""
function nominal_source_induced_velocity_sample(
    t, node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
)

    # Transform from (0,1) to actual position on panel
    z = linear_transform(nondimrange, [node1[1]; node2[1]], t)
    r = linear_transform(nondimrange, [node1[2]; node2[2]], t)

    # get relative geometry
    xi, rho, m, r0 = calculate_xrm(controlpoint, [z; r])

    # Get velocity components at sample points
    vz = source_ring_vz(xi, rho, m, r0)
    vr = source_ring_vr(xi, rho, m, r0)

    #=
    assemble output components in the format:
        [x_j r_j; x_{j+1} r_{j+1}]
    and scale by influence panel length
        (due to transformation of integration range to/from range=(0,1))
     =#
    # return StaticArrays.SMatrix{2,2}([vz * [1.0 - t; t] vr * [1.0 - t; t]] * influence_length)
    return StaticArrays.SVector{4}(
        [vz * (1.0 - t); vz * t; vr * (1.0 - t); vr * t] * influence_length
    )
end

"""
`V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: [vz_j vr_j; vz_{j+1} vr_{j+1}]
"""
function nominal_source_panel_integration(
    node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0], debug=false
)

    # Define function to integrate
    function fsample(t)
        return nominal_source_induced_velocity_sample(
            t, node1, node2, influence_length, controlpoint; nondimrange=nondimrange
        )
    end

    V, err = quadgk(fsample, 0.0, 1.0; maxevals=1e2)

    if debug
        return reshape(V, (2, 2)), err
        # return V, err
    else
        return reshape(V, (2, 2))
        # return V
    end
end

#---------------------------------#
#             Singular            #
#---------------------------------#

"""
"""
function subtracted_singular_source_influence(node, controlpoint)
    rmag2 = (controlpoint[1] - node[1])^2 + (controlpoint[2] - node[2])^2
    #TODO: write up math and check signs here.
    num1z = controlpoint[1] - node[1]
    num1r = controlpoint[2] - node[2]
    den1 = 2.0 * pi * rmag2
    # den2 = 64.0 * controlpoint[2]^2
    den2 = controlpoint[2]^2 #DFDC has no 64 here somehow

    radial = num1r / den1 - log(rmag2 / den2) / (8.0 * pi * controlpoint[2])
    axial = num1z / den1

    return axial, radial
end

"""
"""
function analytically_integrated_source_influence(r, influence_length)
    radial = (influence_length / (4.0 * pi * r)) * (1.0 + log(2.0 * r / influence_length))
    axial = 0.0
    return axial, radial
end

"""

**Arguments:**
- `t::Float` : sample point in range (0,1) selected by quadrature.
"""
function self_source_induced_velocity_sample(
    t, node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
)

    # Transform from (0,1) to actual position on panel
    z0 = linear_transform(nondimrange, [node1[1]; node2[1]], t)
    r0 = linear_transform(nondimrange, [node1[2]; node2[2]], t)

    # get relative geometry
    xi, rho, m, _ = calculate_xrm(controlpoint, [z0; r0])

    # Get velocity components at sample points
    vz = source_ring_vz(xi, rho, m, r0)
    vr = source_ring_vr(xi, rho, m, r0)

    # Get singular piece to subtract
    vzs, vrs = subtracted_singular_source_influence([z0; r0], controlpoint)

    #=
    assemble output components in the format:
        [x_j r_j; x_{j+1} r_{j+1}]
    and scale by influence panel length
        (due to transformation of integration range to/from range=(0,1))
     =#
    return [
        vz * (1.0 - t) .- vzs / 2.0
        vz * t .- vzs / 2.0
        vr * (1.0 - t) .- vrs / 2.0
        vr * t .- vrs / 2.0
    ]
end

"""
`V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: [vz_j vr_j; vz_{j+1} vr_{j+1}]
"""
function self_source_panel_integration(
    node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0], debug=false
)

    # Define function to integrate
    function fsample(t)
        return self_source_induced_velocity_sample(
            t, node1, node2, influence_length, controlpoint; nondimrange=nondimrange
        )
    end

    V, err = quadgk(fsample, 0.0, 0.5, 1.0; maxevals=1e2)

    vza, vra = analytically_integrated_source_influence(controlpoint[2], influence_length)

    V .*= influence_length
    V[3:4] .+= vra / 2.0

    if debug
        return reshape(V, (2, 2)), err
    else
        return reshape(V, (2, 2))
    end
end
