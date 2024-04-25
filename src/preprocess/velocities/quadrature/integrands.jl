######################################################################
#                                                                    #
#                             VORTEX                                 #
#                                                                    #
######################################################################

#---------------------------------#
#             Nominal             #
#---------------------------------#
"""
    nominal_vortex_induced_velocity_sample(
        t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0, 1.0)
    )

Calculate the velocity induced by a linear vortex panel on a point.

# Arguments
- `t::Float` : sample point in range (0,1) selected by quadrature method.
- `node1::Vector{Float}` : first panel node (edge) position.
- `node2::Vector{Float}` : second panel node (edge) position.
- `influence_length::Float` : dimensional length of panel.
- `controlpoint::Vector{Float}` : controlpoint position
- `cache_vec::Vector{Float}` : cache for intermediate calculations

# Keyword Arguments
- `nondimrange::Tuple=(0.0,1.0)` : Non-dimensional range describing the panel length.  Do not change excpet for debugging purposes. Note, can also be a vector.

# Returns
- `V::Matrix{Float}` : 2x2 matrix of axial and radial induced velocities from each of the nodes.
"""
function nominal_vortex_induced_velocity_sample(
    t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0, 1.0)
)

    # Transform from (0,1) to actual position on panel
    cache_vec[1] = linear_transform(nondimrange, (node1[1], node2[1]), t) # z coordinate
    cache_vec[2] = linear_transform(nondimrange, (node1[2], node2[2]), t) # r coordinate

    # get relative geometry: xi, rho, m, r0 = calculate_xrm(controlpoint, [z; r])
    calculate_xrm!(view(cache_vec, 3:6), controlpoint, view(cache_vec, 1:2))

    # Get velocity components at sample points
    cache_vec[7] = vortex_ring_vz!(
        cache_vec[3],#xi
        cache_vec[4],#rho
        cache_vec[5],#m
        cache_vec[6],#rj
        1.0,
        view(cache_vec, 11:15),
    ) #vz

    cache_vec[8] = vortex_ring_vr!(
        cache_vec[3], cache_vec[4], cache_vec[5], cache_vec[6], view(cache_vec, 11:16)
    ) #vr

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
        [
            cache_vec[7] * (1.0 - t)
            cache_vec[7] * t
            cache_vec[8] * (1.0 - t)
            cache_vec[8] * t
        ] * influence_length,
    )
    #  return [vz * [1.0 - t; t] vr * [1.0 - t; t]] * influence_length
end

"""
    nominal_vortex_induced_velocity_sample!(
        V, t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0, 1.0)
    )

In-place version of `nominal_vortex_induced_velocity_sample`.
"""
function nominal_vortex_induced_velocity_sample!(
    V, t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0, 1.0)
)

    # Transform from (0,1) to actual position on panel
    cache_vec[1] = linear_transform(nondimrange, (node1[1], node2[1]), t) # z coordinate
    cache_vec[2] = linear_transform(nondimrange, (node1[2], node2[2]), t) # r coordinate

    # get relative geometry: xi, rho, m, r0 = calculate_xrm(controlpoint, [z; r])
    calculate_xrm!(view(cache_vec, 3:6), controlpoint, view(cache_vec, 1:2))

    # Get velocity components at sample points
    cache_vec[7] = vortex_ring_vz!(
        cache_vec[3],#xi
        cache_vec[4],#rho
        cache_vec[5],#m
        cache_vec[6],#rj
        1.0,
        view(cache_vec, 11:15),
    ) #vz

    cache_vec[8] = vortex_ring_vr!(
        cache_vec[3], cache_vec[4], cache_vec[5], cache_vec[6], view(cache_vec, 11:16)
    ) #vr

    #=
    assemble output components in the format:
        [x_j r_j; x_{j+1} r_{j+1}]
    and scale by influence panel length
        (due to transformation of integration range to/from range=(0,1))
     =#
    V[1] += cache_vec[7] * (1.0 - t) * influence_length
    V[2] += cache_vec[7] * t * influence_length
    V[3] += cache_vec[8] * (1.0 - t) * influence_length
    V[4] += cache_vec[8] * t * influence_length
    return V
end

#---------------------------------#
#             Singular            #
#---------------------------------#

"""
    subtracted_singular_vortex_influence(node, controlpoint)

Calculate the singular portions of the self-induced vortex panel influence to subtract off the integral in the separation of singularity method.

# Arguments
- `node::Vector{Float}` : node position
- `controlpoint::Vector{Float}` : controlpoint position

# Returns
- `axial::Float` : axial direction influence
- `radial::Float` : radial direction influence
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
    subtracted_singular_vortex_influence!(node, controlpoint, cache_vec)

Somewhat in-place `version of subtracted_singular_vortex_influence`.

# Arguments
- `node::Vector{Float}` : node position
- `controlpoint::Vector{Float}` : controlpoint position
- `cache_vec::Vector{Float}` : used to store intermediate values.

# Returns
- `axial::Float` : axial direction influence
- `radial::Float` : radial direction influence
"""
function subtracted_singular_vortex_influence!(node, controlpoint, cache_vec)
    cache_vec[1] = (controlpoint[1] - node[1])^2 + (controlpoint[2] - node[2])^2
    cache_vec[2] = controlpoint[1] - node[1]
    cache_vec[3] = node[2] - controlpoint[2]
    cache_vec[4] = 2.0 * pi * cache_vec[1]
    cache_vec[5] = 64.0 * controlpoint[2]^2

    if isapprox(cache_vec[1], 0.0) || isapprox(controlpoint[2], 0.0)
        axial = 0.0
    else
        axial =
            cache_vec[3] / cache_vec[4] -
            log(cache_vec[1] / cache_vec[5]) / (8.0 * pi * controlpoint[2])
    end

    if isapprox(controlpoint[2], 0.0)
        radial = 0.0
    else
        radial = cache_vec[2] / cache_vec[4]
    end

    return axial, radial
end

"""
    analytically_integrated_vortex_influence(r, influence_length)

Analytical approximation of the singular portions of the self-induced vortex panel velocities to be added back in as part of the separation of singularity method.

# Arguments
- `r::Float` : radial position of self-induced control point
- `influence_length::Float` : dimensional length of the panel

# Returns
- `V::Vector{Float}` : axial and radial induced velocities
"""
function analytically_integrated_vortex_influence(r, influence_length)
    axial = (influence_length / (4.0 * pi * r)) * (1.0 + log(16.0 * r / influence_length))
    radial = 0.0
    return [axial; radial]
end

"""
    analytically_integrated_vortex_influence!(V, r, influence_length)

In-place version of `analytically_integrated_vortex_influence`.
"""
function analytically_integrated_vortex_influence!(V, r, influence_length)
    V[1] = (influence_length / (4.0 * pi * r)) * (1.0 + log(16.0 * r / influence_length))
    V[2] = 0.0
    return V
end

"""
    self_vortex_induced_velocity_sample(
        t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0, 1.0)
    )

Calculate the velocity induced by a linear vortex panel on a point at the midpoint between the panel edges.

# Arguments
- `t::Float` : sample point in range (0,1) selected by quadrature method.
- `node1::Vector{Float}` : first panel node (edge) position.
- `node2::Vector{Float}` : second panel node (edge) position.
- `influence_length::Float` : dimensional length of panel.
- `controlpoint::Vector{Float}` : controlpoint position
- `cache_vec::Vector{Float}` : cache for intermediate calculations

# Keyword Arguments
- `nondimrange::Tuple=(0.0,1.0)` : Non-dimensional range describing the panel length.  Do not change excpet for debugging purposes. Note, can also be a vector.

# Returns
- `V::Matrix{Float}` : 2x2 matrix of axial and radial induced velocities from each of the nodes.
"""
function self_vortex_induced_velocity_sample(
    t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0, 1.0)
)

    # Transform from (0,1) to actual position on panel
    cache_vec[1] = linear_transform(nondimrange, (node1[1], node2[1]), t)
    cache_vec[2] = linear_transform(nondimrange, (node1[2], node2[2]), t)

    # get relative geometry
    calculate_xrm!(view(cache_vec, 3:6), controlpoint, view(cache_vec, 1:2))

    # Get velocity components at sample points
    cache_vec[7] = vortex_ring_vz!(
        cache_vec[3],#xi
        cache_vec[4],#rho
        cache_vec[5],#m
        cache_vec[2],#rj
        1.0,
        view(cache_vec, 11:15),
    ) #vz
    cache_vec[8] = vortex_ring_vr!(
        cache_vec[3], cache_vec[4], cache_vec[5], cache_vec[2], view(cache_vec, 11:16)
    ) #vr

    # Get singular piece to subtract
    cache_vec[9], cache_vec[10] = subtracted_singular_vortex_influence!(
        (cache_vec[1], cache_vec[2]), controlpoint, view(cache_vec, 11:15)
    )

    #=
    assemble output components in the format:
        [x_j r_j; x_{j+1} r_{j+1}]
    and scale by influence panel length
        (due to transformation of integration range to/from range=(0,1))
     =#
    return [
        cache_vec[7] * (1.0 - t) - cache_vec[9] / 2.0
        cache_vec[7] * t .- cache_vec[9] / 2.0
        cache_vec[8] * (1.0 - t) .- cache_vec[10] / 2.0
        cache_vec[8] * t .- cache_vec[10] / 2.0
    ]
end

"""
    self_vortex_induced_velocity_sample!(
        V, t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0, 1.0)
    )

In-place version of `self_vortex_induced_velocity_sample`.
"""
function self_vortex_induced_velocity_sample!(
    V, t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0, 1.0)
)

    # Transform from (0,1) to actual position on panel
    cache_vec[1] = linear_transform(nondimrange, (node1[1], node2[1]), t)
    cache_vec[2] = linear_transform(nondimrange, (node1[2], node2[2]), t)

    # get relative geometry
    calculate_xrm!(view(cache_vec, 3:6), controlpoint, view(cache_vec, 1:2))

    # Get velocity components at sample points
    cache_vec[7] = vortex_ring_vz!(
        cache_vec[3],#xi
        cache_vec[4],#rho
        cache_vec[5],#m
        cache_vec[2],#rj
        1.0,
        view(cache_vec, 11:15),
    ) #vz
    cache_vec[8] = vortex_ring_vr!(
        cache_vec[3], cache_vec[4], cache_vec[5], cache_vec[2], view(cache_vec, 11:16)
    ) #vr

    # Get singular piece to subtract
    cache_vec[9], cache_vec[10] = subtracted_singular_vortex_influence!(
        (cache_vec[1], cache_vec[2]), controlpoint, view(cache_vec, 11:15)
    )

    #=
    assemble output components in the format:
        [x_j r_j; x_{j+1} r_{j+1}]
    and scale by influence panel length
        (due to transformation of integration range to/from range=(0,1))
     =#
    V[1] += cache_vec[7] * (1.0 - t) - cache_vec[9] / 2.0
    V[2] += cache_vec[7] * t .- cache_vec[9] / 2.0
    V[3] += cache_vec[8] * (1.0 - t) .- cache_vec[10] / 2.0
    V[4] += cache_vec[8] * t .- cache_vec[10] / 2.0

    return V
end

######################################################################
#                                                                    #
#                             SOURCE                                 #
#                                                                    #
######################################################################

#---------------------------------#
#             Nominal             #
#---------------------------------#

"""
    nominal_source_induced_velocity_sample(
        t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0, 1.0)
    )

Calculate the velocity induced by a source panel on a point.

# Arguments
- `t::Float` : sample point in range (0,1) selected by quadrature method.
- `node1::Vector{Float}` : first panel node (edge) position.
- `node2::Vector{Float}` : second panel node (edge) position.
- `influence_length::Float` : dimensional length of panel.
- `controlpoint::Vector{Float}` : controlpoint position
- `cache_vec::Vector{Float}` : cache for intermediate calculations

# Keyword Arguments
- `nondimrange::Tuple=(0.0,1.0)` : Non-dimensional range describing the panel length.  Do not change excpet for debugging purposes. Note, can also be a vector.

# Returns
- `V::Matrix{Float}` : 2x2 matrix of axial and radial induced velocities from each of the nodes.
"""
function nominal_source_induced_velocity_sample(
    t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0, 1.0)
)

    # Transform from (0,1) to actual position on panel
    cache_vec[1] = linear_transform(nondimrange, (node1[1], node2[1]), t)
    cache_vec[2] = linear_transform(nondimrange, (node1[2], node2[2]), t)

    # get relative geometry
    calculate_xrm!(view(cache_vec, 3:6), controlpoint, view(cache_vec, 1:2))

    # Get velocity components at sample points
    cache_vec[7] = source_ring_vz!(
        cache_vec[3],#xi
        cache_vec[4],#rho
        cache_vec[5],#m
        cache_vec[6],#rj
        view(cache_vec, 11:15),
    )
    cache_vec[8] = source_ring_vr!(
        cache_vec[3], cache_vec[4], cache_vec[5], cache_vec[6], view(cache_vec, 11:16)
    )

    #=
    assemble output components in the format:
        [x_j r_j; x_{j+1} r_{j+1}]
    and scale by influence panel length
        (due to transformation of integration range to/from range=(0,1))
     =#
    # return StaticArrays.SMatrix{2,2}([cache_vec[7] * [1.0 - t; t] cache_vec[8] * [1.0 - t; t]] * influence_length)
    return StaticArrays.SVector{4}(
        [
            cache_vec[7] * (1.0 - t)
            cache_vec[7] * t
            cache_vec[8] * (1.0 - t)
            cache_vec[8] * t
        ] * influence_length,
    )
end

"""
    nominal_source_induced_velocity_sample!(
        V, t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0,1.0)
    )

In-place version of `nominal_source_induced_velocity_sample`.
"""
function nominal_source_induced_velocity_sample!(
    V, t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0, 1.0)
)

    # Transform from (0,1) to actual position on panel
    cache_vec[1] = linear_transform(nondimrange, (node1[1], node2[1]), t)
    cache_vec[2] = linear_transform(nondimrange, (node1[2], node2[2]), t)

    # get relative geometry
    calculate_xrm!(view(cache_vec, 3:6), controlpoint, view(cache_vec, 1:2))

    # Get velocity components at sample points
    cache_vec[7] = source_ring_vz!(
        cache_vec[3],#xi
        cache_vec[4],#rho
        cache_vec[5],#m
        cache_vec[6],#rj
        view(cache_vec, 11:15),
    )
    cache_vec[8] = source_ring_vr!(
        cache_vec[3], cache_vec[4], cache_vec[5], cache_vec[6], view(cache_vec, 11:16)
    )

    #=
    assemble output components in the format:
        [x_j r_j; x_{j+1} r_{j+1}]
    and scale by influence panel length
        (due to transformation of integration range to/from range=(0,1))
     =#
    V[1] += cache_vec[7] * (1.0 - t) * influence_length
    V[2] += cache_vec[7] * t * influence_length
    V[3] += cache_vec[8] * (1.0 - t) * influence_length
    V[4] += cache_vec[8] * t * influence_length

    return V
end

#---------------------------------#
#             Singular            #
#---------------------------------#

"""
    subtracted_singular_source_influence(node, controlpoint)

Calculate the singular portions of the self-induced source panel influence to subtract off the integral in the separation of singularity method.

# Arguments
- `node::Vector{Float}` : node position
- `controlpoint::Vector{Float}` : controlpoint position

# Returns
- `axial::Float` : axial direction influence
- `radial::Float` : radial direction influence
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
    subtracted_singular_source_influence!(node, controlpoint, cache_vec)

In-place version of `subtracted_singular_source_influence`.
"""
function subtracted_singular_source_influence!(node, controlpoint, cache_vec)
    cache_vec[1] = (controlpoint[1] - node[1])^2 + (controlpoint[2] - node[2])^2
    #TODO: write up math and check signs here.
    cache_vec[2] = controlpoint[1] - node[1]
    cache_vec[3] = controlpoint[2] - node[2]
    cache_vec[4] = 2.0 * pi * cache_vec[1]
    # den2 = 64.0 * controlpoint[2]^2
    cache_vec[5] = controlpoint[2]^2 #DFDC has no 64 here somehow

    radial =
        cache_vec[3] / cache_vec[4] -
        log(cache_vec[1] / cache_vec[5]) / (8.0 * pi * controlpoint[2])
    axial = cache_vec[2] / cache_vec[4]

    return axial, radial
end

"""
    analytically_integrated_source_influence(r, influence_length)

Analytical approximation of the singular portions of the self-induced source panel velocities to be added back in as part of the separation of singularity method.

# Arguments
- `r::Float` : radial position of self-induced control point
- `influence_length::Float` : dimensional length of the panel

# Returns
- `V::Vector{Float}` : axial and radial induced velocities
"""
function analytically_integrated_source_influence!(V, r, influence_length)
    V[2] = (influence_length / (4.0 * pi * r)) * (1.0 + log(2.0 * r / influence_length))
    V[1] = 0.0
    return V
end

"""
    analytically_integrated_source_influence(r, influence_length)

In-place version of `analytically_integrated_source_influence`.
"""
function analytically_integrated_source_influence(r, influence_length)
    radial = (influence_length / (4.0 * pi * r)) * (1.0 + log(2.0 * r / influence_length))
    axial = 0.0
    return axial, radial
end

"""
    self_source_induced_velocity_sample(
        t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0, 1.0)
    )

Calculate the velocity induced by a linear source panel on a point at the midpoint between the panel edges.

# Arguments
- `t::Float` : sample point in range (0,1) selected by quadrature method.
- `node1::Vector{Float}` : first panel node (edge) position.
- `node2::Vector{Float}` : second panel node (edge) position.
- `influence_length::Float` : dimensional length of panel.
- `controlpoint::Vector{Float}` : controlpoint position
- `cache_vec::Vector{Float}` : cache for intermediate calculations

# Keyword Arguments
- `nondimrange::Tuple=(0.0,1.0)` : Non-dimensional range describing the panel length.  Do not change excpet for debugging purposes. Note, can also be a vector.

# Returns
- `V::Matrix{Float}` : 2x2 matrix of axial and radial induced velocities from each of the nodes.
"""
function self_source_induced_velocity_sample(
    t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0, 1.0)
)

    # Transform from (0,1) to actual position on panel
    cache_vec[1] = linear_transform(nondimrange, (node1[1], node2[1]), t)
    cache_vec[2] = linear_transform(nondimrange, (node1[2], node2[2]), t)

    # get relative geometry
    calculate_xrm!(view(cache_vec, 3:6), controlpoint, view(cache_vec, 1:2))

    # Get velocity components at sample points
    cache_vec[7] = source_ring_vz!(
        cache_vec[3],#xi
        cache_vec[4],#rho
        cache_vec[5],#m
        cache_vec[2],#rj
        view(cache_vec, 11:15),
    )
    cache_vec[8] = source_ring_vr!(
        cache_vec[3], cache_vec[4], cache_vec[5], cache_vec[2], view(cache_vec, 11:16)
    )

    # Get singular piece to subtract
    cache_vec[9], cache_vec[10] = subtracted_singular_source_influence!(
        (cache_vec[1], cache_vec[2]), controlpoint, view(cache_vec, 11:15)
    )

    #=
    assemble output components in the format:
        [x_j r_j; x_{j+1} r_{j+1}]
    and scale by influence panel length
        (due to transformation of integration range to/from range=(0,1))
     =#
    return [
        cache_vec[7] * (1.0 - t) .- cache_vec[9] / 2.0
        cache_vec[7] * t .- cache_vec[9] / 2.0
        cache_vec[8] * (1.0 - t) .- cache_vec[10] / 2.0
        cache_vec[8] * t .- cache_vec[10] / 2.0
    ]
end

"""
    self_source_induced_velocity_sample!(
        V, t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0, 1.0)
    )

In-place version of `self_source_induced_velocity_sample`.
"""
function self_source_induced_velocity_sample!(
    V, t, node1, node2, influence_length, controlpoint, cache_vec; nondimrange=(0.0, 1.0)
)

    # Transform from (0,1) to actual position on panel
    cache_vec[1] = linear_transform(nondimrange, (node1[1], node2[1]), t)
    cache_vec[2] = linear_transform(nondimrange, (node1[2], node2[2]), t)

    # get relative geometry
    calculate_xrm!(view(cache_vec, 3:6), controlpoint, view(cache_vec, 1:2))

    # Get velocity components at sample points
    cache_vec[7] = source_ring_vz!(
        cache_vec[3],#xi
        cache_vec[4],#rho
        cache_vec[5],#m
        cache_vec[2],#rj
        view(cache_vec, 11:15),
    )
    cache_vec[8] = source_ring_vr!(
        cache_vec[3], cache_vec[4], cache_vec[5], cache_vec[2], view(cache_vec, 11:16)
    )

    # Get singular piece to subtract
    cache_vec[9], cache_vec[10] = subtracted_singular_source_influence!(
        (cache_vec[1], cache_vec[2]), controlpoint, view(cache_vec, 11:15)
    )

    #=
    assemble output components in the format:
        [x_j r_j; x_{j+1} r_{j+1}]
    and scale by influence panel length
        (due to transformation of integration range to/from range=(0,1))
     =#
    V[1] += cache_vec[7] * (1.0 - t) .- cache_vec[9] / 2.0
    V[2] += cache_vec[7] * t .- cache_vec[9] / 2.0
    V[3] += cache_vec[8] * (1.0 - t) .- cache_vec[10] / 2.0
    V[4] += cache_vec[8] * t .- cache_vec[10] / 2.0

    return V
end

