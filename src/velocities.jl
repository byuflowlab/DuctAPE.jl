######################################################################
#                                                                    #
#                         "MESH" GENERATION                          #
#                                                                    #
######################################################################
"""

calculate "mesh" geometry without creating a mesh object
"""
function calculate_xrm(node, fieldpoint)
    if isapprox(node[2], 0.0)
        return 0.0, 0.0, 0.0, 0.0
    else
        # normalized axial distance
        xi = (fieldpoint[1] - node[1]) / node[2]

        # normalized radial distance
        rho = fieldpoint[2] / node[2]

        # elliptic integral parameter
        m = (4.0 * rho) / (xi^2 + (rho + 1)^2)

        # influence point radial position
        rj = node[2]

        return xi, rho, m, rj
    end
end

######################################################################
#                                                                    #
#                         INDUCED VELOCITIES                         #
#                                                                    #
######################################################################
#---------------------------------#
#    Unit Induced Velocities      #
#---------------------------------#

##### ----- Vortex ----- #####
"""
"""
function vortex_ring_vx(xi, rho, m, r_influence, len)

    # check panel locations
    if isapprox(r_influence, 0.0)
        # if influence on the axis, the influence is set to zero
        return 0.0
    elseif (xi^2 + (rho - 1.0)^2 <= eps())
        # set self-induced case is "smoke ring" self influence in axial direction only.
        return -1.0 / (4.0 * pi * rho) * (log(8.0 * pi * rho / len) - 0.25)
    else
        #get the first denominator
        den1 = 2.0 * pi * r_influence * sqrt(xi^2 + (rho + 1.0)^2)

        #get numerator and denominator of second fraction
        num2 = 2.0 * (rho - 1.0)
        den2 = xi^2 + (rho - 1.0)^2

        #get values for elliptic integrals
        K, E = get_elliptics(m)

        # negative here is due to our convention that the vortex is postive clockwise
        return -1.0 / den1 * (K - (1.0 + num2 / den2) * E)
    end
end

"""
"""
function vortex_ring_vr(xi, rho, m, r_influence)

    # return 0.0 for self-induced, influence on axis, or target on axis cases
    if (xi^2 + (rho - 1.0)^2 <= eps()) || isapprox(r_influence, 0.0) || isapprox(rho, 0.0)
        return 0.0
    else
        #get numerator and denominator of first fraction
        num1 = xi / rho
        den1 = 2.0 * pi * r_influence * sqrt(xi^2 + (rho + 1.0)^2)

        num2 = 2.0 * rho
        den2 = xi^2 + (rho - 1.0)^2

        #get values for elliptic integrals
        K, E = get_elliptics(m)

        return num1 / den1 * (K - (1.0 + num2 / den2) * E)
    end
end

##### ----- Source ----- #####
"""
"""
function source_ring_vx(xi, rho, m, r_influence)

    # return zero for the self-induced off body case
    if (xi^2 + (rho - 1.0)^2 <= eps()) || isapprox(r_influence, 0.0) || isapprox(rho, 0.0)
        return 0.0
    else

        #get values for elliptic integrals
        K, E = get_elliptics(m)

        #get the first denominator
        den1 = 2.0 * pi * r_influence * sqrt(xi^2 + (rho + 1.0)^2)

        #get numerator and denominator of second fraction
        num2 = 2 * xi * E
        den2 = xi^2 + (rho - 1)^2

        return 1.0 / den1 * (num2 / den2)
    end
end

"""
"""
function source_ring_vr(xi, rho, m, r_influence)

    # return zero for the self-induced off-body case
    if (xi^2 + (rho - 1.0)^2 <= eps()) || isapprox(r_influence, 0.0) || isapprox(rho, 0.0)
        return 0.0
    else

        #get values for elliptic integrals
        K, E = get_elliptics(m)

        #get numerator and denominator of first fraction
        num1 = 1.0 / rho
        den1 = 2.0 * pi * r_influence * sqrt(xi^2 + (rho + 1.0)^2)

        #get numerator and denominator of second fraction
        num2 = 2 * rho * (rho - 1.0)
        den2 = xi^2 + (rho - 1)^2

        return num1 / den1 * (K - (1.0 - num2 / den2) * E)
    end
end

#---------------------------------#
#   "Panel" Induced Velocities    #
#---------------------------------#

##### ----- Vortex ----- #####
"""
"""
function constant_vortex_band_induced_velocity!(
    vel, controlpoint, influence_panel_length, fieldpoint, gamma=1.0
)

    # get relative geometry
    xi, rho, m, rj = calculate_xrm(controlpoint, fieldpoint)

    # get unit induced velocities of vortex ring
    vel[1] += vortex_ring_vx(xi, rho, m, rj, influence_panel_length)#length needed here in case of self-induced case
    vel[2] += vortex_ring_vr(xi, rho, m, rj)

    # multiply by panel length and strength
    vel[:] .*= gamma * influence_panel_length

    return nothing
end

"""
"""
function constant_vortex_band_induced_velocity(
    controlpoint::AbstractVector{T1},
    influence_panel_length::T2,
    fieldpoint::AbstractVector{T3},
    gamma::T4=1.0,
) where {T1,T2,T3,T4}

    # Initialize
    T = promote_type(T1, T2, T3, T4)
    vel = zeros(T, 2)

    # get velocities
    constant_vortex_band_induced_velocity!(
        vel, controlpoint, influence_panel_length, fieldpoint, gamma
    )

    return vel
end

"""
"""
function influencefromvortexpanels(
    fieldpoints::AbstractMatrix{T1},
    controlpoints::AbstractArray{T2},
    influence_panel_lengths::AbstractVector{T3},
    strengths::AbstractArray{T4},
) where {T1,T2,T3,T4}
    T = promote_type(T1, T2, T3, T4)
    AIC = zeros(T, size(fieldpoints, 1), size(controlpoints, 1), 2)

    influencefromvortexpanels!(
        AIC, fieldpoints, controlpoints, influence_panel_lengths, strengths
    )

    return AIC
end

"""
"""
function influencefromvortexpanels!(
    AIC, fieldpoints, controlpoints, influence_panel_lengths, strengths
)
    for (i, fp) in enumerate(eachrow(fieldpoints))
        # loop through panels doing the influencing
        for (j, (gamma, cpj, lj)) in
            enumerate(zip(strengths, eachrow(controlpoints), influence_panel_lengths))

            # get unit induced velocity from the panel onto the control point
            constant_vortex_band_induced_velocity!(view(AIC, i, j, :), cpj, lj, fp, gamma)
        end
    end

    return nothing
end

"""
"""
function vfromvortexpanels(
    fieldpoints::AbstractMatrix{T1},
    controlpoints::AbstractArray{T2},
    influence_panel_lengths::AbstractVector{T3},
    strengths::AbstractArray{T4},
) where {T1,T2,T3,T4}
    T = promote_type(T1, T2, T3, T4)
    V = zeros(T, size(fieldpoints, 1), 2)

    vfromvortexpanels!(V, fieldpoints, controlpoints, influence_panel_lengths, strengths)

    return V
end

"""
"""
function vfromvortexpanels!(
    V, fieldpoints, controlpoints, influence_panel_lengths, strengths
)
    for (i, (fp, vel)) in enumerate(zip(eachrow(fieldpoints), eachrow(V)))
        # loop through panels doing the influencing
        for (j, (gamma, cpj, lj)) in
            enumerate(zip(strengths, eachrow(controlpoints), influence_panel_lengths))

            # get unit induced velocity from the panel onto the control point
            constant_vortex_band_induced_velocity!(vel, cpj, lj, fp, gamma)
        end
    end

    return nothing
end

##### ----- Source ----- #####
"""
"""
function constant_source_band_induced_velocity!(
    vel, node, influence_panel_length, fieldpoint, sigma=1.0
)

    # get relative geometry
    xi, rho, m, rj = calculate_xrm(node, fieldpoint)

    # get unit induced velocities of source ring
    vel[1] += source_ring_vx(xi, rho, m, rj)
    vel[2] += source_ring_vr(xi, rho, m, rj)

    # multiply by panel length and strength
    vel[:] .*= sigma * influence_panel_length

    return nothing
end

"""
"""
function constant_source_band_induced_velocity(
    controlpoint::AbstractVector{T1},
    influence_panel_length::T2,
    fieldpoint::AbstractVector{T3},
    sigma::T4=1.0,
) where {T1,T2,T3,T4}

    # Initialize
    T = promote_type(T1, T2, T3, T4)
    vel = zeros(T, 2)

    # get velocities
    constant_source_band_induced_velocity!(
        vel, controlpoint, influence_panel_length, fieldpoint, sigma
    )

    return vel
end

"""
"""
function influencefromsourcepanels(
    fieldpoints::AbstractMatrix{T1},
    controlpoints::AbstractArray{T2},
    influence_panel_lengths::AbstractVector{T3},
    strengths::AbstractArray{T4},
) where {T1,T2,T3,T4}
    T = promote_type(T1, T2, T3, T4)
    AIC = zeros(T, size(fieldpoints, 1), size(controlpoints, 1), 2)

    influencefromsourcepanels!(
        AIC, fieldpoints, controlpoints, influence_panel_lengths, strengths
    )

    return AIC
end

"""
"""
function influencefromsourcepanels!(
    AIC, fieldpoints, controlpoints, influence_panel_length, strengths
)
    for (i, fp) in enumerate(eachrow(fieldpoints))
        # loop through panels doing the influencing
        for (j, (sigma, cpj, lj)) in
            enumerate(zip(strengths, eachrow(controlpoints), influence_panel_length))

            # get unit induced velocity from the panel onto the control point
            constant_source_band_induced_velocity!(view(AIC, i, j, :), cpj, lj, fp, sigma)
        end
    end

    return nothing
end

"""
"""
function vfromsourcepanels(
    fieldpoints::AbstractMatrix{T1},
    controlpoints::AbstractArray{T2},
    influence_panel_lengths::AbstractVector{T3},
    strengths::AbstractArray{T4},
) where {T1,T2,T3,T4}
    T = promote_type(T1, T2, T3, T4)
    V = zeros(T, size(fieldpoints, 1), 2)

    vfromsourcepanels!(V, fieldpoints, controlpoints, influence_panel_lengths, strengths)

    return V
end

"""
"""
function vfromsourcepanels!(
    V, fieldpoints, controlpoints, influence_panel_length, strengths
)
    for (i, (fp, vel)) in enumerate(zip(eachrow(fieldpoints), eachrow(V)))
        # loop through panels doing the influencing
        for (j, (sigma, cpj, lj)) in
            enumerate(zip(strengths, eachrow(controlpoints), influence_panel_length))

            # get unit induced velocity from the panel onto the control point
            constant_source_band_induced_velocity!(vel, cpj, lj, fp, sigma)
        end
    end

    return nothing
end

##### ----- Doublet ----- #####
"""
"""
function constant_doublet_band_induced_velocity!(vel, node1, node2, fieldpoint, mu=1.0)

    # - Get Relative Geometry - #
    xi1, rho1, m1, rj1 = calculate_xrm(node1, fieldpoint)
    xi2, rho2, m2, rj2 = calculate_xrm(node2, fieldpoint)

    #NOTE: this seems a little sketchy.  things don't work if we set the whole induced velocity to zero if one edge of the panel is on the axis, but if we don't do that, then the contribution comes from only one of the panel edges, which seems to violate the definition of the doublet panel in the first place...
    # Check location. If neither panel edge lies on axis, continue, otherwise, keep panel influence as zero
    # if !isapprox(rj1, 0.0) && !isapprox(rj2, 0.0)
    # - influence due to node1 - #
    vel[1] += -mu * vortex_ring_vx(xi1, rho1, m1, rj1, 19.5733 * rho1)#length shouldn't be needed here, since there shouldn't be any self-inducement.  Set length such that log(x) = 0.25.
    vel[2] += -mu * vortex_ring_vr(xi1, rho1, m1, rj1)

    # - influence due to node2 - #
    vel[1] += mu * vortex_ring_vx(xi2, rho2, m2, rj2, 19.5733 * rho2)#length shouldn't be needed here, since there shouldn't be any self-inducement.  Set length such that log(x) = 0.25.
    vel[2] += mu * vortex_ring_vr(xi2, rho2, m2, rj2)
    # end

    return nothing
end

"""
"""
function constant_doublet_band_induced_velocity(
    node1::AbstractVector{T1},
    node2::AbstractVector{T2},
    fieldpoint::AbstractVector{T3},
    mu::T4=1.0,
) where {T1,T2,T3,T4}

    # initialize
    T = promote_type(T1, T2, T3, T4)
    vel = zeros(T, 2)

    # get velocities
    constant_doublet_band_induced_velocity!(vel, node1, node2, fieldpoint, mu)

    return vel
end

"""
"""
function influencefromdoubletpanels(
    fieldpoints::AbstractMatrix{T1}, nodes::AbstractArray{T2}, strengths::AbstractArray{T3}
) where {T1,T2,T3}
    T = promote_type(T1, T2, T3)
    AIC = zeros(T, size(fieldpoints, 1), size(nodes, 1), 2)

    influencefromdoubletpanels!(AIC, fieldpoints, nodes, strengths)

    return AIC
end

"""
"""
function influencefromdoubletpanels!(AIC, fieldpoints, nodes, strengths)

    # Loop through field points
    for (i, fp) in enumerate(eachrow(fieldpoints))
        # loop through panels doing the influencing
        for (j, (mu, p1, p2)) in
            enumerate(zip(strengths, eachrow(nodes[:, 1, :]), eachrow(nodes[:, 2, :])))

            # get unit induced velocity from the panel onto the control point
            constant_doublet_band_induced_velocity!(view(AIC, i, j, :), p1, p2, fp, mu)
        end
    end

    return nothing
end

"""
"""
function vfromdoubletpanels(
    fieldpoints::AbstractMatrix{T1}, nodes::AbstractArray{T2}, strengths::AbstractArray{T3}
) where {T1,T2,T3}
    T = promote_type(T1, T2, T3)
    V = zeros(T, size(fieldpoints, 1), 2)

    vfromdoubletpanels!(V, fieldpoints, nodes, strengths)

    return V
end

"""
"""
function vfromdoubletpanels!(V, fieldpoints, nodes, strengths)

    # Loop through field points
    for (i, (fp, vel)) in enumerate(zip(eachrow(fieldpoints), eachrow(V)))
        # loop through panels doing the influencing
        for (j, (mu, p1, p2)) in
            enumerate(zip(strengths, eachrow(nodes[:, 1, :]), eachrow(nodes[:, 2, :])))

            # get unit induced velocity from the panel onto the control point
            constant_doublet_band_induced_velocity!(vel, p1, p2, fp, mu)
        end
    end

    return nothing
end

"""
note that strenghts here is the entire array of body strengths
"""
function influencefromTE!(
    AIC, fieldpoints, TEnodes, strengths; tol=1e1 * eps(), verbose=false
)
    for (j, te) in enumerate(TEnodes)

        # check that trailing edge points are coincident
        mu = te.sign * strengths[te.idx]

        # Loop through control points being influenced
        for (i, fp) in enumerate(eachrow(fieldpoints))

            # influence due to lower TE
            xi, rho, k2, rj = calculate_xrm(te.pos, fp)
            if isapprox(rj, 0.0)
                AIC[i, j, :] += [0.0; 0.0]
            else
                vx = mu * vortex_ring_vx(xi, rho, k2, rj, 19.5733 * rho)#length shouldn't be needed, set it such that self-induced term would be zero
                vr = mu * vortex_ring_vr(xi, rho, k2, rj)
                AIC[i, j, :] += [vx; vr]
            end
        end
    end

    return nothing
end

"""
note that strenghts here is the entire array of body strengths
"""
function influencefromTE(
    fieldpoints::AbstractMatrix{T1},
    TEnodes,
    strengths::AbstractVector{T2};
    tol=1e1 * eps(),
    verbose=false,
) where {T1,T2}
    T = promote_type(T1, T2)
    AIC = zeros(T, size(fieldpoints, 1), length(TEnodes), 2)

    influencefromTE!(AIC, fieldpoints, TEnodes, strengths; tol=tol, verbose=verbose)

    return AIC
end

"""
"""
function vfromTE!(V, fieldpoints, TEnodes, strengths; tol=1e1 * eps(), verbose=false)
    for (i, te) in enumerate(TEnodes)

        # check that trailing edge points are coincident
        mu = te.sign * strengths[te.idx]

        # Loop through control points being influenced
        for (m, (vel, fp)) in enumerate(zip(eachrow(V), eachrow(fieldpoints)))

            # influence due to TE node
            xi, rho, k2, rj = calculate_xrm(te.pos, fp)
            if isapprox(rj, 0.0)
                vel[:] += [0.0; 0.0]
            else
                vx = mu * vortex_ring_vx(xi, rho, k2, rj, 19.5733 * rho)#length shouldn't be needed, set it such that self-induced term would be zero
                vr = mu * vortex_ring_vr(xi, rho, k2, rj)
                vel[:] += [vx; vr]
            end
        end
    end

    return nothing
end

"""
"""
function vfromTE(
    fieldpoints::AbstractMatrix{T1},
    TEnodes,
    strengths::AbstractVector{T2};
    tol=1e1 * eps(),
    verbose=false,
) where {T1,T2}
    T = promote_type(T1, T2)
    V = zeros(T, size(fieldpoints, 1), 2)

    vfromTE!(V, fieldpoints, TEnodes, strengths; tol=tol, verbose=verbose)

    return V
end

"""
"""
function vfromgradmu(panels, strengths::AbstractVector{T}) where {T}

    # T = promote_type(T1, T2)
    V = zeros(T, panels.npanels, 2)

    vfromgradmu!(V, panels, strengths)

    return V
end

"""
"""
function calc_gradmu(gradmu, len1, cpi, cpi1, mui, mui1)

    # TODO: also somewhat handwavy to set infinite gradients to zero
    if !isapprox(cpi[1] - cpi1[1], 0.0)
        gradmu[1] = len1 * (mui - mui1) ./ (cpi[1] - cpi1[1])
    end
    if !isapprox(cpi[2] - cpi1[2], 0.0)
        gradmu[2] = len1 * (mui - mui1) ./ (cpi[2] - cpi1[2])
    end

    return nothing
end

"""
"""
function vfromgradmu!(V, panels, mu)
    len = panels.len
    cp = panels.controlpoint
    nhat = panels.normal
    that = panels.tangent
    nodes = panels.nodes
    np = panels.npanels
    endpoints = panels.endpoints
    endpointidxs = panels.endpointidxs

    # Loop through panels
    for (i, vel) in enumerate(eachrow(V))

        #check if first panel on a body
        #TODO: does gradient really wrap around trailing edge?
        #TODO: if so, need another check to see if the nodes coincide
        gradmu1 = [0.0; 0.0]
        gradmu2 = [0.0; 0.0]

        # check if on a trailing edge
        if i ∈ endpointidxs
            idx = findall(x -> x == i, endpointidxs)

            if !isapprox(endpoints[idx[1][1], 1, :], endpoints[idx[1][1], 2, :])
                if i ∉ endpointidxs[:, 1]
                    # if !isapprox(cp[i,1]-cp[i-1,1], 0.0)
                    #     gradmu1[1] = len[i-1]*(mu[i]-mu[i-1])./(cp[i,1]-cp[i-1,1])
                    # end
                    # if !isapprox(cp[i,2]-cp[i-1,2], 0.0)
                    #     gradmu1[2] = len[i-1]*(mu[i]-mu[i-1])./(cp[i,2]-cp[i-1,2])
                    # end
                    calc_gradmu(
                        gradmu1, len[i - 1], cp[i, :], cp[i - 1, :], mu[i], mu[i - 1]
                    )

                    #check if last panel on a body
                elseif i ∉ endpointidxs[:, 2]
                    # if !isapprox(cp[i,1]-cp[i+1,1], 0.0)
                    #     gradmu2[1] = len[i+1]*(mu[i]-mu[i+1])./(cp[i,1]-cp[i+1,1])
                    # end
                    # if !isapprox(cp[i,2]-cp[i+1,2], 0.0)
                    #     gradmu2[2] = len[i+1]*(mu[i]-mu[i+1])./(cp[i,2]-cp[i+1,2])
                    # end
                    calc_gradmu(
                        gradmu2, len[i + 1], cp[i, :], cp[i + 1, :], mu[i], mu[i + 1]
                    )
                end
            end

        else
            calc_gradmu(gradmu1, len[i - 1], cp[i, :], cp[i - 1, :], mu[i], mu[i - 1])
            calc_gradmu(gradmu2, len[i + 1], cp[i, :], cp[i + 1, :], mu[i], mu[i + 1])
        end

        # get total weighting
        wt = if i in endpointidxs[:, 1]
            len[i + 1]
        elseif i in endpointidxs[:, 2]
            len[i - 1]
        else
            len[i - 1] + len[i + 1]
        end

        # calculate weighted gradient
        # TODO: this is handwavy solution to large gradients...
        # gradmu = (gradmu1+gradmu2).*(1.0 .- nhat[i,:].^2)/wt
        gradmu = (gradmu1 + gradmu2) .* (that[i, :] .^ 2) / wt

        if norm(gradmu) > 1e6
            gradmu = [0.0; 0.0]
        end

        vel[:] -= gradmu / 2.0
    end

    return nothing
end

#---------------------------------#
# Influence Coefficient Matrices  #
#---------------------------------#

##### ----- Vortex ----- #####
"""
"""
function vortex_panel_influence_matrix(influencers, targets)

    # extract influence panel data
    controlpoint = influencers.controlpoint
    controllength = influencers.len

    # extract target panel data
    fieldpoint = targets.controlpoint
    fieldnormal = targets.normal

    T = promote_type(eltype(controlpoint), eltype(fieldpoint))
    M = size(fieldpoint, 1)
    N = size(controlpoint, 1)

    AIC = zeros(T, M, N)

    vortex_panel_influence_matrix!(
        AIC, controlpoint, controllength, fieldpoint, fieldnormal
    )

    return AIC
end

"""
"""
function vortex_panel_influence_matrix!(
    AIC, controlpoint, controllength, fieldpoint, fieldnormal
)

    # Loop through control points being influenced
    for (i, (fp, nhat)) in enumerate(zip(eachrow(fieldpoint), eachrow(fieldnormal)))
        # loop through panels doing the influencing
        for (j, (cpj, lj)) in enumerate(zip(eachrow(controlpoint), controllength))

            # get unit induced velocity from the panel onto the control point
            vel = constant_vortex_band_induced_velocity(cpj, lj, fp)

            # fill the matrix
            AIC[i, j] += dot(vel, nhat)
        end
    end

    return nothing
end

##### ----- Source ----- #####
"""
"""
function source_panel_influence_matrix(influencers, targets)

    # extract influence panel data
    controlpoint = influencers.controlpoint
    controllength = influencers.len

    # extract target panel data
    fieldpoint = targets.controlpoint
    fieldnormal = targets.normal

    T = promote_type(eltype(controlpoint), eltype(fieldpoint))
    M = size(fieldpoint, 1)
    N = size(controlpoint, 1)

    AIC = zeros(T, M, N)

    source_panel_influence_matrix!(
        AIC, controlpoint, controllength, fieldpoint, fieldnormal
    )

    return AIC
end

"""
"""
function source_panel_influence_matrix!(
    AIC, controlpoint, controllength, fieldpoint, fieldnormal
)

    # Loop through control points being influenced
    for (i, (fp, nhat)) in enumerate(zip(eachrow(fieldpoint), eachrow(fieldnormal)))
        # loop through panels doing the influencing
        for (j, (cpj, lj)) in enumerate(zip(eachrow(controlpoint), controllength))

            # get unit induced velocity from the panel onto the control point
            vel = constant_source_band_induced_velocity(cpj, lj, fp)

            # fill the matrix
            AIC[i, j] += dot(vel, nhat)
        end
    end

    return nothing
end

##### ----- Doublet ----- #####
"""
"""
function doublet_panel_influence_matrix(nodes, targets)
    (; controlpoint, normal) = targets

    T = promote_type(eltype(nodes), eltype(controlpoint))
    M = size(controlpoint, 1)
    N = size(nodes, 1)

    AIC = zeros(T, M, N)

    doublet_panel_influence_matrix!(AIC, nodes, controlpoint, normal)

    return AIC
end

"""
"""
function doublet_panel_influence_matrix!(AIC, nodes, controlpoint, normal)

    # Loop through control points being influenced
    for (i, (cp, nhat)) in enumerate(zip(eachrow(controlpoint), eachrow(normal)))
        # loop through panels doing the influencing
        for (j, (p1, p2)) in
            enumerate(zip(eachrow(nodes[:, 1, :]), eachrow(nodes[:, 2, :])))

            # get unit induced velocity from the panel onto the control point
            vel = constant_doublet_band_induced_velocity(p1, p2, cp)

            # fill the matrix
            AIC[i, j] += dot(vel, nhat)
        end
    end

    return nothing
end

##### ----- Freestream ----- #####
"""
"""
function freestream_influence_vector(
    normals::AbstractMatrix{T1}, Vinfmat::AbstractMatrix{T2}
) where {T1,T2}
    T = promote_type(T1, T2)
    N = size(normals, 1)

    RHS = zeros(T, N)

    freestream_influence_vector!(RHS, normals, Vinfmat)

    return RHS
end

"""
"""
function freestream_influence_vector!(RHS, normals, Vinfmat)
    for (i, (n, v)) in enumerate(zip(eachrow(normals), eachrow(Vinfmat)))
        RHS[i] -= dot(v, n)
    end

    return nothing
end
