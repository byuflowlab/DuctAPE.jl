######################################################################
#                                                                    #
#                         "MESH" GENERATION                          #
#                                                                    #
######################################################################
"""

calculate "mesh" geometry without creating a mesh object
"""
function calculate_xrm(node, fieldpoint)

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

######################################################################
#                                                                    #
#                         INDUCED VELOCITIES                         #
#                                                                    #
######################################################################
#---------------------------------#
#    Unit Induced Velocities      #
#---------------------------------#
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


#---------------------------------#
#   "Panel" Induced Velocities    #
#---------------------------------#
"""
"""
function constant_vortex_band_induced_velocity!(vel, node, influence_panel_length, fieldpoint, gamma=1.0)

    # get relative geometry
    xi, rho, m, rj = calculate_xrm(node, fieldpoint)

    # get unit induced velocities of vortex ring
    vel[1] += vortex_ring_vx(xi, rho, m, rj)
    vel[2] += vortex_ring_vr(xi, rho, m, rj)

    # multiply by panel length and strength
    vel[:] .*= gamma*influence_panel_length

    return nothing

end

"""
"""
function constant_vortex_band_induced_velocity(node::AbstractVector{T1}, influence_panel_length::T2, fieldpoint::AbstractVector{T3}, gamma::T4=1.0) where {T1, T2, T3, T4}

    # Initialize
    T = promote_type(T1, T2, T3, T4)
    vel = zeros(T, 2)

    # get velocities
    constant_vortex_band_induced_velocity!(vel, node, influence_panel_length, fieldpoint, gamma)

    return vel

end

"""
"""
function constant_doublet_band_induced_velocity!(vel, node1, node2, fieldpoint, mu=1.0)

    # influence due to node1
    xi, rho, m, rj = calculate_xrm(node1, fieldpoint)
    if !isapprox(rj, 0.0)
        vel[1] += -mu*vortex_ring_vx(xi, rho, m, rj)
        vel[2] += -mu*vortex_ring_vr(xi, rho, m, rj)
    end

    # influence due to node1
    xi, rho, m, rj = calculate_xrm(node2, fieldpoint)
    if !isapprox(rj, 0.0)
        vel[1] += mu*vortex_ring_vx(xi, rho, m, rj)
        vel[2] += mu*vortex_ring_vr(xi, rho, m, rj)
    end

    return nothing

end

"""
"""
function constant_doublet_band_induced_velocity(node1::AbstractVector{T1}, node2::AbstractVector{T2}, fieldpoint::AbstractVector{T3}, mu::T4=1.0) where {T1, T2, T3, T4}

    # initialize
    T = promote_type(T1, T2, T3, T4)
    vel = zeros(T, 2)

    # get velocities
    constant_doublet_band_induced_velocity!(vel, node1, node2, fieldpoint, mu)

    return vel

end


function vfromdoubletpanels(fieldpoints::AbstractMatrix{T1}, nodes::AbstractArray{T2}, strengths::AbstractArray{T3}) where {T1, T2, T3}

    T = promote_type(T1, T2, T3)
    V = zeros(T, size(fieldpoints, 1), 2)

    vfromdoubletpanels!(V, fieldpoints, nodes, strengths)

    return V

end

function vfromdoubletpanels!(V, fieldpoints, nodes, strengths)

    for (i, (fp, vel)) in enumerate(zip(eachrow(fieldpoints), eachrow(V)))
        # loop through panels doing the influencing
        for (j, (mu, p1, p2)) in enumerate(zip(strengths, eachrow(nodes[:,1,:]), eachrow(nodes[:,2,:])))

            # get unit induced velocity from the panel onto the control point
            constant_doublet_band_induced_velocity!(vel, p1, p2, fp, mu)

        end
    end

    return nothing

end

function vfromTE!(V, fieldpoints, TEnodes, TEidxs, strengths; tol = 1e1*eps(), verbose=false)

    for (i, (te1, te2)) in enumerate(zip(eachrow(TEnodes[:,1,:]), eachrow(TEnodes[:,2,:])))

        # check that trailing edge points are coincident
        if norm(te2-te1) < tol

            mulow = strengths[TEidxs[i,1]]
            muupp = strengths[TEidxs[i,2]]

            # Loop through control points being influenced
            for (m, (vel, fp)) in enumerate(zip(eachrow(V), eachrow(fieldpoints)))

                # influence due to lower TE
                xi, rho, k2, rj = calculate_xrm(te1, fp)
                if isapprox(rj, 0.0)
                    vel[:] += [0.0;0.0]
                else
                    vx = mulow*vortex_ring_vx(xi, rho, k2, rj)
                    vr = mulow*vortex_ring_vr(xi, rho, k2, rj)
                    vel[:] += [vx; vr]
                end

                # influence due to upper TE
                xi, rho, k2, rj = calculate_xrm(te2, fp)
                if isapprox(rj, 0.0)
                    vel[:] += [0.0;0.0]
                else
                    vx = mulow*vortex_ring_vx(xi, rho, k2, rj)
                    vr = mulow*vortex_ring_vr(xi, rho, k2, rj)
                    vel[:] += [vx; vr]
                end

            end

        elseif verbose
            @warn "Trailing edge points for body $i are not coincident, no wake velocity added."
        end
    end

    return nothing
end

function vfromTE(fieldpoints::AbstractMatrix{T1}, TEnodes::AbstractArray{T2}, TEidxs, strengths::AbstractVector{T3}; tol = 1e1*eps(), verbose=false) where {T1,T2,T3}

    T = promote_type(T1, T2, T3)
    V = zeros(T, size(fieldpoints, 1), 2)

    vfromTE!(V, fieldpoints, TEnodes, TEidxs, strengths; tol=tol, verbose=verbose)

    return V
end

function vfromgradmu(panels, strengths::AbstractVector{T}) where {T}

    # T = promote_type(T1, T2)
    V = zeros(T, panels.npanels, 2)

    vfromgradmu!(V, panels, strengths)

    return V
end

function calc_gradmu(gradmu, len1, cpi, cpi1, mui, mui1)

        # TODO: also somewhat handwavy to set infinite gradients to zero
    if !isapprox(cpi[1]-cpi1[1], 0.0)
        gradmu[1] = len1*(mui-mui1)./(cpi[1]-cpi1[1])
    end
    if !isapprox(cpi[2]-cpi1[2], 0.0)
        gradmu[2] = len1*(mui-mui1)./(cpi[2]-cpi1[2])
    end

    return nothing
end

function vfromgradmu!(V, panels, mu)

    len = panels.length
    cp = panels.control_point
    nhat = panels.normal
    that = panels.tangent
    nodes = panels.nodes
    np = panels.npanels
    TEnodes = panels.TEnodes
    TEidxs = panels.TEidxs

    # Loop through panels
    for (i, vel) in enumerate(eachrow(V))

        #check if first panel on a body
        #TODO: does gradient really wrap around trailing edge?
        #TODO: if so, need another check to see if the nodes coincide
        gradmu1 = [0.0; 0.0]
        gradmu2 = [0.0; 0.0]

        # check if on a trailing edge
        if i ∈ TEidxs

            idx = findall(x -> x == i, TEidxs)

            if !isapprox(TEnodes[idx[1][1],1,:], TEnodes[idx[1][1],2,:])

                if i ∉ TEidxs[:,1]
                    # if !isapprox(cp[i,1]-cp[i-1,1], 0.0)
                    #     gradmu1[1] = len[i-1]*(mu[i]-mu[i-1])./(cp[i,1]-cp[i-1,1])
                    # end
                    # if !isapprox(cp[i,2]-cp[i-1,2], 0.0)
                    #     gradmu1[2] = len[i-1]*(mu[i]-mu[i-1])./(cp[i,2]-cp[i-1,2])
                    # end
                    calc_gradmu(gradmu1, len[i-1], cp[i,:], cp[i-1,:], mu[i], mu[i-1])

                    #check if last panel on a body
                elseif i ∉ TEidxs[:,2]
                    # if !isapprox(cp[i,1]-cp[i+1,1], 0.0)
                    #     gradmu2[1] = len[i+1]*(mu[i]-mu[i+1])./(cp[i,1]-cp[i+1,1])
                    # end
                    # if !isapprox(cp[i,2]-cp[i+1,2], 0.0)
                    #     gradmu2[2] = len[i+1]*(mu[i]-mu[i+1])./(cp[i,2]-cp[i+1,2])
                    # end
                    calc_gradmu(gradmu2, len[i+1], cp[i,:], cp[i+1,:], mu[i], mu[i+1])
                end
            end

        else
            calc_gradmu(gradmu1, len[i-1], cp[i,:], cp[i-1,:], mu[i], mu[i-1])
            calc_gradmu(gradmu2, len[i+1], cp[i,:], cp[i+1,:], mu[i], mu[i+1])
        end

        # get total weighting
        wt = i in TEidxs[:,1] ? len[i+1] : i in TEidxs[:,2] ? len[i-1] : len[i-1] + len[i+1]

        # calculate weighted gradient
        # TODO: this is handwavy solution to large gradients...
        # gradmu = (gradmu1+gradmu2).*(1.0 .- nhat[i,:].^2)/wt
        gradmu = (gradmu1+gradmu2).*(that[i,:].^2)/wt

        if norm(gradmu) > 1e6
            gradmu = [0.0; 0.0]
        end

        vel[:] -= gradmu/2.0

    end

    return nothing

end
