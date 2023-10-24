#---------------------------------#
#     Influence Coefficients      #
#---------------------------------#

##### ----- Vortex ----- #####
"""
out of place calculation of panel method influence coefficients (V dot nhat) for a set of control points (on panels) due to a set of axisymmetric vortex rings (also on body surface)

Used for constructing the LHS influence Matrix for the panel method system, as well as RHS due to wake influences.

**Arguments:**
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `node::Matrix{Float}` : [z r] coordinates of panel nodes (edges)
- `nodemap::Matrix{Int}` : [1 2] node indices for each panel
- `influence_length::Vector{Float}` : lengths of influencing panels

**Returns:**
- `AICn::Matrix{Float}` : N controlpoint x N+1 node  array of V dot nhat values
- `AICt::Matrix{Float}` : N controlpoint x N+1 node  array of V dot that values
"""
function vortex_aic_boundary_on_boundary(
    controlpoint, normal, tangent, node, nodemap, influence_length
)
    T = promote_type(eltype(node), eltype(controlpoint))
    M = size(controlpoint, 1)
    N = size(node, 1)

    AICn = zeros(T, M, N)
    AICt = zeros(T, M, N)

    vortex_aic_boundary_on_boundary!(
        AICn, AICt, controlpoint, normal, tangent, node, nodemap, influence_length
    )

    return AICn, AICt
end

"""
in place calculation of panel method influence coefficients (V dot nhat) for a set of control points (on panels) due to a set of axisymmetric vortex rings (also on body surface)

Used for constructing the LHS influence Matrix for the panel method system, as well as RHS due to wake influences.

**Arguments:**
- `AIC::Matrix{Float}` : N-controlpoint x N-node  array of V dot nhat values for
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `node::Matrix{Float}` : [z r] coordinates of vortex rings
- `influence_length::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
"""
function vortex_aic_boundary_on_boundary!(
    AICn, AICt, controlpoint, normal, tangent, node, nodemap, influence_length
)

    # Loop through control points being influenced
    for (i, (cpi, nhat, that)) in
        enumerate(zip(eachrow(controlpoint), eachrow(normal), eachrow(tangent)))
        # loop through panels doing the influencing
        for (j, (nmap, lj)) in enumerate(zip(eachrow(nodemap), influence_length))

            # get unit induced velocity from the panel onto the control point
            # TODO: this allocates the vel's put 2x2 vel object in eventual cache and update the integration functions to be in place.
            # println("i:", i, " j:", j)
            if i != j
                # @time "nominal" begin
                vel = nominal_vortex_panel_integration(
                    node[nmap[1], :], node[nmap[2], :], lj, cpi
                )
                # end
            else
                # @time "self" begin
                vel = self_vortex_panel_integration(
                    node[nmap[1], :], node[nmap[2], :], lj, cpi
                )
                # end
            end

            for k in 1:2
                # fill the Matrix
                AICn[i, nmap[k]] += dot(vel[k, :], nhat)
                AICt[i, nmap[k]] += dot(vel[k, :], that)
            end #for k
        end #for j
    end #for i

    return nothing
end

"""
out of place calculation of panel method influence coefficients (V dot nhat) for a set of control points (NOT on panels) due to a set of axisymmetric vortex rings (on body surface)

Used for constructing the LHS influence Matrix for the panel method system, as well as RHS due to wake influences.

**Arguments:**
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `node::Matrix{Float}` : [z r] coordinates of panel nodes (edges)
- `nodemap::Matrix{Int}` : [1 2] node indices for each panel
- `influence_length::Vector{Float}` : lengths of influencing panels

**Returns:**
- `AICn::Matrix{Float}` : N controlpoint x N+1 node  array of V dot nhat values
- `AICt::Matrix{Float}` : N controlpoint x N+1 node  array of V dot that values
"""
function vortex_aic_boundary_on_field(
    controlpoint, normal, tangent, node, nodemap, influence_length
)
    T = promote_type(eltype(node), eltype(controlpoint))
    M = size(controlpoint, 1)
    N = size(node, 1)

    AICn = zeros(T, M, N)
    AICt = zeros(T, M, N)

    vortex_aic_boundary_on_field!(
        AICn, AICt, controlpoint, normal, tangent, node, nodemap, influence_length
    )

    return AICn, AICt
end

"""
in place calculation of panel method influence coefficients (V dot nhat) for a set of control points (NOT on panels) due to a set of axisymmetric vortex rings (on body surface)

Used for constructing the LHS influence Matrix for the panel method system, as well as RHS due to wake influences.

**Arguments:**
- `AIC::Matrix{Float}` : N-controlpoint x N-node  array of V dot nhat values for
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `node::Matrix{Float}` : [z r] coordinates of vortex rings
- `influence_length::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
"""
function vortex_aic_boundary_on_field!(
    AICn, AICt, controlpoint, normal, tangent, node, nodemap, influence_length
)

    # Loop through control points being influenced
    for (i, (cpi, nhat, that)) in
        enumerate(zip(eachrow(controlpoint), eachrow(normal), eachrow(tangent)))
        # loop through panels doing the influencing
        for (j, (nmap, lj)) in enumerate(zip(eachrow(nodemap), influence_length))

            # get unit induced velocity from the panel onto the control point
            # TODO: this allocates the vel's put a single 2x2 vel object in eventual cache and update the integration functions to be in place.
            vel = nominal_vortex_panel_integration(
                node[nmap[1], :], node[nmap[2], :], lj, cpi
            )

            for k in 1:2
                # fill the Matrix
                AICn[i, nmap[k]] += dot(vel[k, :], nhat)
                AICt[i, nmap[k]] += dot(vel[k, :], that)
            end #for k
        end #for j
    end #for i

    return nothing
end

##### ----- Add Kutta ----- #####
"""
LHS is pre-allocated (zeros) full size left-hand side matrix
AICn are computed influence coefficients for panels/nodes
- `kids::Vector{Int}` : [1 2] indices of where to put 1's for kutta condition
"""
function add_kutta!(LHS, AICn, kids)

    # zero out LHS just in case
    LHS .= 0.0

    # get size of non-kutta matrix
    ni, nj = size(AICn)

    # plug in non-kutta into kutta matrix
    LHS[1:ni, 1:nj] = AICn

    # add kutta condition
    for kid in eachrow(kids)
        LHS[kid[1], kid[2]] = 1.0
    end

    return LHS
end

##### ----- Trailing Edge Gap Panels ----- #####
"""
"""
function add_te_gap_aic!(
    AICn,
    AICt,
    controlpoint,
    normal,
    tangent,
    tenode,
    teinfluence_length,
    tendotn,
    tencrossn,
    teadjnodeidxs,
)

    # Loop through control points being influenced
    for (i, (cpi, nhat, that)) in
        enumerate(zip(eachrow(controlpoint), eachrow(normal), eachrow(tangent)))
        # loop through bodies
        for (j, (lj, ndn, ncn, nmap)) in enumerate(
            zip(
                teinfluence_length,
                eachrow(tendotn),
                eachrow(tencrossn),
                eachrow(teadjnodeidxs),
            ),
        )

            # get unit induced velocity from the panel onto the control point
            # TODO: this allocates the vel's put a single 2x2 vel object in eventual cache and update the integration functions to be in place.
            vvel = nominal_vortex_panel_integration(
                tenode[j, 1, :], tenode[j, 2, :], lj, cpi
            )
            svel = nominal_source_panel_integration(
                tenode[j, 1, :], tenode[j, 2, :], lj, cpi
            )

            for k in 1:2
                # fill the Matrix
                AICn[i, nmap[k]] += dot(ndn[k] * vvel[k, :] + ncn[k] * svel[k, :], nhat)
                # AICn[i, nmap[k]] += dot(-ncn[k] * svel[k, :], nhat)
                AICt[i, nmap[k]] += dot(ndn[k] * vvel[k, :] + ncn[k] * svel[k, :], that)
                # AICt[i, nmap[k]] += dot(-ncn[k] * svel[k, :], that)
            end #for k
        end #for j
    end #for i

    return AICn
end

#TODO: need to update these maybe, if moving to linear sources.
##### ----- Source ----- #####
"""
out of place calculation of panel method influence coefficients (V dot nhat) for a set of control points (on panels) due to a set of axisymmetric source rings

Used for constructing the RHS influence Matrix for the panel method system (rotor source induced velocity on boundaries)

**Arguments:**
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `node::Matrix{Float}` : [z r] coordinates of vortex rings
- `influence_length::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.

**Returns:**
- `AIC::Matrix{Float}` : N-controlpoint x N-node  array of V dot nhat values for
"""
function source_influence_matrix(controlpoint, normal, node, influence_length)
    T = promote_type(eltype(node), eltype(controlpoint))
    M = size(controlpoint, 1)
    N = size(node, 1)

    AIC = zeros(T, M, N)

    source_influence_matrix!(AIC, controlpoint, normal, node, influence_length)

    return AIC
end

"""
in place calculation of panel method influence coefficients (V dot nhat) for a set of control points (on panels) due to a set of axisymmetric source rings

Used for constructing the RHS influence Matrix for the panel method system (rotor source induced velocity on boundaries)

**Arguments:**
- `AIC::Matrix{Float}` : N-controlpoint x N-node  array of V dot nhat values for
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `node::Matrix{Float}` : [z r] coordinates of vortex rings
- `influence_length::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
"""
function source_influence_matrix!(AIC, controlpoint, normal, node, influence_length)

    # Loop through control points being influenced
    for (i, (cpi, nhat)) in enumerate(zip(eachrow(controlpoint), eachrow(normal)))
        # loop through panels doing the influencing
        for (j, (nj, lj)) in enumerate(zip(eachrow(node), influence_length))

            # get unit induced velocity from the panel onto the control point
            # TODO: this allocates
            vel = source_induced_velocity(nj, lj, cpi) #note input strength is default unity

            # fill the Matrix
            AIC[i, j] += dot(vel, nhat)
        end
    end

    return nothing
end

##### ----- Freestream ----- #####
"""
out of place calculation of RHS contributions due to freestream.
Note that the freestream is assumed to have zero radial component in the underlying theory, but here we allow an arbitrary 2D vector for velocity for taking the dot product easier.

**Arguments:**
- `normals::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `Vinfmat::Matrix{Float}` : [z r] components of freestream velocity (r's should be zero)

**Returns:**
- `RHS::Vector{Float}` : vector of normal components of freestream velocity on input panels
"""
function freestream_influence_vector(normals, Vinfmat)
    T = promote_type(eltype(normals), eltype(Vinfmat))
    N = size(normals, 1)

    RHS = zeros(T, N)

    freestream_influence_vector!(RHS, normals, Vinfmat)

    return RHS
end

"""
in place calculation of RHS contributions due to freestream.
Note that the freestream is assumed to have zero radial component in the underlying theory, but here we allow an arbitrary 2D vector for velocity for taking the dot product easier.

**Arguments:**
- `RHS::Vector{Float}` : vector of normal components of freestream velocity on input panels
- `normals::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `Vinfmat::Matrix{Float}` : [z r] components of freestream velocity (r's should be zero)
"""
function freestream_influence_vector!(RHS, normals, Vinfmat)
    for (i, (n, v)) in enumerate(zip(eachrow(normals), eachrow(Vinfmat)))
        RHS[i] -= dot(v, n)
    end

    return nothing
end

#---------------------------------#
#      LHS Matrix Assembly        #
#---------------------------------#
"""
"""
function assemble_lhs_matrix(AICn, AICpcp, panels; dummyval=1.0)

    # get type
    TF = promote_type(eltype(AICn), eltype(AICpcp))

    # initialize LHS matrix
    LHS = zeros(TF, panels.totnode + 2, panels.totnode + 2)

    assemble_lhs_matrix!(LHS, AICn, AICpcp, panels; dummyval=dummyval)

    return LHS
end

"""
AICn is nominal v dot n
AICpcp is v dot n where n is for psuedo control point
AICg is v dot n where v is from trailing edge gap panels
"""
function assemble_lhs_matrix!(LHS, AICn, AICpcp, panels; dummyval=1.0)

    # extract stuff from panels
    (; npanel, nnode, totpanel, totnode, prescribednodeidxs) = panels

    # - Place standard influence coefficients - #
    # 1:totpanel, 1:totnode are standard AIC values
    LHS[1:totpanel, 1:totnode] .= AICn

    # - Place Dummy influences on right-most columns - #
    # Duct Dummy's
    LHS[1:npanel[1], end - 1] .+= dummyval
    # Hub Dummy's
    LHS[(npanel[1] + 1):totpanel, end] .+= dummyval

    # - Place internal duct pseudo control point row - #
    LHS[totpanel + 1, 1:totnode] .+= AICpcp[1, :]

    # - Place Kutta condition Row - #
    LHS[(totpanel + 2), 1] = LHS[totpanel + 2, nnode[1]] = 1.0

    # - Place hub LE prescribed panel row - #
    LHS[totpanel + 3, prescribednodeidxs[1]] = 1.0

    # - Place hub TE prescribed panel row OR hub internal pseudo control point row - #
    if length(prescribednodeidxs) > 1
        # prescribed TE node
        LHS[totpanel + 4, prescribednodeidxs[2]] = 1.0
    else
        # internal pseudo control point
        LHS[totpanel + 4, 1:totnode] .= AICpcp[2, :]
    end

    return LHS
end

#---------------------------------#
#      RHS Matrix Assembly        #
#---------------------------------#
"""
"""
function assemble_rhs_matrix(vdnb, vdnpcp, panels)

    # get type
    TF = promote_type(eltype(vdnb), eltype(vdnpcp))

    # initialize RHS matrix
    RHS = zeros(TF, panels.totnode + 2)

    assemble_rhs_matrix!(RHS, vdnb, vdnpcp, panels)

    return RHS
end

"""
"""
function assemble_rhs_matrix!(RHS, vdnb, vdnpcp, panels)

    # extract stuff from panels
    (; npanel, nnode, totpanel, totnode, prescribednodeidxs) = panels

    # - Place standard influence coefficients - #
    # 1:totpanel, 1:totnode are standard AIC values
    RHS[1:totpanel] .= vdnb

    # - Place Duct pseudo control point freestream influence - #
    RHS[totnode - 1] = vdnpcp[1]

    # - Place hub pseudo control point freestream influence - #
    if length(prescribednodeidxs) == 1
        RHS[totnode + 2] = vdnpcp[2]
    end

    return RHS
end

######################################################################
#                                                                    #
#               Augment with influence on Internal Panels            #
#                                                                    #
######################################################################

###### ----- Vortex ----- #####

#function vortex_influence_on_internal_panels(RHS, itpanels, influencepanels)
#    M, N = size(RHS)
#    nit = length(itpanels.itcontrolpoint[:, 1])

#    augRHS = zeros(eltype(RHS), M + nit, N)
#    augRHS[1:M, 1:N] .= RHS

#    vortex_influence_on_internal_panels!(augRHS, itpanels, influencepanels)

#    return augRHS
#end

#function vortex_influence_on_internal_panels!(augRHS, itpanels, influencepanels)
#    M = itpanels.totpanel

#    for (j, (cp, lj)) in
#        enumerate(zip(eachrow(influencepanels.node), influencepanels.influence_length))
#        for (i, (cpi, nhat)) in
#            enumerate(zip(eachrow(itpanels.node), eachrow(itpanels.itnormal)))

#            # get unit induced velocity from the panel onto the control point
#            vel = vortex_induced_velocity(cp, lj, cpi)

#            #TODO: double check if this should be + or - or if it matters
#            augRHS[M + i, j] += dot(vel, nhat)
#        end
#    end

#    return nothing
#end

###### ----- Source ----- #####

#function source_influence_on_internal_panels(RHS, itpanels, influencepanels)
#    M, N = size(RHS)
#    nit = length(itpanels.itcontrolpoint[:, 1])

#    augRHS = zeros(eltype(RHS), M + nit, N)
#    augRHS[1:M, 1:N] .= RHS

#    source_influence_on_internal_panels!(augRHS, itpanels, influencepanels)

#    return augRHS
#end

#function source_influence_on_internal_panels!(augRHS, itpanels, influencepanels)
#    M = itpanels.totpanel

#    for (j, (cp, lj)) in
#        enumerate(zip(eachrow(influencepanels.node), influencepanels.influence_length))
#        for (i, (cpi, nhat)) in
#            enumerate(zip(eachrow(itpanels.node), eachrow(itpanels.itnormal)))

#            # get unit induced velocity from the panel onto the control point
#            vel = source_induced_velocity(cp, lj, cpi)

#            #TODO: double check if this should be + or - or if it matters
#            augRHS[M + i, j] += dot(vel, nhat)
#        end
#    end

#    return nothing
#end

###### ----- Doublet ----- #####

#function doublet_panel_influence_on_internal_panels(
#    LHS,
#    itpanels,
#    influencepanels, #prescribedpanels
#)
#    M, N = size(LHS)
#    nit = length(itpanels.itcontrolpoint[:, 1])
#    # pid = (p -> p[1]).(prescribedpanels)

#    augLHS = zeros(eltype(LHS), M + nit, N + nit)
#    augLHS[1:M, 1:N] .= LHS

#    for (i, eid) in enumerate(eachrow(itpanels.endpointidxs))
#        augLHS[eid[1]:eid[2], N + i] .= -1.0
#        # for ip in pid
#        #     augLHS[ip, N + i] = 0.0
#        # end
#    end

#    doublet_panel_influence_on_internal_panels!(augLHS, influencepanels, itpanels)

#    return augLHS
#end

#function doublet_panel_influence_on_internal_panels!(augLHS, influencepanels, itpanels)

#    # rename for convenience
#    nodes = influencepanels.nodes
#    M = itpanels.totpanel

#    for (i, (cpi, nhat)) in
#        enumerate(zip(eachrow(itpanels.itcontrolpoint), eachrow(itpanels.itnormal)))
#        for (j, (p1, p2)) in
#            enumerate(zip(eachrow(nodes[:, 1, :]), eachrow(nodes[:, 2, :])))

#            # get unit induced velocity from the panel onto the control point
#            vel = constant_doublet_band_induced_velocity(p1, p2, cpi)

#            # fill the Matrix
#            augLHS[M + i, j] += dot(vel, nhat)
#        end
#    end

#    return nothing
#end

###### ----- Freestream ----- #####

#function freestream_influence_on_internal_panels(RHS, panels, Vinfvec)
#    M = length(RHS)
#    nit = length(panels.itcontrolpoint[:, 1])

#    augRHS = zeros(eltype(RHS), M + nit)
#    augRHS[1:M] .= RHS

#    freestream_influence_on_internal_panels!(augRHS, panels, Vinfvec)

#    return augRHS
#end

#function freestream_influence_on_internal_panels!(augRHS, panels, Vinfvec)
#    M = panels.totpanel

#    for (i, nhat) in enumerate(eachrow(panels.itnormal))
#        augRHS[M + i] -= dot(Vinfvec, nhat)
#    end

#    return nothing
#end
