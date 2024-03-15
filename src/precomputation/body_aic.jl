#=
Body Linear System vectors and matrices
Influence matrices are given in terms of normal or tangential velocities at control points.
=#

#---------------------------------#
#     Influence Coefficients      #
#---------------------------------#

##### ----- Vortex ----- #####
"""
out of place calculation of panel method influence coefficients (V dot nhat) for a set of control points (on panels) due to a set of axisymmetric vortex rings (also on body surface)

Used for constructing the LHS influence Matrix for the panel method system, as well as RHS due to wake influences.

# Arguments:
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `node::Matrix{Float}` : [z r] coordinates of panel nodes (edges)
- `nodemap::Matrix{Int}` : [1 2] node indices for each panel
- `influence_length::Vector{Float}` : lengths of influencing panels

# Returns:
- `AICn::Matrix{Float}` : N controlpoint x N+1 node  array of V dot nhat values
"""
function vortex_aic_boundary_on_boundary(
    controlpoint, normal, node, nodemap, influence_length
)
    T = promote_type(eltype(node), eltype(controlpoint))
    M = size(controlpoint, 2)
    N = size(node, 2)

    AICn = zeros(T, M, N)

    vortex_aic_boundary_on_boundary!(
        AICn, controlpoint, normal, node, nodemap, influence_length
    )

    return AICn
end

"""
    vortex_aic_boundary_on_boundary!(
    AICn, controlpoint, normal, node, nodemap, influence_length
)

in place calculation of panel method influence coefficients (V dot nhat) for a set of control points (on panels) due to a set of axisymmetric vortex rings (also on body surface)

Used for constructing the LHS influence Matrix for the panel method system, as well as RHS due to wake influences.

# Arguments:
- `AIC::Matrix{Float}` : N-controlpoint x N-node  array of V dot nhat values for
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `node::Matrix{Float}` : [z r] coordinates of vortex rings
- `nodemap::Matrix{Int}` : [id1 id2] id numbers for node on panel used for bookkeeping due to integration across panels.
- `influence_length::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
"""
function vortex_aic_boundary_on_boundary!(
    AICn, controlpoint, normal, node, nodemap, influence_length; cache_vec=nothing
)

    # NOTE: it is slighlty faster/fewer allocations to just define a new static array in the loop than to preallocate such a small matrix.
    # vel = zeros(eltype(AICn), 2, 2)
    if isnothing(cache_vec)
        cache_vec = zeros(eltype(node), 20)
    end

    # loop through panels doing the influencing
    for (j, (nmap, lj)) in enumerate(zip(eachcol(nodemap), influence_length))
        # Loop through control points being influenced
        for (i, (cpi, nhat)) in enumerate(zip(eachcol(controlpoint), eachcol(normal)))
            n1 = view(node, :, nmap[1])
            n2 = view(node, :, nmap[2])

            # get unit induced velocity from the panel onto the control point
            if i != j
                # vel .= nominal_vortex_panel_integration(n1, n2, lj, cpi)
                vel = StaticArrays.SMatrix{2,2}(
                    nominal_vortex_panel_integration(n1, n2, lj, cpi, cache_vec)
                )
            else
                # vel .= self_vortex_panel_integration(n1, n2, lj, cpi)
                vel = StaticArrays.SMatrix{2,2}(
                    self_vortex_panel_integration(n1, n2, lj, cpi, cache_vec)
                )
            end

            for k in 1:2
                # fill the Matrix
                AICn[i, nmap[k]] += dot(vel[k, :], nhat)
            end #for k
        end #for i
    end #for j

    return nothing
end

"""
out of place calculation of panel method influence coefficients (V dot nhat) for a set of control points (NOT on panels) due to a set of axisymmetric vortex rings (on body surface)

Used for constructing the LHS influence Matrix for the panel method system, as well as RHS due to wake influences.

# Arguments:
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `node::Matrix{Float}` : [z r] coordinates of panel nodes (edges)
- `nodemap::Matrix{Int}` : [1 2] node indices for each panel
- `influence_length::Vector{Float}` : lengths of influencing panels

# Returns:
- `AICn::Matrix{Float}` : N controlpoint x N+1 node  array of V dot nhat values
"""
function vortex_aic_boundary_on_field(
    controlpoint, normal, node, nodemap, influence_length; cache_vec=nothing
)
    T = promote_type(eltype(node), eltype(controlpoint))
    M = size(controlpoint, 2)
    N = size(node, 2)

    AICn = zeros(T, M, N)

    vortex_aic_boundary_on_field!(
        AICn, controlpoint, normal, node, nodemap, influence_length; cache_vec=cache_vec
    )

    return AICn
end

"""
in place calculation of panel method influence coefficients (V dot nhat) for a set of control points (NOT on panels) due to a set of axisymmetric vortex rings (on body surface)

Used for constructing the LHS influence Matrix for the panel method system, as well as RHS due to wake influences.

# Arguments:
- `AIC::Matrix{Float}` : N-controlpoint x N-node  array of V dot nhat values for
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `node::Matrix{Float}` : [z r] coordinates of vortex rings
- `influence_length::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
"""
function vortex_aic_boundary_on_field!(
    AICn, controlpoint, normal, node, nodemap, influence_length; cache_vec=nothing
)
    # vel = zeros(eltype(AICn), 2, 2)
    if isnothing(cache_vec)
        cache_vec = zeros(eltype(node), 20)
    end

    # Loop through control points being influenced
    for (i, (cpi, nhat)) in enumerate(zip(eachcol(controlpoint), eachcol(normal)))
        # loop through panels doing the influencing
        for (j, (nmap, lj)) in enumerate(zip(eachcol(nodemap), influence_length))
            n1 = view(node, :, nmap[1])
            n2 = view(node, :, nmap[2])

            # check of self-induced:
            if isapprox(cpi, 0.5 * (n1 .+ n2))
                # if so:
                # vel .= self_vortex_panel_integration(n1, n2, lj, cpi)
                vel = StaticArrays.SMatrix{2,2}(
                    self_vortex_panel_integration(n1, n2, lj, cpi, cache_vec)
                )
            else
                # if not:
                # vel .= nominal_vortex_panel_integration(n1, n2, lj, cpi)
                vel = StaticArrays.SMatrix{2,2}(
                    nominal_vortex_panel_integration(n1, n2, lj, cpi, cache_vec)
                )
            end

            # # get unit induced velocity from the panel onto the control point
            # # vel .= nominal_vortex_panel_integration(n1, n2, lj, cpi)
            # vel = StaticArrays.SMatrix{2,2}(
            #     nominal_vortex_panel_integration(n1, n2, lj, cpi)
            # )

            for k in 1:2
                # fill the Matrix
                AICn[i, nmap[k]] += dot(vel[k, :], nhat)
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
    for kid in eachcol(kids)
        LHS[kid[1], kid[2]] = 1.0
    end

    return LHS
end

##### ----- Trailing Edge Gap Panels ----- #####
"""
"""
function add_te_gap_aic!(
    AICn,
    controlpoint,
    normal,
    tenode,
    teinfluence_length,
    tendotn,
    tencrossn,
    teadjnodeidxs;
    wake=false,
    cache_vec=nothing,
)
    if isnothing(cache_vec)
        cache_vec = zeros(eltype(controlpoint), 20)
    end

    # Loop through control points being influenced
    for (i, (cpi, nhat)) in enumerate(zip(eachcol(controlpoint), eachcol(normal)))
        # loop through bodies
        for (j, (lj, ndn, ncn, nmap)) in enumerate(
            zip(
                teinfluence_length,
                eachcol(tendotn),
                eachcol(tencrossn),
                eachcol(teadjnodeidxs),
            ),
        )

            # get unit induced velocity from the panel onto the control point
            vvel = StaticArrays.SMatrix{2,2}(
                nominal_vortex_panel_integration(
                    tenode[j, 1, :], tenode[j, 2, :], lj, cpi, cache_vec
                ),
            )
            svel = StaticArrays.SMatrix{2,2}(
                nominal_source_panel_integration(
                    tenode[j, 1, :], tenode[j, 2, :], lj, cpi, cache_vec
                ),
            )

            for k in 1:2
                # fill the Matrix
                AICn[i, nmap[k]] += dot(ndn[k] * vvel[k, :] + ncn[k] * svel[k, :], nhat)
                if wake
                    # wake "TE Panels" only have the vortex influence
                    AICn[i, nmap[k]] += dot(ndn[k] * vvel[k, :], nhat)
                end
            end #for k
        end #for j
    end #for i

    return AICn
end

# TODO: need to test these
##### ----- Source ----- #####
"""
out of place calculation of panel method influence coefficients (V dot nhat) for a set of control points (on panels) due to a set of axisymmetric source rings not on body surface.

Used for constructing the RHS boundary conditions due to rotor source panels.

# Arguments:
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `node::Matrix{Float}` : [z r] coordinates of panel nodes (edges)
- `nodemap::Matrix{Int}` : [1 2] node indices for each panel
- `influence_length::Vector{Float}` : lengths of influencing panels

# Returns:
- `AICn::Matrix{Float}` : N controlpoint x N+1 node  array of V dot nhat values
"""
function source_aic(
    controlpoint, normal, node, nodemap, influence_length; cache_vec=nothing
)
    T = promote_type(eltype(node), eltype(controlpoint))
    M = size(controlpoint, 2)
    N = size(node, 2)

    AICn = zeros(T, M, N)

    source_aic!(
        AICn, controlpoint, normal, node, nodemap, influence_length; cache_vec=cache_vec
    )

    return AICn
end

"""
    source_aic!(
    AICn, controlpoint, normal, node, nodemap, influence_length
)

# Arguments:
- `AIC::Matrix{Float}` : N-controlpoint x N-node  array of V dot nhat values for
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `node::Matrix{Float}` : [z r] coordinates of vortex rings
- `nodemap::Matrix{Int}` : [id1 id2] id numbers for node on panel used for bookkeeping due to integration across panels.
- `influence_length::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
"""
function source_aic!(
    AICn, controlpoint, normal, node, nodemap, influence_length; cache_vec=nothing
)
    # vel = zeros(eltype(AICn), 2, 2)
    if isnothing(cache_vec)
        cache_vec = zeros(eltype(node), 20)
    end

    # Loop through control points being influenced
    for (i, (cpi, nhat)) in enumerate(zip(eachcol(controlpoint), eachcol(normal)))
        # loop through panels doing the influencing
        for (j, (nmap, lj)) in enumerate(zip(eachcol(nodemap), influence_length))
            n1 = view(node, :, nmap[1])
            n2 = view(node, :, nmap[2])

            # get unit induced velocity from the panel onto the control point
            # vel .= nominal_vortex_panel_integration(n1, n2, lj, cpi)
            vel = StaticArrays.SMatrix{2,2}(
                nominal_source_panel_integration(n1, n2, lj, cpi, cache_vec)
            )

            for k in 1:2
                # fill the Matrix
                AICn[i, nmap[k]] += dot(vel[k, :], nhat)
            end #for k
        end #for j
    end #for i

    return AICn
end

##### ----- Freestream ----- #####
"""
out of place calculation of RHS contributions due to freestream.
Note that the freestream is assumed to have zero radial component in the underlying theory, but here we allow an arbitrary 2D vector for velocity for taking the dot product easier.

# Arguments:
- `normals::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `Vinfmat::Matrix{Float}` : [z r] components of freestream velocity (r's should be zero)

# Returns:
- `RHS::Vector{Float}` : vector of normal components of freestream velocity on input panels
"""
function freestream_influence_vector(normals, Vinfmat)
    T = promote_type(eltype(normals), eltype(Vinfmat))
    N = size(normals, 2)

    RHS = zeros(T, N)

    freestream_influence_vector!(RHS, normals, Vinfmat)

    return RHS
end

"""
in place calculation of RHS contributions due to freestream.
Note that the freestream is assumed to have zero radial component in the underlying theory, but here we allow an arbitrary 2D vector for velocity for taking the dot product easier.

# Arguments:
- `RHS::Vector{Float}` : vector of normal components of freestream velocity on input panels
- `normals::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `Vinfmat::Matrix{Float}` : [z r] components of freestream velocity (r's should be zero)
"""
function freestream_influence_vector!(RHS, normals, Vinfmat)
    for (i, (n, v)) in enumerate(zip(eachcol(normals), eachcol(Vinfmat)))
        RHS[i] -= dot(v, n)
    end

    return nothing
end

#---------------------------------#
#      LHS Matrix Assembly        #
#---------------------------------#
"""
"""
function assemble_lhs_matrix(
    AICn, AICpcp, npanel, nnode, totpanel, totnode, prescribednodeidxs; dummyval=1.0
)

    # get type
    TF = promote_type(eltype(AICn), eltype(AICpcp))

    # initialize LHS matrix
    LHS = zeros(TF, totnode + 2, totnode + 2)

    assemble_lhs_matrix!(
        LHS,
        AICn,
        AICpcp,
        npanel,
        nnode,
        totpanel,
        totnode,
        prescribednodeidxs;
        dummyval=dummyval,
    )

    return LHS
end

"""
AICn is nominal v dot n
AICpcp is v dot n where n is for psuedo control point
AICg is v dot n where v is from trailing edge gap panels
"""
function assemble_lhs_matrix!(
    LHS, AICn, AICpcp, npanel, nnode, totpanel, totnode, prescribednodeidxs; dummyval=1.0
)

    # - Place standard influence coefficients - #
    # 1:totpanel, 1:totnode are standard AIC values
    LHS[1:totpanel, 1:totnode] .= AICn

    # - Place Dummy influences on right-most columns - #
    # Duct Dummy's
    LHS[1:npanel[1], end - 1] .+= dummyval
    # Hub Dummy's
    LHS[(npanel[1] + 1):totpanel, end] .+= dummyval

    # - Place internal duct pseudo control point row - #
    LHS[totpanel + 1, 1:totnode] .+= @view(AICpcp[1, :])

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
        LHS[totpanel + 4, 1:totnode] .= @view(AICpcp[2, :])
    end

    return LHS
end

#---------------------------------#
#      RHS Matrix Assembly        #
#---------------------------------#
"""
"""
function assemble_rhs_matrix(
    vdnb, vdnpcp, npanel, nnode, totpanel, totnode, prescribednodeidxs
)

    # get type
    TF = promote_type(eltype(vdnb), eltype(vdnpcp))

    # initialize RHS matrix
    RHS = zeros(TF, totnode + 2)

    assemble_rhs_matrix!(
        RHS, vdnb, vdnpcp, npanel, nnode, totpanel, totnode, prescribednodeidxs
    )

    return RHS
end

"""
"""
function assemble_rhs_matrix!(
    RHS, vdnb, vdnpcp, npanel, nnode, totpanel, totnode, prescribednodeidxs
)

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

function assemble_rhs_matrix!(
    RHS, Vinf, normal, itnormal, totpanel, totnode, prescribednodeidxs
)

    # - Place standard influence coefficients - #
    # 1:totpanel, 1:totnode are standard AIC values
    RHS[1:totpanel] .= vdnb
    freestream_influence_vector!(view(RHS, 1:totpanel), normal, Vinf)

    # - Place Duct pseudo control point freestream influence - #
    vdnpcp = freestream_influence_vector(itnormal, Vinf)
    RHS[totnode - 1] = vdnpcp[1]

    # - Place hub pseudo control point freestream influence - #
    if length(prescribednodeidxs) == 1
        RHS[totnode + 2] = vdnpcp[2]
    end

    return RHS
end

"""
"""
function calculate_normal_velocity(velocity_arrray, normal)
    AIC = zeros(eltype(velocity_arrray), size(velocity_arrray, 1), size(velocity_arrray, 2))

    calculate_normal_velocity!(AIC, velocity_arrray, normal)

    return AIC
end

"""
"""
function calculate_normal_velocity!(AIC, velocity_arrray, normal)
    for j in 1:size(AIC, 2)
        for i in 1:size(AIC, 1)
            AIC[i, j] = dot(velocity_arrray[i, j, :], normal[:, i])
        end
    end
    return AIC
end
