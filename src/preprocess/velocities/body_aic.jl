#---------------------------------#
#     Influence Coefficients      #
#---------------------------------#

##### ----- Vortex ----- #####
"""
    vortex_aic_boundary_on_boundary(
        controlpoint, normal, node, nodemap, influence_length, integration_options
    )

Calculate panel method influence coefficients (V dot nhat) for a set of control points (on panels) due to a set of axisymmetric vortex rings (also on body surface)

Can be used for constructing the LHS influence Matrix for the panel method system.

# Arguments
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `node::Matrix{Float}` : [z r] coordinates of panel nodes (edges)
- `nodemap::Matrix{Int}` : [1 2] node indices for each panel
- `influence_length::Vector{Float}` : lengths of influencing panels
- `integration_options::IntegrationOptions` : integration options

# Returns
- `AICn::Matrix{Float}` : N controlpoint x N+1 node  array of V dot nhat values
"""
function vortex_aic_boundary_on_boundary(
    controlpoint, normal, node, nodemap, influence_length, integration_options
)
    T = promote_type(eltype(node), eltype(controlpoint))
    M = size(controlpoint, 2)
    N = size(node, 2)

    AICn = zeros(T, M, N)

    vortex_aic_boundary_on_boundary!(
        AICn, controlpoint, normal, node, nodemap, influence_length, integration_options
    )

    return AICn
end

"""
    vortex_aic_boundary_on_boundary!(
        AICn,
        controlpoint,
        normal,
        node,
        nodemap,
        influence_length,
        integration_options;
        integration_caches=nothing,
    )

In-place verion of `vortex_aic_boundary_on_boundary`.

integration_caches is a named tuple containing caching for intermediate calculation values.
"""
function vortex_aic_boundary_on_boundary!(
    AICn,
    controlpoint,
    normal,
    node,
    nodemap,
    influence_length,
    integration_options;
    integration_caches=nothing,
)

    # NOTE: it is slighlty faster/fewer allocations to just define a new static array in the loop than to preallocate such a small matrix.
    # vel = zeros(eltype(AICn), 2, 2)
    if isnothing(integration_caches)
        # integration_cache = zeros(eltype(controlpoint), 20)
        nominal_integration_cache = allocate_integration_containers(
            integration_options.nominal, eltype(AICn)
        )
        singular_integration_cache = allocate_integration_containers(
            integration_options.singular, eltype(AICn)
        )
    else
        nominal_integration_cache = integration_caches.nominal
        singular_integration_cache = integration_caches.singular
    end

    # loop through panels doing the influencing
    for (j, (nmap, lj)) in enumerate(zip(eachcol(nodemap), influence_length))
        # Loop through control points being influenced
        for (i, (cpi, nhat)) in enumerate(zip(eachcol(controlpoint), eachcol(normal)))
            n1 = view(node, :, Int(nmap[1]))
            n2 = view(node, :, Int(nmap[2]))

            # get unit induced velocity from the panel onto the control point
            if i != j
                # vel .= nominal_vortex_panel_integration(n1, n2, lj, cpi)
                vel = StaticArrays.SMatrix{2,2}(
                    nominal_vortex_panel_integration(
                        integration_options.nominal,
                        n1,
                        n2,
                        lj,
                        cpi,
                        nominal_integration_cache,
                    ),
                )
            else
                # vel .= self_vortex_panel_integration(n1, n2, lj, cpi)
                vel = StaticArrays.SMatrix{2,2}(
                    self_vortex_panel_integration(
                        integration_options.singular,
                        n1,
                        n2,
                        lj,
                        cpi,
                        singular_integration_cache,
                    ),
                )
            end

            for k in 1:2
                # fill the Matrix
                AICn[i, Int(nmap[k])] += dot(vel[k, :], nhat)
            end #for k
        end #for i
    end #for j

    return nothing
end

"""
    vortex_aic_boundary_on_field(
        controlpoint,
        normal,
        node,
        nodemap,
        influence_length,
        integration_options;
        integration_caches=nothing,
    )

Calculate panel method influence coefficients (V dot nhat) for a set of control points (NOT on panels) due to a set of axisymmetric vortex rings (on body surface)

Used for constructing portions of the panel method LHS matrix related to the pseudo control points in the bodies.

# Arguments
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `node::Matrix{Float}` : [z r] coordinates of panel nodes (edges)
- `nodemap::Matrix{Int}` : [1 2] node indices for each panel
- `influence_length::Vector{Float}` : lengths of influencing panels
- `integration_options::IntegrationOptions` : integration options

# Keyword Arguments
- `integration_caches::NamedTuple=nothing` : caches for intermediate values in integration.

# Returns
- `AICn::Matrix{Float}` : N controlpoint x N+1 node  array of V dot nhat values
"""
function vortex_aic_boundary_on_field(
    controlpoint,
    normal,
    node,
    nodemap,
    influence_length,
    integration_options;
    integration_caches=nothing,
)
    T = promote_type(eltype(node), eltype(controlpoint))
    M = size(controlpoint, 2)
    N = size(node, 2)

    AICn = zeros(T, M, N)

    vortex_aic_boundary_on_field!(
        AICn,
        controlpoint,
        normal,
        node,
        nodemap,
        influence_length,
        integration_options;
        integration_caches=integration_caches,
    )

    return AICn
end

"""
    vortex_aic_boundary_on_field!(
        AICn,
        controlpoint,
        normal,
        node,
        nodemap,
        influence_length,
        integration_options;
        integration_caches=nothing,
    )

In-place version of `vortex_aic_boundary_on_field`.
"""
function vortex_aic_boundary_on_field!(
    AICn,
    controlpoint,
    normal,
    node,
    nodemap,
    influence_length,
    integration_options;
    integration_caches=nothing,
)

    # vel = zeros(eltype(AICn), 2, 2)
    if isnothing(integration_caches)
        # integration_cache = zeros(eltype(controlpoint), 20)
        nominal_integration_cache = allocate_integration_containers(
            integration_options.nominal, eltype(AICn)
        )
        singular_integration_cache = allocate_integration_containers(
            integration_options.singular, eltype(AICn)
        )
    else
        nominal_integration_cache = integration_caches.nominal
        singular_integration_cache = integration_caches.singular
    end
    # Loop through control points being influenced
    for (i, (cpi, nhat)) in enumerate(zip(eachcol(controlpoint), eachcol(normal)))
        # loop through panels doing the influencing
        for (j, (nmap, lj)) in enumerate(zip(eachcol(nodemap), influence_length))
            n1 = view(node, :, Int(nmap[1]))
            n2 = view(node, :, Int(nmap[2]))

            # check of self-induced:
            if isapprox(cpi, 0.5 * (n1 .+ n2))
                # if so:
                # vel .= self_vortex_panel_integration(n1, n2, lj, cpi)
                vel = StaticArrays.SMatrix{2,2}(
                    self_vortex_panel_integration(
                        integration_options.singular,
                        n1,
                        n2,
                        lj,
                        cpi,
                        singular_integration_cache,
                    ),
                )
            else
                # if not:
                # vel .= nominal_vortex_panel_integration(n1, n2, lj, cpi)
                vel = StaticArrays.SMatrix{2,2}(
                    nominal_vortex_panel_integration(
                        integration_options.nominal,
                        n1,
                        n2,
                        lj,
                        cpi,
                        nominal_integration_cache,
                    ),
                )
            end

            # # get unit induced velocity from the panel onto the control point
            # # vel .= nominal_vortex_panel_integration(n1, n2, lj, cpi)
            # vel = StaticArrays.SMatrix{2,2}(
            #     nominal_vortex_panel_integration(n1, n2, lj, cpi)
            # )

            for k in 1:2
                # fill the Matrix
                AICn[i, Int(nmap[k])] += dot(vel[k, :], nhat)
            end #for k
        end #for j
    end #for i

    return nothing
end

##### ----- Add Kutta ----- #####
"""
    add_kutta!(LHS, AICn, kids)

Add Kutta condition (γ_1 + γ_N = 0) to LHS matrix.

# Arguments
- `LHS::Matrix{Float}` : a pre-allocated (zeros) full size left-hand side matrix
- `AICn::Matrix{Float}` :  influence coefficients for panels/nodes
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
    add_te_gap_aic!(
        AICn,
        controlpoint,
        normal,
        tenode,
        teinfluence_length,
        tendotn,
        tencrossn,
        teadjnodeidxs,
        integration_options;
        wake=false,
        integration_caches=nothing,
    )

Add trailing edge gap aerodynmic influence coefficient contributions to the AIC matrix.

# Arguments
- `AICn::Matrix{Float}` : N controlpoint x N+1 node  array of V dot nhat values
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `tenode::Matrix{Float}` : [z r] coordinates of trailing edge panel nodes (edges)
- `teinfluence_length::Vector{Float}` : lengths of influencing trailing edge panels
- `tendotn::Matrix{Float}` : nhat of trailing edge panel dotted with nhat of adjacent panels
- `tencrossn::Matrix{Float}` : nhat of trailing edge panel crossed with nhat of adjacent panels
- `teadjnodeidxs::Matrix{Float}` : indices of nodes adjacent to trailing edge panel
- `integration_options::IntegrationOptions` : integration options

# Keyword Arguments
- `wake::Bool=false` : flag as to whether this function is being applied to a wake sheet.
- `integration_caches::NamedTuple=nothing` : caches for intermediate values in integration.
"""
function add_te_gap_aic!(
    AICn,
    controlpoint,
    normal,
    tenode,
    teinfluence_length,
    tendotn,
    tencrossn,
    teadjnodeidxs,
    integration_options;
    wake=false,
    integration_caches=nothing,
)
    if isnothing(integration_caches)
        # integration_cache = zeros(eltype(controlpoint), 20)
        nominal_integration_cache = allocate_integration_containers(
            integration_options.nominal, eltype(AICn)
        )
    else
        nominal_integration_cache = integration_caches
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
                    integration_options.nominal,
                    tenode[j, 1, :],
                    tenode[j, 2, :],
                    lj,
                    cpi,
                    nominal_integration_cache,
                ),
            )
            svel = StaticArrays.SMatrix{2,2}(
                nominal_source_panel_integration(
                    integration_options.nominal,
                    tenode[j, 1, :],
                    tenode[j, 2, :],
                    lj,
                    cpi,
                    nominal_integration_cache,
                ),
            )

            for k in 1:2
                # fill the Matrix
                AICn[i, Int(nmap[k])] += dot(
                    ndn[k] * vvel[k, :] + ncn[k] * svel[k, :], nhat
                )
                #TODO: is this a bug? seems like the above line should be in an else statement
                if wake
                    # wake "TE Panels" only have the vortex influence
                    AICn[i, Int(nmap[k])] += dot(ndn[k] * vvel[k, :], nhat)
                end
            end #for k
        end #for j
    end #for i

    return AICn
end

##### ----- Source ----- #####
"""
    source_aic(
        controlpoint,
        normal,
        node,
        nodemap,
        influence_length,
        integration_options;
        integration_caches=nothing,
    )

Calculate panel method influence coefficients (V dot nhat) for a set of control points (on panels) due to a set of axisymmetric source rings not on body surface.

Can be used for constructing the RHS boundary conditions due to rotor source panels.

# Arguments
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `normal::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `node::Matrix{Float}` : [z r] coordinates of panel nodes (edges)
- `nodemap::Matrix{Int}` : [1 2] node indices for each panel
- `influence_length::Vector{Float}` : lengths of influencing panels
- `integration_options::IntegrationOptions` : integration options

# Keyword Arguments
- `integration_caches::NamedTuple=nothing` : caches for intermediate values in integration.

# Returns
- `AICn::Matrix{Float}` : N controlpoint x N+1 node  array of V dot nhat values
"""
function source_aic(
    controlpoint,
    normal,
    node,
    nodemap,
    influence_length,
    integration_options;
    integration_caches=nothing,
)
    T = promote_type(eltype(node), eltype(controlpoint))
    M = size(controlpoint, 2)
    N = size(node, 2)

    AICn = zeros(T, M, N)

    source_aic!(
        AICn,
        controlpoint,
        normal,
        node,
        nodemap,
        influence_length,
        integration_options;
        integration_caches=integration_caches,
    )

    return AICn
end

"""
    source_aic!(
        AICn,
        controlpoint,
        normal,
        node,
        nodemap,
        influence_length,
        integration_options;
        integration_caches=nothing,
    )

In-place version of `source_aic`.
"""
function source_aic!(
    AICn,
    controlpoint,
    normal,
    node,
    nodemap,
    influence_length,
    integration_options;
    integration_caches=nothing,
)
    # vel = zeros(eltype(AICn), 2, 2)
    if isnothing(integration_caches)
        # integration_cache = zeros(eltype(controlpoint), 20)
        nominal_integration_cache = allocate_integration_containers(
            integration_options.nominal, eltype(AICn)
        )
    else
        nominal_integration_cache = integration_caches
    end

    # Loop through control points being influenced
    for (i, (cpi, nhat)) in enumerate(zip(eachcol(controlpoint), eachcol(normal)))
        # loop through panels doing the influencing
        for (j, (nmap, lj)) in enumerate(zip(eachcol(nodemap), influence_length))
            n1 = view(node, :, Int(nmap[1]))
            n2 = view(node, :, Int(nmap[2]))

            # get unit induced velocity from the panel onto the control point
            # vel .= nominal_vortex_panel_integration(n1, n2, lj, cpi)
            vel = StaticArrays.SMatrix{2,2}(
                nominal_source_panel_integration(
                    integration_options.nominal, n1, n2, lj, cpi, nominal_integration_cache
                ),
            )

            for k in 1:2
                # fill the Matrix
                AICn[i, Int(nmap[k])] += dot(vel[k, :], nhat)
            end #for k
        end #for j
    end #for i

    return AICn
end

##### ----- Freestream ----- #####
"""
    freestream_influence_vector(normals, Vinfmat)

Calculate RHS contributions due to freestream.

Note that the freestream is assumed to have zero radial component in the underlying theory, but here we allow an arbitrary 2D vector for velocity for taking the dot product easier.

# Arguments
- `normals::Matrix{Float}` : unit normal vectors of the panels on which the control points are situated.
- `Vinfmat::Matrix{Float}` : [z r] components of freestream velocity (r's should be zero)

# Returns
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
    freestream_influence_vector!(RHS, normals, Vinfmat)

In-place version of `freestream_influence_vector`.
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
    assemble_lhs_matrix(
        AICn, AICpcp, npanel, nnode, totpanel, totnode, prescribednodeidxs; dummyval=1.0
    )

Assemble the LHS matrix of the panel method.

# Arguments
- `AICn::Matrix{Float}` : N controlpoint x N+1 node  array of V dot nhat values
- `AICpcp::Matrix{Float}` : Nbodies controlpoint x N+1 node  array of V dot nhat values (influence on psuedo control points)
- `npanel::Vector{Int}` : number of panels comprising each body
- `nnode::Vector{Int}` : number of nodes comprising each body
- `totpanel::Int` : total number of panels
- `totnode::Int` : total number of nodes
- `prescribednodeidxs::Vector{Int}` : indices of nodes with prescribed strengths (those on the axis of rotation).

# Keyword Arguments
- `dummyval::Float=1.0` : value for dummy input for prescribed and internal control points in the system. Do not change except for debugging purposes.

# Returns
- `LHS::Matrix{Float}` : The full LHS matrix for the panel method.
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
    assemble_lhs_matrix!(
        LHS, AICn, AICpcp, npanel, nnode, totpanel, totnode, prescribednodeidxs; dummyval=1.0
    )

In-place version of `assemble_lhs_matrix`.
"""
function assemble_lhs_matrix!(
    LHS, AICn, AICpcp, npanel, nnode, totpanel, totnode, prescribednodeidxs; dummyval=1.0
)

    # - Place standard influence coefficients - #
    # 1:totpanel, 1:totnode are standard AIC values
    LHS[1:Int(totpanel[]), 1:Int(totnode[])] .= AICn

    # - Place Dummy influences on right-most columns - #
    # Duct Dummy's
    LHS[1:Int(npanel[1]), end - 1] .= dummyval
    # Hub Dummy's
    LHS[(Int(npanel[1]) + 1):Int(totpanel[]), end] .= dummyval

    # - Place internal duct pseudo control point row - #
    LHS[Int(totpanel[]) + 1, 1:Int(totnode[])] .= @view(AICpcp[1, :])

    # - Place Kutta condition Row - #
    LHS[(Int(totpanel[]) + 2), 1] = LHS[Int(totpanel[]) + 2, Int(nnode[1])] = 1.0

    # - Place hub LE prescribed panel row - #
    LHS[Int(totpanel[]) + 3, Int(prescribednodeidxs[1])] = 1.0

    # - Place hub TE prescribed panel row OR hub internal pseudo control point row - #
    if iszero(Int(prescribednodeidxs[2]))
        # internal pseudo control point
        LHS[Int(totpanel[]) + 4, 1:Int(totnode[])] .= @view(AICpcp[2, :])
    else
        # prescribed TE node
        LHS[Int(totpanel[]) + 4, Int(prescribednodeidxs[2])] = 1.0
    end

    return LHS
end

"""
    factorize_LHS(A::AbstractMatrix)

Returns the LU decomposition of `A`.
"""
function factorize_LHS(A::AbstractMatrix{T}) where {T}

    # Allocate memory for pivot
    if T <: ForwardDiff.Dual #|| T<:RD.TrackedReal  # Automatic differentiation case
        Tprimal = T.parameters[T <: ForwardDiff.Dual ? 2 : 1]
        Apivot = zeros(Tprimal, size(A))

    else
        Apivot = zeros(T, size(A))
    end

    # LU decomposition
    return factorize_LHS!(Apivot, A)
end

"""
    factorize_LHS!(Apivot::AbstractMatrix, A::AbstractMatrix)

Returns the LU decomposition of `A` using `Apivot` as storage memory to pivot
leaving `A` unchanged.
"""
function factorize_LHS!(Apivot, A::AbstractMatrix{T}) where {T}

    # Prepare pivot array
    extract_primals!(Apivot, A)

    # LU decomposition
    Alu = LinearAlgebra.lu!(Apivot, NoPivot(); check=false)

    return Alu
end

#---------------------------------#
#      RHS Matrix Assembly        #
#---------------------------------#
"""
    assemble_rhs_matrix(
        vdnb, vdnpcp, npanel, nnode, totpanel, totnode, prescribednodeidxs
    )

# Arguments
- `vdnb::Vector{Float}` : V dot nhat for body panels
- `vdnpcp::Vector{Float}` : V dot nhat for pseudo control points
- `npanel::Vector{Int}` : number of panels comprising each body
- `nnode::Vector{Int}` : number of nodes comprising each body
- `totpanel::Int` : total number of body panels
- `totnode::Int` : total number of body nodes
- `prescribednodeidxs::Vector{Int}` : indices of nodes with prescribed strengths (those on the axis of rotation)

# Returns
- `RHS::Vector{Float}` : the RHS vector of the panel method.
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
    assemble_rhs_matrix!(
        RHS, vdnb, vdnpcp, npanel, nnode, totpanel, totnode, prescribednodeidxs
    )

In-place version of `assemble_rhs_matrix`.
"""
function assemble_rhs_matrix!(
    RHS, vdnb, vdnpcp, npanel, nnode, totpanel, totnode, prescribednodeidxs
)

    # - Place standard influence coefficients - #
    # 1:totpanel, 1:totnode are standard AIC values
    RHS[1:Int(totpanel[])] .= vdnb

    # - Place Duct pseudo control point freestream influence - #
    RHS[Int(totnode[]) - 1] = vdnpcp[1]

    # - Place hub pseudo control point freestream influence - #
    if iszero(Int(prescribednodeidxs[2]))
        RHS[Int(totnode[]) + 2] = vdnpcp[2]
    end

    return RHS
end

"""
    calculate_normal_velocity(velocity_vector, normal)

Calculate normal velocity_vector (V dot nhat).

# Arguments
- `velocity_vector::Matrix{Float}` : velocity vector [z r] on each panel
- `normal::Matrix{Float}` : the panel unit normals

# Returns
- `AIC::Matrix{Float}` : V dot n on each panel
"""
function calculate_normal_velocity(velocity_vector, normal)

    # initialize
    AIC = zeros(eltype(velocity_vector), size(velocity_vector, 1), size(velocity_vector, 2))

    calculate_normal_velocity!(AIC, velocity_vector, normal)

    return AIC
end

"""
    calculate_normal_velocity!(AIC, velocity_vector, normal)

In-place version of `calculate_normal_velocity`.
"""
function calculate_normal_velocity!(AIC, velocity_vector, normal)
    for j in 1:size(AIC, 2)
        for i in 1:size(AIC, 1)
            AIC[i, j] = dot(velocity_vector[i, j, :], normal[:, i])
        end
    end
    return AIC
end
