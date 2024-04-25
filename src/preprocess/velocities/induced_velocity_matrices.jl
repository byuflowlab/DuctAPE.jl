######################################################################
#                                                                    #
#          Matrix of velocities induced by panels on points          #
#                                                                    #
######################################################################

#---------------------------------#
#     Vortex Panels on Points     #
#---------------------------------#

"""
    induced_velocities_from_vortex_panels_on_points(
        controlpoints,
        nodes,
        nodemap,
        influence_lengths,
        strengths,
        integration_options;
        integration_caches=nothing,
    )

Calculate axial and radial components of induced velocity for a set of control points due to a set of axisymmetric vortex panels (bands).

Used for getting the unit induced velocities due to the body panels on the rotor/wake as well as the unit induced velocity due to the wake on the body/rotor.

# Arguments
- `controlpoints::Matrix{Float}` [z r] coordinates of points being influenced
- `nodes::Matrix{Float}` : [z r] coordinates of vortex rings
- `nodemap::Matrix{Int}` : mapping from panel index to associated node indices
- `influence_lengths::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
- `strengths::Matrix{Float}` : vortex constant circulation values
- `integration_options::IntegrationOptions` : integration options

# Keyword Arguments
- `integration_caches::NamedTuple=nothing` : cache used in in-place quadrature functions.

# Returns
- `VEL::Array{Float}` : N-controlpoint x N-node x [vz, vr] array of induced velocity components
"""
function induced_velocities_from_vortex_panels_on_points(
    controlpoints,
    nodes,
    nodemap,
    influence_lengths,
    strengths,
    integration_options;
    integration_caches=nothing,
)

    # Initialize
    T = promote_type(eltype(controlpoints), eltype(nodes), eltype(strengths))
    VEL = zeros(T, size(controlpoints, 2), size(nodes, 2), 2)

    induced_velocities_from_vortex_panels_on_points!(
        VEL,
        controlpoints,
        nodes,
        nodemap,
        influence_lengths,
        strengths,
        integration_options;
        integration_caches=integration_caches,
    )

    return VEL
end

"""
    induced_velocities_from_vortex_panels_on_points!(
        VEL,
        controlpoint,
        node,
        nodemap,
        influence_length,
        strength,
        integration_options;
        integration_caches=nothing,
    )

In-place version of `induced_velocities_from_vortex_panels_on_points`.
"""
function induced_velocities_from_vortex_panels_on_points!(
    VEL,
    controlpoint,
    node,
    nodemap,
    influence_length,
    strength,
    integration_options;
    integration_caches=nothing,
)

    # vel = zeros(VEL, 2, 2)
    if isnothing(integration_caches)
        #integration_cache = zeros(eltype(node), 20)
        nominal_integration_cache = allocate_integration_containers(
            integration_options.nominal, VEL
        )
        singular_integration_cache = allocate_integration_containers(
            integration_options.singular, VEL
        )
    else
        nominal_integration_cache = integration_caches.nominal
        singular_integration_cache = integration_caches.singular
    end

    # loop through panels doing the influencing
    for (j, (nmap, lj, gammaj)) in
        enumerate(zip(eachcol(nodemap), influence_length, eachcol(strength)))
        # Loop through control points being influenced
        for (i, cpi) in enumerate(eachcol(controlpoint))
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

            for k in 1:2
                # fill the Matrix
                # TODO: is having this gammaj here even necessary, or is it always just 1? If so, really need to just remove it since we're allocating large arrays every time we call this function.
                VEL[i, Int(nmap[k]), :] += gammaj[k] * vel[k, :]
            end #for k
        end #for i
    end #for j

    return VEL
end

#---------------------------------#
#     Source Panels on Points     #
#---------------------------------#

"""
    induced_velocities_from_source_panels_on_points(
        controlpoints,
        nodes,
        nodemap,
        influence_lengths,
        strengths,
        integration_options;
        integration_caches=nothing,
    )


Calculate axial and radial components of induced velocity for a set of control points due to a set of axisymmetric source panels (bands)

Used for getting the unit induced velocities due to the body panels on the rotor/wake as well as the unit induced velocity due to the wake on the body/rotor.

# Arguments
- `controlpoints::Matrix{Float}` [z r] coordinates of points being influenced
- `nodes::Matrix{Float}` : [z r] coordinates of vortex rings
- `nodemap::Matrix{Int}` : mapping from panel index to associated node indices
- `influence_lengths::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
- `strengths::Matrix{Float}` : source constant strength values
- `integration_options::IntegrationOptions` : integration options

# Returns:
- `VEL::Array{Float}` : N-controlpoint x N-node x [vz, vr] array of induced velocity components
"""
function induced_velocities_from_source_panels_on_points(
    controlpoints,
    nodes,
    nodemap,
    influence_lengths,
    strengths,
    integration_options;
    integration_caches=nothing,
)

    # Initialize
    T = promote_type(eltype(controlpoints), eltype(nodes), eltype(strengths))
    VEL = zeros(T, size(controlpoints, 2), size(nodes, 2), 2)

    induced_velocities_from_source_panels_on_points!(
        VEL,
        controlpoints,
        nodes,
        nodemap,
        influence_lengths,
        strengths,
        integration_options;
        integration_caches=integration_caches,
    )

    return VEL
end

"""
    induced_velocities_from_source_panels_on_points!(
        VEL,
        controlpoint,
        node,
        nodemap,
        influence_length,
        strength,
        integration_options;
        integration_caches=nothing,
    )

In-place version of `induced_velocities_from_source_panels_on_points`.
"""
function induced_velocities_from_source_panels_on_points!(
    VEL,
    controlpoint,
    node,
    nodemap,
    influence_length,
    strength,
    integration_options;
    integration_caches=nothing,
)
    # vel = zeros(VEL, 2, 2)
    if isnothing(integration_caches)
        #integration_cache = zeros(eltype(node), 20)
        nominal_integration_cache = allocate_integration_containers(
            integration_options.nominal, VEL
        )
        singular_integration_cache = allocate_integration_containers(
            integration_options.singular, VEL
        )
    else
        nominal_integration_cache = integration_caches.nominal
        singular_integration_cache = integration_caches.singular
    end

    # loop through panels doing the influencing
    for (j, (nmap, lj, sigmaj)) in
        enumerate(zip(eachcol(nodemap), influence_length, eachcol(strength)))

        # Loop through control points being influenced
        for (i, cpi) in enumerate(eachcol(controlpoint))
            n1 = view(node, :, Int(nmap[1]))
            n2 = view(node, :, Int(nmap[2]))

            # check of self-induced:
            if isapprox(cpi, 0.5 * (n1 .+ n2))
                # if so:
                # vel .= self_source_panel_integration(n1, n2, lj, cpi)
                vel = StaticArrays.SMatrix{2,2}(
                    self_source_panel_integration(
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
                # vel .= nominal_source_panel_integration(n1, n2, lj, cpi)
                vel = StaticArrays.SMatrix{2,2}(
                    nominal_source_panel_integration(
                        integration_options.nominal,
                        n1,
                        n2,
                        lj,
                        cpi,
                        nominal_integration_cache,
                    ),
                )
            end

            for k in 1:2
                # fill the Matrix
                VEL[i, Int(nmap[k]), :] += sigmaj[k] * vel[k, :]
            end #for k
        end #for i
    end #for j

    return VEL
end

#---------------------------------#
#     Trailing Edge Gap Panel     #
#---------------------------------#

"""
    induced_velocities_from_trailing_edge_gap_panel!(
        VEL,
        controlpoint,
        tenode,
        teinfluence_length,
        tendotn,
        tencrossn,
        teadjnodeidxs,
        integration_options;
        wake=false,
        integration_caches=nothing,
    )


Calculate axial and radial components of induced velocity for a set of control points due to any trailing edge gap panels.

Used for getting the unit induced velocities due to the body body trailing edge gap panels on the body/rotor/wake.

Note, this function is also used to calculate the influence of the wake ends rather than modeling a semi-infinite fortex sheet.

# Arguments
- `VEL::Array{Float}` : N-controlpoint x N-node x [vz, vr] array of induced velocity components (modified in place)
- `controlpoints::Matrix{Float}` [z r] coordinates of points being influenced
- `nodes::Matrix{Float}` : [z r] coordinates of vortex rings
- `nodemap::Matrix{Int}` : mapping from panel index to associated node indices
- `influence_lengths::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
- `strengths::Matrix{Float}` : vortex constant circulation values
- `integration_options::IntegrationOptions` : integration options

# Keyword Arguments
- `wake::Bool=false` : flag to indicate if this is being used for a wake sheet.
- `integration_caches::NamedTuple=nothing` : cache used in in-place quadrature functions.
"""
function induced_velocities_from_trailing_edge_gap_panel!(
    VEL,
    controlpoint,
    tenode,
    teinfluence_length,
    tendotn,
    tencrossn,
    teadjnodeidxs,
    integration_options;
    wake=false,
    integration_caches=nothing,
)

    # vvel = zeros(eltype(AICn), 2, 2)
    # svel = zeros(eltype(AICn), 2, 2)
    if isnothing(integration_caches)
        # integration_cache = zeros(eltype(controlpoint), 20)
        nominal_integration_cache = allocate_integration_containers(
            integration_options.nominal, VEL
        )
        singular_integration_cache = allocate_integration_containers(
            integration_options.singular, VEL
        )
    else
        nominal_integration_cache = integration_caches.nominal
        singular_integration_cache = integration_caches.singular
    end

    # Loop through control points being influenced
    for (i, cpi) in enumerate(eachcol(controlpoint))
        # loop through bodies
        for (j, (lj, ndn, ncn, nmap)) in enumerate(
            zip(
                teinfluence_length,
                eachcol(tendotn),
                eachcol(tencrossn),
                eachcol(teadjnodeidxs),
            ),
        )
            n1 = view(tenode, j, 1, :)
            n2 = view(tenode, j, 2, :)

            # get unit induced velocity from the panel onto the control point
            if isapprox(cpi, 0.5 * (n1 .+ n2))
                # if so:
                vvel = StaticArrays.SMatrix{2,2}(
                    self_vortex_panel_integration(
                        integration_options.singular,
                        n1,
                        n2,
                        lj,
                        cpi,
                        singular_integration_cache,
                    ),
                )
                svel = StaticArrays.SMatrix{2,2}(
                    self_source_panel_integration(
                        integration_options.singular,
                        n1,
                        n2,
                        lj,
                        cpi,
                        nominal_integration_cache,
                    ),
                )
            else
                # if not:
                vvel = StaticArrays.SMatrix{2,2}(
                    nominal_vortex_panel_integration(
                        integration_options.nominal,
                        n1,
                        n2,
                        lj,
                        cpi,
                        singular_integration_cache,
                    ),
                )
                svel = StaticArrays.SMatrix{2,2}(
                    nominal_source_panel_integration(
                        integration_options.nominal,
                        n1,
                        n2,
                        lj,
                        cpi,
                        nominal_integration_cache,
                    ),
                )
            end

            for k in 1:2
                # fill the Matrix
                VEL[i, Int(nmap[k]), :] += ndn[k] * vvel[k, :] + ncn[k] * svel[k, :]
                if wake
                    # wake "TE Panels" only have the vortex influence
                    VEL[i, Int(nmap[k]), :] += ndn[k] * vvel[k, :]
                end
            end #for k
        end #for j
    end #for i

    return VEL
end
