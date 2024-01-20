#=
Influence matrices for panels on points
given in terms of raw velocity components
=#

######################################################################
#                                                                    #
#          Matrix of velocities induced by panels on points          #
#                                                                    #
######################################################################

#---------------------------------#
#     Vortex Panels on Points     #
#---------------------------------#

"""
out of place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric vortex panels (bands)

Used for getting the unit induced velocities due to the body panels on the rotor/wake as well as the unit induced velocity due to the wake on the body/rotor.

# Arguments:
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `node::Matrix{Float}` : [z r] coordinates of vortex rings
- `influence_length::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
- `gamma::Vector{Float}` : vortex constant circulation values

# Returns:
- `AIC::Array{Float}` : N-controlpoint x N-node x [vz, vr] array of induced velocity components
"""
function induced_velocities_from_vortex_panels_on_points(
    controlpoints, nodes, nodemap, influence_lengths, strengths
)

    # Initialize
    T = promote_type(eltype(controlpoints), eltype(nodes), eltype(strengths))
    VEL = zeros(T, size(controlpoints, 2), size(nodes, 2), 2)

    induced_velocities_from_vortex_panels_on_points!(
        VEL, controlpoints, nodes, nodemap, influence_lengths, strengths
    )

    return VEL
end

"""
in place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric vortex panels (bands)

Used for getting the unit induced velocities due to the body panels on the rotor/wake as well as the unit induced velocity due to the wake on the body/rotor.

# Arguments:
- `VEL::Array{Float}` : N-controlpoint x N-node x [vz, vr] array of induced velocity components
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `node::Matrix{Float}` : [z r] coordinates of vortex rings
- `influence_length::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
- `gamma::Vector{Float}` : vortex constant circulation values
"""
function induced_velocities_from_vortex_panels_on_points!(
    VEL, controlpoint, node, nodemap, influence_length, strength
)
    # vel = zeros(eltype(VEL), 2, 2)

    # loop through panels doing the influencing
    for (j, (nmap, lj, gammaj)) in
        enumerate(zip(eachcol(nodemap), influence_length, eachcol(strength)))
        # Loop through control points being influenced
        for (i, cpi) in enumerate(eachcol(controlpoint))
            n1 = view(node, :, nmap[1])
            n2 = view(node, :, nmap[2])

            # check of self-induced:
            if isapprox(cpi, 0.5 * (n1 .+ n2))
                # if so:
                # vel .= self_vortex_panel_integration(n1, n2, lj, cpi)
                vel = StaticArrays.SMatrix{2,2}(
                    self_vortex_panel_integration(n1, n2, lj, cpi)
                )
            else
                # if not:
                # vel .= nominal_vortex_panel_integration(n1, n2, lj, cpi)
                vel = StaticArrays.SMatrix{2,2}(
                    nominal_vortex_panel_integration(n1, n2, lj, cpi)
                )
            end

            for k in 1:2
                # fill the Matrix
                VEL[i, nmap[k], :] += gammaj[k] * vel[k, :]
            end #for k
        end #for i
    end #for j

    return VEL
end

#---------------------------------#
#     Source Panels on Points     #
#---------------------------------#

"""
out of place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric source panels (bands)

Used for getting the unit induced velocities due to the body panels on the rotor/wake as well as the unit induced velocity due to the wake on the body/rotor.

# Arguments:
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `node::Matrix{Float}` : [z r] coordinates of source rings
- `influence_length::Vector{Float}` : lengths over which source ring influence is applied on the surface.
- `gamma::Vector{Float}` : source constant circulation values

# Returns:
- `AIC::Array{Float}` : N-controlpoint x N-node x [vz, vr] array of induced velocity components
"""
function induced_velocities_from_source_panels_on_points(
    controlpoints, nodes, nodemap, influence_lengths, strengths
)

    # Initialize
    T = promote_type(eltype(controlpoints), eltype(nodes), eltype(strengths))
    VEL = zeros(T, size(controlpoints, 2), size(nodes, 2), 2)

    induced_velocities_from_source_panels_on_points!(
        VEL, controlpoints, nodes, nodemap, influence_lengths, strengths
    )

    return VEL
end

"""
in place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric source panels (bands)

Used for getting the unit induced velocities due to the body panels on the rotor/wake as well as the unit induced velocity due to the wake on the body/rotor.

# Arguments:
- `VEL::Array{Float}` : N-controlpoint x N-node x [vz, vr] array of induced velocity components
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `node::Matrix{Float}` : [z r] coordinates of source rings
- `influence_length::Vector{Float}` : lengths over which source ring influence is applied on the surface.
- `gamma::Vector{Float}` : source constant circulation values
"""
function induced_velocities_from_source_panels_on_points!(
    VEL, controlpoint, node, nodemap, influence_length, strength
)
    # vel = zeros(eltype(VEL), 2, 2)

    #TODO; for speedups, update panel initialization to flip rows and columns such that these functions use eachcol rather than eachrow

    # loop through panels doing the influencing
    for (j, (nmap, lj, sigmaj)) in
        enumerate(zip(eachcol(nodemap), influence_length, eachcol(strength)))
        # Loop through control points being influenced
        for (i, cpi) in enumerate(eachcol(controlpoint))
            n1 = view(node, :, nmap[1])
            n2 = view(node, :, nmap[2])

            # check of self-induced:
            if isapprox(cpi, 0.5 * (n1 .+ n2))
                # if so:
                # vel .= self_source_panel_integration(n1, n2, lj, cpi)
                vel = StaticArrays.SMatrix{2,2}(
                    self_source_panel_integration(n1, n2, lj, cpi)
                )
            else
                # if not:
                # vel .= nominal_source_panel_integration(n1, n2, lj, cpi)
                vel = StaticArrays.SMatrix{2,2}(
                    nominal_source_panel_integration(n1, n2, lj, cpi)
                )
            end

            for k in 1:2
                # fill the Matrix
                VEL[i, nmap[k], :] += sigmaj[k] * vel[k, :]
            end #for k
        end #for i
    end #for j

    return VEL
end

#---------------------------------#
#     Trailing Edge Gap Panel     #
#---------------------------------#
"""
"""
function induced_velocities_from_trailing_edge_gap_panel!(
    VEL, controlpoint, tenode, teinfluence_length, tendotn, tencrossn, teadjnodeidxs
)

    # vvel = zeros(eltype(AICn), 2, 2)
    # svel = zeros(eltype(AICn), 2, 2)

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
                    self_vortex_panel_integration(n1, n2, lj, cpi)
                )
                svel = StaticArrays.SMatrix{2,2}(
                    self_source_panel_integration(n1, n2, lj, cpi)
                )
            else
                # if not:
                vvel = StaticArrays.SMatrix{2,2}(
                    nominal_vortex_panel_integration(n1, n2, lj, cpi)
                )
                svel = StaticArrays.SMatrix{2,2}(
                    nominal_source_panel_integration(n1, n2, lj, cpi)
                )
            end

            for k in 1:2
                # fill the Matrix
                VEL[i, nmap[k], :] += ndn[k] * vvel[k, :] + ncn[k] * svel[k, :]
            end #for k
        end #for j
    end #for i

    return VEL
end

######################################################################
#                                                                    #
#    Vector of total (sum of) velocity induced by panels on points   #
#                                                                    #
######################################################################

#TODO: The rest of these need to be updated still

##---------------------------------#
##     Vortex Panels on Points     #
##---------------------------------#

#"""
#out of place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric vortex rings but only returning sum of velocities on each control point

#Used in calculating velocities on body surfaces in preparation to obtain tangential components and eventually pressure distributions.
#Note: there is probably a more efficient way to achieve this functionality.

## Arguments:
#- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
#- `node::Matrix{Float}` : [z r] coordinates of vortex rings
#- `influence_length::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
#- `gamma::Vector{Float}` : vortex constant circulation values

## Returns:
#- `V::Array{Float}` : N-controlpoint x [vz, vr] array of summed induced velocity components due to vortex rings
#"""
#function total_velocities_induced_by_vortex_panels(
#    controlpoints::AbstractMatrix{T1},
#    nodes::AbstractMatrix{T2},
#    influence_lengths::AbstractVector{T3},
#    strengthal_velocities_induced_by_vortex_panelss::AbstractArray{T4},
#) where {T1,T2,T3,T4}

#    # Initialize
#    T = promote_type(T1, T2, T3, T4)
#    V = zeros(T, size(controlpoints, 1), 2)

#    total_velocities_induced_by_vortex_panels!(
#        V, controlpoints, nodes, influence_lengths, strengths
#    )

#    return V
#end

#"""
#In place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric vortex rings but only returning sum of velocities on each control point

#Used in calculating velocities on body surfaces in preparation to obtain tangential components and eventually pressure distributions.
#Note: there is probably a more efficient way to achieve this functionality.

## Arguments:
#- `V::Array{Float}` : N-controlpoint x [vz, vr] array of summed induced velocity components due to vortex rings
#- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
#- `node::Matrix{Float}` : [z r] coordinates of vortex rings
#- `influence_length::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
#- `gamma::Vector{Float}` : vortex constant circulation values
#"""
#function total_velocities_induced_by_vortex_panels!(
#    V, controlpoints, nodes, influence_lengths, strengths
#)

#    # Loop through control points
#    for (i, (cpi, vel)) in enumerate(zip(eachcol(controlpoints), eachcol(V)))
#        # loop through panels doing the influencing
#        for (j, (gamma, nj, lj)) in
#            enumerate(zip(strengths, eachcol(nodes), influence_lengths))

#            # get unit induced velocity from the panel onto the control point
#            vortex_induced_velocity!(vel, cpi, nj, lj, gamma)
#        end
#    end

#    return nothing
#end

##---------------------------------#
##     Source Panels on Points     #
##---------------------------------#

#"""
#out of place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric source rings but only returning sum of velocities on each control point

#Used in calculating velocities on body surfaces in preparation to obtain tangential components and eventually pressure distributions.
#Note: there is probably a more efficient way to achieve this functionality.

## Arguments:
#- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
#- `node::Matrix{Float}` : [z r] coordinates of source rings
#- `influence_length::Vector{Float}` : lengths over which source ring influence is applied on the surface.
#- `sigma::Vector{Float}` : source constant circulation values

## Returns:
#- `V::Array{Float}` : N-controlpoint x [vz, vr] array of summed induced velocity components due to source rings
#"""
#function total_velocities_induced_by_source_panels(
#    controlpoints::AbstractMatrix{T1},
#    nodes::AbstractArray{T2},
#    influence_lengths::AbstractVector{T3},
#    strengths::AbstractArray{T4},
#) where {T1,T2,T3,T4}

#    # Initialize
#    T = promote_type(T1, T2, T3, T4)
#    V = zeros(T, size(controlpoints, 1), 2)

#    total_velocities_induced_by_source_panels!(
#        V, controlpoints, nodes, influence_lengths, strengths
#    )

#    return V
#end

#"""
#in place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric source rings but only returning sum of velocities on each control point

#Used in calculating velocities on body surfaces in preparation to obtain tangential components and eventually pressure distributions.
#Note: there is probably a more efficient way to achieve this functionality.

## Arguments:
#- `V::Array{Float}` : N-controlpoint x [vz, vr] array of summed induced velocity components due to source rings
#- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
#- `node::Matrix{Float}` : [z r] coordinates of source rings
#- `influence_length::Vector{Float}` : lengths over which source ring influence is applied on the surface.
#- `sigma::Vector{Float}` : source constant circulation values
#"""
#function total_velocities_induced_by_source_panels!(
#    V, controlpoints, nodes, influence_length, strengths
#)
#    for (i, (cpi, vel)) in enumerate(zip(eachcol(controlpoints), eachcol(V)))
#        # loop through panels doing the influencing
#        for (j, (sigma, nj, lj)) in
#            enumerate(zip(strengths, eachcol(nodes), influence_length))

#            # get unit induced velocity from the panel onto the control point
#            source_induced_velocity!(vel, cpi, nj, lj, sigma)
#        end
#    end

#    return nothing
#end

