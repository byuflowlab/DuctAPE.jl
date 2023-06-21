#########################################
#                                       #
#           Panel Generation            #
#                                       #
#########################################
"""
"""
function generate_panels( coordinates::Matrix{TF}) where {TF}

    ## -- SETUP -- ##

    # Separate coordinates
    x = coordinates[:, 1]
    r = coordinates[:, 2]

    # Check if any r coordinates are negative (not allowed in axisymmetric method)
    @assert all(r -> r >= 0.0, r)

    # - Rename for Convenience - #
    npanel = length(x) - 1

    # - Initialize Outputs - #
    TF = eltype(coordinates)
    panel_center = zeros(TF, npanel, 2)
    panel_edges = zeros(TF, npanel, 2, 2) # panel, edge, x-r
    panel_length = zeros(TF, npanel)
    panel_normal = zeros(TF, npanel, 2)
    panel_tangent = zeros(TF, npanel, 2)

    ## -- Loop Through Coordinates -- ##
    for ip in 1:npanel

        # Get nodes (panel edges)
        panel_edges[ip,:,:] = [x[ip] r[ip]; x[ip+1] r[ip+1]]

        # Calculate control point (panel center)
        panel_center[ip, :] = [0.5 * (x[ip] + x[ip + 1]); 0.5 * (r[ip] + r[ip + 1])]

        # Calculate panel length
        panel_vector, panel_length[ip] = get_r(panel_edges[ip,1,:], panel_edges[ip,2,:])

        # Calculate panel unit normal
        panel_normal[ip, :] = get_panel_normal(panel_vector, panel_length[ip])

        # Calculate panel unit tangent
        panel_tangent[ip, :] = get_panel_tangent(panel_vector, panel_length[ip])

    end

    return (; panel_center, panel_edges, panel_length, panel_normal, panel_tangent)
end

#########################################
#                                       #
#          Geometry Functions           #
#                                       #
#########################################
"""
    function get_r(node,point)

Calculate the vector, \$\\mathbf{r}\$, and distance, \$|r|\$, from the node to the evaluation point

**Arguments:**
 - `node::Array{Float}` : [x y] position of node
 - `point::Array{Float}` : [x y] position of point.

**Returns**
 - `r::Vector{Float}` : vector from node to evaluation point
 - `rmag::Float` : length of panel between node and evaluation point
"""
function get_r(node, point)

    # Need to make adjustments for sqrt(0) cases
    if isapprox(point, node)
        TF = eltype(node)
        r = zeros(TF, 2)
        rmag = TF(0.0)

        return r, rmag

    else
        # Calculate vector
        r = point .- node

        # Calculate magnitude
        rmag = sqrt(r[1]^2 + r[2]^2)

        return r, rmag
    end
end

"""
    get_panel_normal(d, dmag)

Get unit normal to panel.

**Arguments:**
 - `d::Vector{Float}` : vector from node1 to node2.
 - `dmag::Float` : panel length

"""
function get_panel_normal(d, dmag)

    # get unit tangent
    that = get_panel_tangent(d, dmag)

    # use fancy trick to rotate to be unit normal
    nhat = [-that[2]; that[1]]

    return nhat
end

"""
    get_panel_tangent(d, dmag)

Get unit tangent to panel.

**Arguments:**
 - `d::Vector{Float}` : vector from node1 to node2.
 - `dmag::Float` : panel length

"""
function get_panel_tangent(d, dmag)
    return (dmag == 0.0) ? [0.0; 0.0] : (d / dmag)
end

#########################################
#                                       #
#          "Mesh" Generation            #
#                                       #
#########################################
"""

calculate "mesh" geometry without creating a mesh object
"""
function calculate_xrm(influencing_point, affected_point)
    xi = (affected_point[1] - influencing_point[1]) / influencing_point[2]
    rho = affected_point[2] / influencing_point[2]
    m = (4.0 * rho) / (xi^2 + (rho + 1)^2)
    rj = influencing_point[2]

    return xi, rho, m, rj
end

"""
    get_relative_geometry(influence_panels, affect_panels; kwargs)

Function similar to FLOWFoil's meshing function, but only creates the mesh one way rather
than both ways (i.e. doesn't loop through both inputs as both sources and targets).

 # Arguments
 - `influence_panels::Vector{FLOWFoil.AxisymmetricPanel}` : vector of panel objects doing the influencing
 - `affect_panels::Vector{FLOWFoil.AxisymmetricPanel}` : vector of panel objects being affected

Multiple Dispatch allows for single panel objects as one or both inputs as well if there is only one body influencing and/or being affected.

 # Returns:
 - `mesh::OneWayMesh` : OneWayMesh object with relative geometry from influence to affected panels.
"""
function get_relative_geometry(influence_panels, affect_panels)

    ### --- Convenience Variables --- ###
    nbodies_i = length(influence_panels)
    nbodies_a = length(affect_panels)
    npanels_i = [influence_panels[i].npanels for i in 1:nbodies_i]
    npanels_a = [affect_panels[i].npanels for i in 1:nbodies_a]
    total_influence_panels = sum(npanels_i)
    total_affect_panels = sum(npanels_a)

    # - Define Body Indexing - #

    #find starting indices for each body
    cspanels_i = cumsum(npanels_i)
    cspanels_a = cumsum(npanels_a)

    # put together index ranges of panels for each body
    panel_indices_i = [
        (1 + (i == 1 ? 0 : cspanels_i[i - 1])):(cspanels_i[i]) for i in 1:nbodies_i
    ]
    panel_indices_a = [
        (1 + (i == 1 ? 0 : cspanels_a[i - 1])):(cspanels_a[i]) for i in 1:nbodies_a
    ]

    # - Map indices - #
    mesh2panel_i = reduce(vcat, [1:npanels_i[i] for i in 1:nbodies_i])
    mesh2panel_a = reduce(vcat, [1:npanels_a[i] for i in 1:nbodies_a])

    ### --- Initialize Vectors --- ###
    TF = typeof(sum([influence_panels[i].panel_length[1] for i in 1:nbodies_i]))

    ### --- General Mesh Fields --- ###

    # x-component of normalized distance from influencing panel center to field point
    xi = zeros(TF, (total_affect_panels, total_influence_panels))

    # r-component of normalized distance from influencing panel center to field point
    rho = zeros(TF, (total_affect_panels, total_influence_panels))

    # variable used in elliptic function calculations
    k2 = zeros(TF, (total_affect_panels, total_influence_panels))

    ### --- Loop through bodies --- ###
    for m in 1:nbodies_a
        for n in 1:nbodies_i
            ### --- Loop through panels --- ###
            for i in panel_indices_a[m]
                for j in panel_indices_i[n]

                    #TODO: this will likely change depending on how the panel method gets set up.
                    # Get xi-locations of influencing and influenced panels
                    xi = affect_panels[m].panel_center[mesh2panel_a[i], 1]
                    xj = influence_panels[n].panel_center[mesh2panel_i[j], 1]

                    # Get rho-locations of influencing and influenced panels
                    ri = affect_panels[m].panel_center[mesh2panel_a[i], 2]
                    rj = influence_panels[n].panel_center[mesh2panel_i[j], 2]

                    # Get normalized relative geometry
                    xi[i,j], rho[i,j], k2[i,j], _ = calculate_xrm(influencing_point, affected_point)

                end #for jth influencing panel
            end #for ith influenced panel
        end #for nth influencing body
    end #for mth influenced body

    # Return Mesh
    # return AxisymmetricMesh(nbodies, panel_indices, mesh2panel,xi, rho, k2)
    return (;
        nbodies_i,
        nbodies_a,
        panel_indices_i,
        panel_indices_a,
        mesh2panel_i,
        mesh2panel_a,
        xi,
        rho,
        k2,
    )
end

function generate_one_way_mesh(influence_panels::TP, affect_panels) where {TP<:NamedTuple}
    return generate_one_way_mesh([influence_panels], affect_panels)
end

function generate_one_way_mesh(influence_panels, affect_panels::TP) where {TP<:NamedTuple}
    return generate_one_way_mesh(influence_panels, [affect_panels])
end

function generate_one_way_mesh(influence_panels::TP, affect_panels::TP) where {TP<:NamedTuple}
    return generate_one_way_mesh([influence_panels], [affect_panels])
end
