#########################################
#                                       #
#           Panel Generation            #
#                                       #
#########################################
"""
avoid issues with single bodies being input not as vectors
"""
function generate_panels(coordinates::Matrix{TF}) where {TF}
    return generate_panels([coordinates])
end

"""
generates NamedTuple of panel geometry items from a vector of matrices of coordinates
assumes annular airfoils are given first in coordinates array (for tracking kutta conditions)
"""
function generate_panels(coordinates::Vector{Matrix{TF}}) where {TF}

    ## -- SETUP -- ##
    # Get total number of panels (sum of nedges-1 for each body)
    npanel = [length(eachrow(c)) - 1 for c in coordinates]
    totpanel = sum(npanel)

    # - Initialize Outputs - #
    controlpoint = zeros(TF, totpanel, 2)
    nodes = zeros(TF, totpanel, 2, 2) # panel, edge, x-r
    endpoints = zeros(TF, length(coordinates), 2, 2) # TE, upper-lower, x-r
    endpointidxs = ones(Int, length(coordinates), 2) # lower idx, upper idx
    panel_length = zeros(TF, totpanel)
    normal = zeros(TF, totpanel, 2)
    tangent = zeros(TF, totpanel, 2)

    # initialize index for entire array
    pidx = 1

    # loop through bodies
    for (ib, c) in enumerate(coordinates)

        # Separate coordinates
        x = c[:, 1]
        r = c[:, 2]

        # Check if any r coordinates are negative (not allowed in axisymmetric method)
        @assert all(r -> r >= 0.0, r)

        ## -- Loop Through Coordinates -- ##
        for ip in 1:npanel[ib]
            if ip == 1
                endpoints[ib, 1, :] = [x[ip] r[ip]]
                endpointidxs[ib, 1] = pidx
            elseif ip == npanel[ib]
                endpoints[ib, 2, :] = [x[ip + 1] r[ip + 1]]
                endpointidxs[ib, 2] = pidx
            end

            # Get nodes (panel edges)
            nodes[pidx, :, :] = [x[ip] r[ip]; x[ip + 1] r[ip + 1]]

            # Calculate control point (panel center)
            controlpoint[pidx, :] = [0.5 * (x[ip] + x[ip + 1]); 0.5 * (r[ip] + r[ip + 1])]

            # Calculate panel length
            panel_vector, panel_length[pidx] = get_r(nodes[pidx, 1, :], nodes[pidx, 2, :])

            # Calculate panel unit normal
            normal[pidx, :] = get_panel_normal(panel_vector, panel_length[pidx])

            # Calculate panel unit tangent
            tangent[pidx, :] = get_panel_tangent(panel_vector, panel_length[pidx])

            # iterate "global" panel index
            pidx += 1
        end
    end

    return (;
        controlpoint,
        nodes,
        len=panel_length,
        normal,
        tangent,
        endpoints,
        endpointidxs,
        npanels=totpanel,
    )
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
