#########################################
#                                       #
#           Panel Generation            #
#                                       #
#########################################
"""
avoid issues with single bodies being input not as vectors
"""
function generate_panels(coordinates::Matrix{TF}; kwargs...) where {TF}
    return generate_panels([coordinates]; kwargs...)
end

"""
generates NamedTuple of panel geometry items from a vector of matrices of coordinates
assumes annular airfoils are given first in coordinates array (for tracking kutta conditions)
"""
function generate_panels(
    coordinates::Vector{Matrix{TF}}; itpanelscale=0.05, axistol=1e-15, tetol=1e1 * eps()
) where {TF}

    ## -- SETUP -- ##
    # number of panels to generate for each body
    npanel = [length(eachrow(c)) - 1 for c in coordinates]
    # number of bodies
    nbodies = length(npanel)
    # total number of panels in system
    totpanel = sum(npanel)
    # total number of nodes in system
    totnode = totpanel + nbodies

    # - Initialize Outputs - #
    # control points
    controlpoint = zeros(TF, totpanel, 2) # panel, x-r
    # nodes
    node = zeros(TF, totnode, 2) # node, x-r
    # node map
    nodemap = zeros(Int, totpanel, 2) # node, x-r
    #TODO: is endpoint information still needed? may be able to clean up.
    # endpoints of bodies, rotor, or wakes
    endpoints = zeros(TF, length(coordinates), 2, 2) # TE, upper-lower, x-r
    # indices of endpoints
    endpointidxs = ones(Int, length(coordinates), 2) # lower idx, upper idx, lower or upper
    # panel lengths
    influence_length = zeros(TF, totpanel)
    # panel unit normals
    normal = zeros(TF, totpanel, 2)
    # panel unit tangents
    tangent = zeros(TF, totpanel, 2)

    # initialize index for entire array
    pidx = 1 # panel index
    nidx = 1 # node index (+1 extra after each body)

    # loop through bodies
    for (ib, c) in enumerate(coordinates)

        # Separate coordinates
        x = view(c, :, 1)
        r = view(c, :, 2)

        # Check if any r coordinates are negative (not allowed in axisymmetric method)
        @assert all(r -> r >= 0.0, r) "Some coordinates have negative radial components."

        ## -- Loop Through Coordinates -- ##
        for ip in 1:npanel[ib]

            # Get nodes (panel edges)
            node[nidx, :] = [x[ip]; r[ip]]
            node[nidx + 1, :] = [x[ip + 1]; r[ip + 1]]
            nodemap[pidx,:] .= [nidx; nidx+1]

            # Calculate control point (panel center)
            controlpoint[pidx, :] .= [0.5 * (x[ip] + x[ip + 1]); 0.5 * (r[ip] + r[ip + 1])]

            # Calculate panel length
            influence_length[pidx] += get_r(node[nidx, :], node[nidx + 1, :])[2]

            # Calculate panel unit normal
            normal[pidx, :] = get_panel_normal(get_r(c[ip, :], c[ip + 1, :])...)

            # Calculate panel unit tangent
            tangent[pidx, :] = get_panel_tangent(get_r(c[ip, :], c[ip + 1, :])...)

            # - Get endpoints - # TODO: still necessary?
            if ip == 1
                #lower
                endpoints[ib, 1, :] = [x[ip] r[ip]]
                endpointidxs[ib, 1] = nidx
            elseif ip == npanel[ib]
                #upper
                endpoints[ib, 2, :] = [x[ip + 1] r[ip + 1]]
                endpointidxs[ib, 2] = nidx + 1
            end

            # iterate "global" panel index
            pidx += 1
            nidx += 1
        end

        # iterate "global" node index, since there is +1 nodes than panels
        nidx += 1
    end

    ## - Internal Panel Stuff - #
    #TODO: DFDC uses this, it may or may not be necessary
    #itcontrolpoint = zeros(TF, nbodies, 2)
    #itnormal = zeros(TF, nbodies, 2)

    #for ib in 1:nbodies

    #    #TODO: for more robust, consider splining the coordinates, then getting the midpoint of the bodies, the normals can probably just be arbitrary, [1,0] should work, then there's no axial flow on the panel.

    #    #rename for convenience
    #    p1id = endpointidxs[ib, 1]
    #    pnid = endpointidxs[ib, 2]

    #    # get node coordinates
    #    n11 = node[p1id, 1, :]
    #    n12 = node[p1id, 2, :]
    #    nn1 = node[pnid, 1, :]
    #    nn2 = node[pnid, 2, :]

    #    xtan = n11[1] - n12[1] + nn2[1] - nn1[1]
    #    rtan = n11[2] - n12[2] + nn2[2] - nn1[2]
    #    stan = sqrt(xtan^2 + rtan^2)

    #    lenbar = 0.5 * (influence_length[p1id] + influence_length[pnid])

    #    itcontrolpoint[ib, 1] =
    #        0.5 * (n11[1] + nn2[1]) - itpanelscale * lenbar * xtan / stan
    #    itcpr = 0.5 * (n11[2] + nn2[2]) - itpanelscale * lenbar * rtan / stan
    #    #TODO: need to update this to be more robust for center bodies
    #    itcontrolpoint[ib, 2] = itcpr < axistol ? 1e-3 : itcpr

    #    itnormal[ib, 1] = xtan / stan
    #    itnormal[ib, 2] = rtan / stan
    #end

    return (;
        nbodies,
        npanels=totpanel,
        controlpoint,
        node,
        nodemap,
        influence_length,
        normal,
        tangent,
        endpoints,
        endpointidxs,
        # itcontrolpoint,
        # itnormal,
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
