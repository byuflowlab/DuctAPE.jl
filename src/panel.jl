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
    coordinates::Vector{Matrix{TF}}; itcpshift=0.05, axistol=1e-15, tegaptol=1e1 * eps()
) where {TF}

    ## -- SETUP -- ##
    # number of nodes according to coordinates
    nnode = [length(eachrow(c)) for c in coordinates]
    # number of panels to generate for each body
    npanel = nnode .- 1
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
    # first and last nodes of bodies, rotor, or wakes
    endnodes = zeros(TF, nbodies, 2, 2) # TE, upper-lower, x-r
    # indices of endnodess
    endnodeidxs = ones(Int, nbodies, 2) # lower idx, upper idx, lower or upper
    # indices of endpanels
    endpanelidxs = ones(Int, nbodies, 2) # lower idx, upper idx, lower or upper
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
            nodemap[pidx, :] .= [nidx; nidx + 1]

            # Calculate control point (panel center)
            controlpoint[pidx, :] .= [0.5 * (x[ip] + x[ip + 1]); 0.5 * (r[ip] + r[ip + 1])]

            # Calculate panel length
            influence_length[pidx] += get_r(node[nidx, :], node[nidx + 1, :])[2]

            # Calculate panel unit normal
            normal[pidx, :] = get_panel_normal(get_r(c[ip, :], c[ip + 1, :])...)

            # Calculate panel unit tangent
            tangent[pidx, :] = get_panel_tangent(get_r(c[ip, :], c[ip + 1, :])...)

            # - Get endpoints - #
            if ip == 1
                #lower
                endnodes[ib, 1, :] = [x[ip] r[ip]]
                endnodeidxs[ib, 1] = nidx
                endpanelidxs[ib, 1] = pidx
            elseif ip == npanel[ib]
                #upper
                endnodes[ib, 2, :] = [x[ip + 1] r[ip + 1]]
                endnodeidxs[ib, 2] = nidx + 1
                endpanelidxs[ib, 2] = pidx
            end

            # iterate "global" panel index
            pidx += 1
            nidx += 1
        end

        # iterate "global" node index, since there is +1 nodes than panels
        nidx += 1
    end

    # - Internal Panel Stuff - #
    itcontrolpoint = zeros(TF, nbodies, 2)
    itnormal = zeros(TF, nbodies, 2)
    #TODO: unnecessary?
    ittangent = zeros(TF, nbodies, 2)

    if size(node, 1) > 2
        #TODO: maybe move this into it's own in-place function
        for ib in 1:nbodies

            #TODO: for more robust, consider splining the coordinates, then getting the midpoint of the bodies, the normals can probably just be arbitrary, [1,0] should work, then there's no axial flow on the panel.

            #rename for convenience
            p1id = endnodeidxs[ib, 1]
            pnid = endnodeidxs[ib, 2]

            # get node coordinates
            n1 = node[p1id, :] #first node on panel 1 (of body ib)
            n2 = node[p1id + 1, :] #second node on panel 1 (of body ib)
            nn = node[pnid - 1, :] #first node on panel N (of body ib)
            nnp1 = node[pnid, :] #second node on pane N (of body ib)

            #TODO: is there a sign error or something here?
            #rmagN - rmag1
            xtan = nnp1[1] - nn[1] - (n2[1] - n1[1])
            # xtan = n2[1] - n1[1] + nn[1] - nnp1[1]
            rtan = nnp1[2] - nn[2] - (n2[2] - n1[2])
            # rtan = n2[2] - n1[2] + nn[2] - nnp1[2]
            stan = sqrt(xtan^2 + rtan^2)

            lenbar = 0.5 * (influence_length[p1id] + influence_length[pnid - ib])

            itcontrolpoint[ib, 1] =
                0.5 * (n1[1] + nnp1[1]) - itcpshift * lenbar * xtan / stan
            itcpr = 0.5 * (n1[2] + nnp1[2]) - itcpshift * lenbar * rtan / stan
            #TODO: need to update this to be more robust for center bodies
            itcontrolpoint[ib, 2] = itcpr < axistol ? 1e-3 : itcpr

            itnormal[ib, 1] = xtan / stan
            itnormal[ib, 2] = rtan / stan
            ittangent[ib, 1] = -rtan / stan
            ittangent[ib, 2] = xtan / stan
        end
    end

    # - Trailing Edge Panel Stuff - #
    # TODO: also consider moving this into it's own function
    tenode = zeros(TF, nbodies, 2, 2) # body, node1-2, x-r
    tenormal = zeros(TF, nbodies, 2) #body, x-r
    teinfluence_length = zeros(TF, nbodies)
    teadjnodeidxs = similar(endnodeidxs) .= 1 #can't be the same as endpoints, because we may have repeated values for non-duct bodies
    tendotn = zeros(TF, nbodies, 2) #bodies, node1,2
    tencrossn = zeros(TF, nbodies, 2) #bodies, node1,2

    for ib in 1:nbodies
        # check if signs of x-tangents are the same
        if sign(tangent[endpanelidxs[ib, 1], 1]) != sign(tangent[endpanelidxs[ib, 2], 1])
            # if not: it's a duct
            # set first node to last endnode,
            tenode[ib, 1, :] .= endnodes[ib, 2, :]
            # and second node to the first end node (keep the direction consistent)
            tenode[ib, 2, :] .= endnodes[ib, 1, :]
            # set normal as above
            tenormal[ib, :] = get_panel_normal(get_r(tenode[ib, 1, :], tenode[ib, 2, :])...)
            # set node id's to the endnode ids
            teadjnodeidxs[ib, 1] = endnodeidxs[ib, 2]
            teadjnodeidxs[ib, 2] = endnodeidxs[ib, 1]
            # set dots and crosses for each adjacent panel
            tendotn[ib, 1] = dot(tenormal[ib, :], normal[endpanelidxs[ib, 2], :])
            tendotn[ib, 2] = dot(tenormal[ib, :], normal[endpanelidxs[ib, 1], :])
            tencrossn[ib, 1] = cross2mag(tenormal[ib, :], normal[endpanelidxs[ib, 2], :])
            tencrossn[ib, 2] = cross2mag(tenormal[ib, :], normal[endpanelidxs[ib, 1], :])

        else
            # if so: it's anything else (assuming the hub nose doesn't curve inward...
            # set first node to last endnode,
            tenode[ib, 1, :] .= endnodes[ib, 2, :]
            # and second node to the axis with same axial position
            tenode[ib, 2, :] .= [endnodes[ib, 2, 1], 0.0]
            # set normal to parallel to axis
            tenormal[ib, :] = [1.0, 0.0]
            # set both adjacent node id's to the last endnode id
            teadjnodeidxs[ib, :] .= endnodeidxs[ib, 2]
            # set both dots and crosses to the same thing based on the single adjacent node.
            tendotn[ib, 1] = dot(tenormal[ib, :], normal[endpanelidxs[ib, 2], :])
            tencrossn[ib, 1] = 0.0 # unnecessary, but just in case initialization changes
            tencrossn[ib, :] .= cross2mag(tenormal[ib, :], normal[endpanelidxs[ib, 2], :])
            # tencrossn[ib, 2] = -1.0
        end
        teinfluence_length[ib] = get_r(tenode[ib, 1, :], tenode[ib, 2, :])[2]
    end

    # - Prescribed Nodes - #
    # Save the node index for nodes that are on the axis and need to be prescribed.
    prescribednodeidxs = findall(x -> abs(x) <= eps(), node[:, 2])

    return (;
        nbodies,
        npanel,
        nnode,
        totpanel,
        totnode,
        controlpoint,
        node,
        nodemap,
        influence_length,
        normal,
        tangent,
        endnodes,
        endnodeidxs,
        endpanelidxs,
        itcontrolpoint,
        itnormal,
        ittangent,
        prescribednodeidxs,
        tenode,
        tenormal,
        teinfluence_length,
        teadjnodeidxs,
        tendotn,
        tencrossn,
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
