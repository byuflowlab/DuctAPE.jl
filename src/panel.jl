#########################################
#                                       #
#           Panel Generation            #
#                                       #
#########################################
"""
avoid issues with single bodies being input not as vectors
"""
function generate_panels(coordinates::Matrix{TF}; kwargs...) where {TF}
    return generate_panels([coordinates], kwargs...)
end

"""
generates NamedTuple of panel geometry items from a vector of matrices of coordinates
assumes annular airfoils are given first in coordinates array (for tracking kutta conditions)
"""
function generate_panels(
    coordinates::Vector{Matrix{TF}};
    body=false,
    itpanelscale=0.05,
    axistol=1e-15,
    tetol=1e1 * eps(),
) where {TF}

# TODO: consider splitting out into multiple dispatches rather than having all the ifs for the body case. they're different enough that it probably makes sense to split them out.

    ## -- SETUP -- ##
    # Get total number of panels (sum of nedges-1 for each body)
    npanel = [length(eachrow(c)) - 1 for c in coordinates]
    nbodies = length(npanel)
    totpanel = sum(npanel)
    totnode = totpanel + nbodies

    # - Initialize Outputs - #
    controlpoint = zeros(TF, body ? totnode : totpanel, 2) # panel, x-r
    node = zeros(TF, totnode, 2) # node, x-r
    endpoints = zeros(TF, length(coordinates), 2, 2) # TE, upper-lower, x-r
    endpointidxs = ones(Int, length(coordinates), 2) # lower idx, upper idx, lower or upper
    influence_length = zeros(TF, totnode)
    normal = zeros(TF, body ? totnode : totpanel, 2)
    tangent = zeros(TF, body ? totnode : totpanel, 2)

    # initialize index for entire array
    pidx = 1 # panel index
    nidx = 1 # node index (+1 extra after each body)

    # loop through bodies
    for (ib, c) in enumerate(coordinates)

        # Separate coordinates
        x = view(c,:, 1)
        r = view(c,:, 2)

        # Check if any r coordinates are negative (not allowed in axisymmetric method)
        @assert all(r -> r >= 0.0, r)

        ## -- Loop Through Coordinates -- ##
        for ip in 1:npanel[ib]

            # Get nodes (panel edges)
            node[nidx, :] = [x[ip]; r[ip]]
            node[nidx + 1, :] = [x[ip + 1]; r[ip + 1]]

            # Calculate control point (panel center)
            controlpoint[pidx, :] .= [0.5 * (x[ip] + x[ip + 1]); 0.5 * (r[ip] + r[ip + 1])]

            # Calculate panel length
            influence_length[nidx] += get_r(node[nidx, :], controlpoint[pidx, :])[2]
            influence_length[nidx + 1] += get_r(node[nidx + 1, :], controlpoint[pidx, :])[2]

            # Calculate panel unit normal
            normal[pidx, :] = get_panel_normal(get_r(c[ip, :], c[ip + 1, :])...)

            # Calculate panel unit tangent
            tangent[pidx, :] = get_panel_tangent(get_r(c[ip, :], c[ip + 1, :])...)

            # - Get endpoints and make any adjustments - #
            if ip == 1
                #lower
                endpoints[ib, 1, :] = [x[ip] r[ip]]
                # endpointidxs[ib, 1, :] = [pidx; -1]
                if body
                    endpointidxs[ib, 1] = pidx - ib + 1
                else
                    endpointidxs[ib, 1] = pidx
                end
            elseif ip == npanel[ib]
                #upper
                endpoints[ib, 2, :] = [x[ip + 1] r[ip + 1]]
                if body
                    endpointidxs[ib, 2] = pidx - ib + 1
                else
                    endpointidxs[ib, 2] = pidx
                end
                # endpointidxs[ib, 2, :] = [pidx; 1]

                # Add extra control point at TE for bodies, and shift any nodes in halfway to control point.
                if body
                    if all(
                        sign.(tangent[(pidx - npanel[ib] + 1):pidx, 1]) .==
                        sign.(abs.(tangent[(pidx - npanel[ib] + 1):pidx, 1])),
                    )
                        # if we don't wrap around, just use last point
                        controlpoint[pidx + 1, :] .= [x[ip + 1]; r[ip + 1]]
                    else
                        # if we wrap around, use average of first and last point
                        controlpoint[pidx + 1, :] .= [
                            0.5 * (x[ip + 1] + x[1])
                            0.5 * (r[ip + 1] + r[1])
                        ]

                        # shift first node location.
                        # Assumes duct coordinates are given first
                        node[1, :] = 0.5 * (node[1, :] .+ controlpoint[1, :])
                    end

                    # - shift last node location.
                    node[nidx + 1, :] = 0.5 * (node[nidx + 1, :] .+ controlpoint[pidx, :])

                    # manually set tangent and normal for extra control point
                    tangent[pidx + 1, :] .= [1.0; 0.0]
                    normal[pidx + 1, :] .= [0.0; 1.0]
                end
            end

            # iterate "global" panel index
            pidx += 1
            nidx += 1
        end
        if body
            pidx += 1
        end
        nidx += 1
    end

    # # - Trailing Edge Nodes - #
    # TEnodes = []
    # for t in 1:length(endpoints[:, 1, 1])
    #     if isapprox(endpoints[t, 1, :], endpoints[t, 2, :]; atol=tetol)
    #         push!(TEnodes, (; pos=endpoints[t, 1, :], idx=endpointidxs[t, 1], sign=1)) # lower is negative
    #         push!(TEnodes, (; pos=endpoints[t, 2, :], idx=endpointidxs[t, 2], sign=-1)) #upper is positive
    #     end
    # end

    ## - Internal Panel Stuff - #
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

    # # - ∇μ prep stuff - #
    # # get min and max points to split up geometry
    # xmm = []
    # rmm = []
    # for ib in 1:nbodies
    #     brange = endpointidxs[ib, 1]:endpointidxs[ib, 2]
    #     append!(xmm, findlocalmax(controlpoint[brange, 1]))
    #     append!(xmm, findlocalmin(controlpoint[brange, 1]))
    #     append!(rmm, findlocalmax(controlpoint[brange, 2]))
    #     append!(rmm, findlocalmin(controlpoint[brange, 2]))
    # end
    # sort!(xmm)
    # sort!(rmm)

    # # put together ranges over which to take gradients
    # dxrange = [[xmm[i] xmm[i + 1] 0] for i in 1:(length(xmm) - 1)]
    # drrange = [[rmm[i] rmm[i + 1] 0] for i in 1:(length(rmm) - 1)]

    # # identify which ranges will need to be reversed
    # for xr in dxrange
    #     if controlpoint[xr[1], 1] > controlpoint[xr[2], 1]
    #         xr[3] = 1
    #     end
    # end
    # for rr in drrange
    #     if controlpoint[rr[1], 2] > controlpoint[rr[2], 2]
    #         rr[3] = 1
    #     end
    # end

    return (;
        controlpoint,
        node,
        influence_length,
        normal,
        tangent,
        endpoints,
        endpointidxs,
        # TEnodes,
        npanels=totpanel,
        # itcontrolpoint,
        # itnormal,
        nbodies,
        # dxrange,
        # drrange,
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
