#=

Functions for generating meshes (relational geometry between objects)

Authors: Judd Mehr,

=#

"""
"""
struct OneWayMesh
    n_influence_bodies
    n_affect_bodies
    influence_panel_indices
    affect_panel_indices
    mesh2panel_influence
    mesh2panel_affect
    x
    r
    m
end

function generate_one_way_mesh(
    influence_panels::TP, affect_panels; singularity="vortex"
) where {TP<:ff.Panel}
    return generate_one_way_mesh([influence_panels], affect_panels; singularity=singularity)
end

function generate_one_way_mesh(
    influence_panels, affect_panels::TP; singularity="vortex"
) where {TP<:ff.Panel}
    return generate_one_way_mesh(influence_panels, [affect_panels]; singularity=singularity)
end

function generate_one_way_mesh(
    influence_panels::TP, affect_panels::TP; singularity="vortex"
) where {TP<:ff.Panel}
    return generate_one_way_mesh(
        [influence_panels], [affect_panels]; singularity=singularity
    )
end
"""
function similar to flowfoil's meshing function, but only creates the mesh one way rather than both ways (i.e. doesn't loop through both inputs
"""
function generate_one_way_mesh(influence_panels, affect_panels; singularity="vortex")

    #TODO: copy over meshing function from flowfoil and adjust

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
    x = zeros(TF, (total_affect_panels, total_influence_panels))

    # r-component of normalized distance from influencing panel center to field point
    r = zeros(TF, (total_affect_panels, total_influence_panels))

    # variable used in elliptic function calculations
    k2 = zeros(TF, (total_affect_panels, total_influence_panels))

    ### --- Loop through bodies --- ###
    for m in 1:nbodies_a
        for n in 1:nbodies_i
            ### --- Loop through panels --- ###
            for i in panel_indices_a[m]
                for j in panel_indices_i[n]

                    # Get x-locations of influencing and influenced panels
                    xi = affect_panels[m].panel_center[mesh2panel_a[i], 1]
                    xj = influence_panels[n].panel_center[mesh2panel_i[j], 1]

                    # Get r-locations of influencing and influenced panels
                    ri = affect_panels[m].panel_center[mesh2panel_a[i], 2]
                    rj = influence_panels[n].panel_center[mesh2panel_i[j], 2]

                    # Calculate normalized distance components for current set of panels
                    if singularity == "vortex"
                        x[i, j] = (xi - xj) / rj
                    elseif singularity == "source"
                        x[i, j] = (xi - xj) / ri
                    else
                        @error "no singularity of type $(singularity)"
                    end
                    r[i, j] = ri / rj

                    # Calculate the k^2 value for the elliptic integrals
                    k2[i, j] = 4.0 * r[i, j] / (x[i, j]^2 + (r[i, j] + 1.0)^2)
                end #for jth influencing panel
            end #for ith influenced panel
        end #for nth influencing body
    end #for mth influenced body

    # Return Mesh
    # return AxisymmetricMesh(nbodies, panel_indices, mesh2panel,x, r, k2)
    return OneWayMesh(
        nbodies_i,
        nbodies_a,
        panel_indices_i,
        panel_indices_a,
        mesh2panel_i,
        mesh2panel_a,
        x,
        r,
        k2,
    )
end
