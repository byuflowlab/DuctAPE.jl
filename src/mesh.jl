#=

Functions for generating meshes (relational geometry between objects)

Authors: Judd Mehr,

=#

"""
    OneWayMesh

**Fields:**
- `n_influence_bodies::Int` : number of bodies of influence, e.g. 2 for the duct and hub as part of the body system.
- `n_affect_bodies::Int` : number of bodies being affected, e.g. 1 for the first rotor.
- `influence_panel_indices::Vector{StepRange{Int64, Int64}}` : vector of index ranges of the panels for each of the bodies of influence.
- `affect_panel_indices::Vector{StepRange{Int64, Int64}}` : vector of index ranges of the panels for each of the affected bodies.
- `mesh2panel_influence::Vector{Int}` : vector mapping the index of the full matrix to the indices of the individual bodies of influence
- `mesh2panel_affect:: Vector{Int}` : vector mapping the index of the full matrix to the indices of the individual affected bodies
- `x::Matrix{Float}` : x values used in coefficient and velocity calculations
- `r::Matrix{Float}` : r values used in coefficient and velocity calculations
- `m::Matrix{Float}` : m values used in coefficient and velocity calculations
"""
struct OneWayMesh{TF}
    n_influence_bodies::Int64
    n_affect_bodies::Int64
    influence_panel_indices::Vector{UnitRange{Int64}}
    affect_panel_indices::Vector{UnitRange{Int64}}
    mesh2panel_influence::Vector{Int64}
    mesh2panel_affect::Vector{Int64}
    x::Matrix{TF}
    r::Matrix{TF}
    m::Matrix{TF}
end

function generate_one_way_mesh(influence_panels::TP, affect_panels) where {TP<:ff.Panel}
    return generate_one_way_mesh([influence_panels], affect_panels; singularity=singularity)
end

function generate_one_way_mesh(influence_panels, affect_panels::TP) where {TP<:ff.Panel}
    return generate_one_way_mesh(influence_panels, [affect_panels]; singularity=singularity)
end

function generate_one_way_mesh(influence_panels::TP, affect_panels::TP) where {TP<:ff.Panel}
    return generate_one_way_mesh([influence_panels], [affect_panels])
end

"""
    generate_one_way_mesh(influence_panels, affect_panels; kwargs)

Function similar to FLOWFoil's meshing function, but only creates the mesh one way rather than both ways (i.e. doesn't loop through both inputs as both sources and targets).

**Arguments:**
- `influence_panels::Vector{FLOWFoil.AxisymmetricPanel}` : vector of panel objects doing the influencing
- `affect_panels::Vector{FLOWFoil.AxisymmetricPanel}` : vector of panel objects being affected

Multiple Dispatch allows for single panel objects as one or both inputs as well if there is only one body influencing and/or being affected.

**Returns:**
- `mesh::OneWayMesh` : OneWayMesh object with relative geometry from influence to affected panels.
"""
function generate_one_way_mesh(influence_panels, affect_panels)

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
                    # TODO: suspect this should be the same for both cases, Lewis seems to have a typo normalizing x by ri instead of rj.
                    x[i, j] = (xi - xj) / rj
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
