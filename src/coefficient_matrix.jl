#=

Additional functions needed to create one-way coefficient matrices

Authors: Judd Mehr,

=#

#---------------------------------#
#       Vortex Coefficients       #
#---------------------------------#
function assemble_one_way_coefficient_matrix(
    mesh, influence_panels::TP, affect_panels; singularity="vortex"
) where {TP<:ff.Panel}
    return assemble_one_way_coefficient_matrix(
        mesh, [influence_panels], affect_panels; singularity=singularity
    )
end

function assemble_one_way_coefficient_matrix(
    mesh, influence_panels, affect_panels::TP; singularity="vortex"
) where {TP<:ff.Panel}
    return assemble_one_way_coefficient_matrix(
        mesh, influence_panels, [affect_panels]; singularity=singularity
    )
end

function assemble_one_way_coefficient_matrix(
    mesh, influence_panels::TP, affect_panels::TP; singularity="vortex"
) where {TP<:ff.Panel}
    return assemble_one_way_coefficient_matrix(
        mesh, [influence_panels], [affect_panels]; singularity=singularity
    )
end

"""
    assemble_one_way_coefficient_matrix(mesh, influence_panels, affect_panels; kwargs)

Function similar to FLOWFoil assemble coefficient functions, but only for one-way effects.

**Arguments:**
- `mesh::OneWayMesh` : OneWayMesh object with relative geometry from influence to affected panels.
- `influence_panels::Vector{FLOWFoil.AxisymmetricPanel}` : vector of panel objects doing the influencing
- `affect_panels::Vector{FLOWFoil.AxisymmetricPanel}` : vector of panel objects being affected

Multiple Dispatch allows for single panel objects as one or both inputs as well if there is only one body influencing and/or being affected.

**Keyword Arguments:**
- `singularity::String` : selects "vortex" or "source" as the singularity for which to calculate the x values.  vortex is default.

**Returns:**
- `A::Matrix{Float}` : Aerodynamic coefficient matrix of influence on affect panels.
"""
function assemble_one_way_coefficient_matrix(
    mesh, influence_panels, affect_panels; singularity="vortex"
)

    ### --- SETUP --- ###

    # - Rename for Convenience - #
    idx_i = mesh.influence_panel_indices
    N = idx_i[end][end]

    idx_a = mesh.affect_panel_indices
    M = idx_a[end][end]

    m2p_i = mesh.mesh2panel_influence
    m2p_a = mesh.mesh2panel_affect

    # initialize coefficient matrix
    TF = eltype(mesh.m)
    amat = zeros(TF, (M, N))

    # Loop through system

    ### --- Loop through bodies --- ###
    for m in 1:(mesh.n_affect_bodies)
        for n in 1:(mesh.n_influence_bodies)
            ### --- Loop through panels --- ###
            for i in idx_a[m]
                for j in idx_i[n]

                    ### --- Calculate influence coefficient --- ###
                    if singularity == "vortex"
                        amat[i, j] = calculate_ring_vortex_influence_off_body(
                            affect_panels[m], influence_panels[n], mesh, i, j
                        )
                    elseif singularity == "source"
                        amat[i, j] = calculate_ring_source_influence_off_body(
                            affect_panels[m], influence_panels[n], mesh, i, j
                        )
                    else
                        @error "no singularity of type $(singularity)"
                    end
                end
            end
        end
    end

    return amat
end

"""
    calculate_ring_vortex_influence_off_body(paneli, panelj, mesh, i, j)

Function simiar to FLOWFoil's calculate_ring_vortex_influence function, but specifically for geometry not located on the body; used in one-way coefficient calculations.

**Arguments:**
- `paneli::FLOWFoil.AxiSymPanel` : the ith panel (the panel being influenced).
- `panelj::FLOWFoil.AxiSymPanel` : the jth panel (the panel doing the influencing).
- `mesh::OneWayMesh` : OneWayMesh object with relative geometry from influence to affected panels.
- `i::Int` : index for ith panel
- `j::Int` : index for jth panel

**Returns:**
- `aij::Float` : Influence of vortex ring strength at panel j onto panel i.
"""
function calculate_ring_vortex_influence_off_body(paneli, panelj, mesh, i, j)
    m2p_i = mesh.mesh2panel_influence
    m2p_a = mesh.mesh2panel_affect

    #calculate unit velocities
    u = ff.get_u_ring_vortex(
        mesh.x[i, j],
        mesh.r[i, j],
        panelj.panel_center[m2p_i[j], 2],
        panelj.panel_length[m2p_i[j]],
        mesh.m[i, j],
    )

    v = ff.get_v_ring_vortex(
        mesh.x[i, j], mesh.r[i, j], panelj.panel_center[m2p_i[j], 2], mesh.m[i, j]
    )

    #return appropriate strength
    # if asin(sqrt(m)) != pi / 2
    if mesh.m[i, j] != 1.0

        #panels are different
        return (
            -u * cos(paneli.panel_angle[m2p_a[i]]) + v * sin(paneli.panel_angle[m2p_a[i]])
        ) * panelj.panel_length[m2p_i[j]]
    else
        #same panel -> self induction equation

        #NOTE: this is not eqn 4.22 in Lewis.  Their code uses this expression which seems to avoid singularities better.  Not sure how they changed the second term (from dj/4piR to -R) though; perhaps the R in the text != the curvature in the code (radiusofcurvature vs curvature).

        # constant used in multiple places to clean things up
        cons = 4.0 * pi * panelj.panel_center[m2p_i[j], 2] / panelj.panel_length[m2p_i[j]]

        # return self inducement coefficient
        return -0.5 - panelj.panel_curvature[m2p_i[j]] -
               (log(2.0 * cons) - 0.25) / cons * cos(panelj.panel_angle[m2p_i[j]])
    end
end

#---------------------------------#
#       Source Coefficients       #
#---------------------------------#

"""
    calculate_ring_source_influence_off_body(paneli, panelj, mesh, i, j)

Function simiar to FLOWFoil's calculate_ring_vortex_influence function, but specifically for geometry not located on the body; used in one-way coefficient calculations, and for sources rather than vortices.

**Arguments:**
- `paneli::FLOWFoil.AxiSymPanel` : the ith panel (the panel being influenced).
- `panelj::FLOWFoil.AxiSymPanel` : the jth panel (the panel doing the influencing).
- `mesh::OneWayMesh` : OneWayMesh object with relative geometry from influence to affected panels.
- `i::Int` : index for ith panel
- `j::Int` : index for jth panel

**Returns:**
- `aij::Float` : Influence of source ring strength at panel j onto panel i.
"""
function calculate_ring_source_influence_off_body(paneli, panelj, mesh, i, j)
    m2p_i = mesh.mesh2panel_influence
    m2p_a = mesh.mesh2panel_affect

    #calculate unit velocities
    u = ff.get_u_ring_source(
        mesh.x[i, j],
        mesh.r[i, j],
        panelj.panel_center[m2p_i[j], 2],
        panelj.panel_length[m2p_i[j]],
        mesh.m[i, j],
    )

    v = ff.get_v_ring_source(
        mesh.x[i, j], mesh.r[i, j], panelj.panel_center[m2p_i[j], 2], mesh.m[i, j]
    )

    #return appropriate strength
    # if asin(sqrt(m)) != pi / 2
    if mesh.m[i, j] != 1.0

        #panels are different
        return (
            -u * sin(paneli.panel_angle[m2p_a[i]]) + v * cos(paneli.panel_angle[m2p_a[i]])
        ) * panelj.panel_length[m2p_i[j]]
    else
        #same panel -> self induction equation

        # return self inducement coefficient
        return 0.5 # +
        # (
        # -u * sin(paneli.panel_angle[m2p_a[i]]) + v * cos(paneli.panel_angle[m2p_a[i]])
        # ) * panelj.panel_length[m2p_i[j]]
    end
end
