#=
in order to generate all the required coefficient matrices, we need to first make sure each geometric item has an associated panels object.
we then need to assemble the required meshes.
we can then put together the coefficient matrices.

needed matrices include:
- body to body (default in FLOWFoil), for system linear solve
- wake to body, for full system linear solve
- rotor to body, for full system linear solve
- body to rotor, for induced velocity calculation
- wake to rotor, for induced velocity calculation
- rotor to rotor, for induced velocity calculation

we only need flowfoil boundary conditions for the body to body case, so we don't want to use the system generation functions from flowfoil, but rather just the vortex and source coefficient functions.

=#

"""
"""
function assemble_one_way_vortex_matrix(mesh)
    return nothing
end

function assemble_one_way_source_matrix(mesh, influence_panels, affect_panels)

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
                    #TODO: need to add or re-write this function to take in i and j already adjusted for mesh2panel.
                    amat[i, j] = ff.calculate_ring_vortex_influence_off_body(
                        affect_panels[m], influence_panels[n], mesh, i, j
                    )
                end
            end
        end
    end

    return amat
end

"""
"""
function calculate_ring_vortex_influence_off_body(paneli, panelj, mesh, i, j)
    m2p_i = mesh.mesh2panel_influence
    m2p_a = mesh.mesh2panel_affect

    #calculate unit velocities
    u = ff.get_u_ring(
        mesh.x[i, j],
        mesh.r[i, j],
        panelj.panel_center[m2p_i[j], 2],
        panelj.panel_length[m2p_i[j]],
        mesh.m[i, j],
    )

    v = ff.get_v_ring(
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
