#=

Additional functions needed to create one-way coefficient matrices

Authors: Judd Mehr,

=#

#= TODO:
need to redo everything here, almost.  We aren't actually wanting the vortex coefficients, but rather just the induced velocities due to a unit singularity which we can multiple by the singularities to get the velocities at the given points.  still keep things as matricies, but we'll want to return a matrix for the axial components and a second matrix for the radial components.  Also need to alter the induced velocity expressions to match Ryall and Collins for these cases where we aren't working with a panel method (probably). BUT need to check a single panel case to make sure that all the velocity directions are correct.
=#

#---------------------------------#
#       Vortex Coefficients       #
#---------------------------------#
function assemble_induced_velocity_matrices(
    mesh, influence_panels::TP, affect_panels; singularity="vortex"
) where {TP<:ff.Panel}
    return assemble_induced_velocity_matrices(
        mesh, [influence_panels], affect_panels; singularity=singularity
    )
end

function assemble_induced_velocity_matrices(
    mesh, influence_panels, affect_panels::TP; singularity="vortex"
) where {TP<:ff.Panel}
    return assemble_induced_velocity_matrices(
        mesh, influence_panels, [affect_panels]; singularity=singularity
    )
end

function assemble_induced_velocity_matrices(
    mesh, influence_panels::TP, affect_panels::TP; singularity="vortex"
) where {TP<:ff.Panel}
    return assemble_induced_velocity_matrices(
        mesh, [influence_panels], [affect_panels]; singularity=singularity
    )
end

"""
    assemble_induced_velocity_matrices(mesh, influence_panels, affect_panels; kwargs)

Function similar to FLOWFoil assemble coefficient functions, but only for one-way effects.

**Arguments:**
- `mesh::OneWayMesh` : OneWayMesh object with relative geometry from influence to affected panels.
- `influence_panels::Vector{FLOWFoil.AxisymmetricPanel}` : vector of panel objects doing the influencing
- `affect_panels::Vector{FLOWFoil.AxisymmetricPanel}` : vector of panel objects being affected

Multiple Dispatch allows for single panel objects as one or both inputs as well if there is only one body influencing and/or being affected.

**Keyword Arguments:**
- `singularity::String` : selects "vortex" or "source" as the singularity for which to calculate the x values.  vortex is default.

**Returns:**
- `vxdmat::Matrix{Float}` : v_x_ij * d_j for all i, j
- `vrdmat::Matrix{Float}` : v_r_ij * d_j for all i, j
"""
function assemble_induced_velocity_matrices(
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
    vxdmat = zeros(TF, (M, N))
    vrdmat = zeros(TF, (M, N))

    # Loop through system

    ### --- Loop through bodies --- ###
    for m in 1:(mesh.n_affect_bodies)
        for n in 1:(mesh.n_influence_bodies)
            ### --- Loop through panels --- ###
            for i in idx_a[m]
                for j in idx_i[n]

                    ### --- Calculate influence coefficient --- ###
                    if singularity == "vortex"
                        vxdmat[i, j], vrdmat[i, j] = calculate_ring_vortex_influence_off_body(
                            affect_panels[m], influence_panels[n], mesh, i, j
                        )
                    elseif singularity == "source"
                        vxdmat[i, j], vrdmat[i, j] = calculate_ring_source_influence_off_body(
                            affect_panels[m], influence_panels[n], mesh, i, j
                        )
                    else
                        @error "no singularity of type $(singularity)"
                    end
                end
            end
        end
    end

    return [[vxdmat] [vrdmat]]
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
    vx = get_vx_ring_vortex_off_body(
        mesh.x[i, j], mesh.r[i, j], panelj.panel_center[m2p_i[j], 2], mesh.m[i, j]
    )

    vr = get_vr_ring_vortex_off_body(
        mesh.x[i, j], mesh.r[i, j], panelj.panel_center[m2p_i[j], 2], mesh.m[i, j]
    )

    return vx * panelj.panel_length[j], vr * panelj.panel_length[j]
end

"""
    get_vx_ring_vortex_off_body(x, r, rj, m)

Calculate x-component of velocity influence of vortex ring.

**Arguments:**
- `x::Float` : ratio of difference of ith and jth panel x-locations and jth panel r-location ( (xi-xj)/rj )
- `r::Float` : ratio of r-locations of ith and jth panels (ri/rj)
- `rj::Float` : r-location of the jth panel control point
- `m::Float` : Elliptic Function parameter

**Returns:**
- `vxij::Float` : x-component of velocity induced by panel j onto panel i
"""
function get_vx_ring_vortex_off_body(x, r, rj, m; probe=false)

    #get the first denominator
    den1 = 2.0 * pi * rj * sqrt(x^2 + (r + 1.0)^2)

    #get numerator and denominator of second fraction
    num2 = 2 * (r - 1)
    den2 = x^2 + (r - 1)^2

    #get values for elliptic integrals
    K, E = get_elliptics(m)

    # set self-induced case to zero for the off-body case
    if x^2 + (r - 1.0)^2 <= eps()
        return 0.0
    else
        # negative here is due to our convention that the vortex is postive clockwise
        return -1.0 / den1 * (K - (1.0 + num2 / den2) * E)
    end
end

"""
    get_vr_ring_vortex_off_body(x, r, rj, m)

Calculate r-component of velocity influence of vortex ring.

**Arguments:**
- `x::Float` : ratio of difference of ith and jth panel x-locations and jth panel r-location ( (xi-xj)/rj )
- `r::Float` : ratio of r-locations of ith and jth panels (ri/rj)
- `rj::Float` : r-location of the jth panel control point
- `m::Float` : Elliptic Function parameter

**Returns:**
- `vrij::Float` : r-component of velocity induced by panel j onto panel i
"""
function get_vr_ring_vortex_off_body(x, r, rj, m; probe=false)

    #get numerator and denominator of first fraction
    num1 = x / r
    den1 = 2.0 * pi * rj * sqrt(x^2 + (r + 1.0)^2)

    num2 = 2 * r
    den2 = x^2 + (r - 1)^2

    #get values for elliptic integrals
    K, E = get_elliptics(m)

    # set self-induced case to zero for the off-body case
    if x^2 + (r - 1.0)^2 <= eps()
        return 0.0
    else
        return num1 / den1 * (K - (1.0 + num2 / den2) * E)
    end
end

"""
    get_elliptics(m)

Calculate value of elliptic functions for the given geometry parameter.

**Arguments:**
- `m::Float` : Elliptic Function parameter

**Returns:**
- `K::Float` : K(m), value of elliptic function of the first kind at m.
- `E::Float` : E(m), value of eeliptic function of the second kind at m.
"""
function get_elliptics(m)
    if m > 1 || isnan(m)
        #m cannot be greater than 1 for elliptic functions, and cannot mathematically be either, but numerically might be infinitesimally larger.
        m = 1.0
    end
    return SpecialFunctions.ellipk(m), SpecialFunctions.ellipe(m)
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
    vx = get_vx_ring_source_off_body(
        mesh.x[i, j],
        mesh.r[i, j],
        panelj.panel_center[m2p_i[j], 2],
        panelj.panel_length[m2p_i[j]],
        mesh.m[i, j],
    )

    vr = get_vr_ring_source_off_body(
        mesh.x[i, j], mesh.r[i, j], panelj.panel_center[m2p_i[j], 2], mesh.m[i, j]
    )

    return vx * panelj.panel_length[j], vr * panelj.panel_length[j]
end

"""
    get_vx_ring_source_off_body(x, r, rj, m)

Calculate x-component of velocity influence of source ring.

**Arguments:**
- `x::Float` : ratio of difference of ith and jth panel x-locations and jth panel r-location ( (xi-xj)/rj )
- `r::Float` : ratio of r-locations of ith and jth panels (ri/rj)
- `rj::Float` : r-location of the jth panel control point
- `m::Float` : Elliptic Function parameter

**Returns:**
- `uij::Float` : x-component of velocity induced by panel j onto panel i
"""
function get_vx_ring_source_off_body(x, r, rj, dj, m; probe=false)

    #get values for elliptic integrals
    K, E = get_elliptics(m)

    #get the first denominator
    den1 = 2.0 * pi * rj * sqrt(x^2 + (r + 1.0)^2)

    #get numerator and denominator of second fraction
    num2 = 2 * x * E
    den2 = x^2 + (r - 1)^2

    # return zero for the self-induced off body case
    if x^2 + (r - 1.0)^2 <= eps()
        return 0.0
    else
        return 1.0 / den1 * (num2 / den2)
    end
end

"""
    get_vr_ring_source_off_body(x, r, rj, m)

Calculate r-component of velocity influence of source ring.

**Arguments:**
- `x::Float` : ratio of difference of ith and jth panel x-locations and jth panel r-location ( (xi-xj)/rj )
- `r::Float` : ratio of r-locations of ith and jth panels (ri/rj)
- `rj::Float` : r-location of the jth panel control point
- `m::Float` : Elliptic Function parameter

**Returns:**
- `vij::Float` : r-component of velocity induced by panel j onto panel i
"""
function get_vr_ring_source_off_body(x, r, rj, m; probe=false)

    #get values for elliptic integrals
    K, E = get_elliptics(m)

    #get numerator and denominator of first fraction
    num1 = 1.0 / r
    den1 = 2.0 * pi * rj * sqrt(x^2 + (r + 1.0)^2)

    #get numerator and denominator of second fraction
    num2 = 2 * r * (r - 1.0)
    den2 = x^2 + (r - 1)^2

    # return zero for the self-induced off-body case
    if x^2 + (r - 1.0)^2 <= eps()
        return 0.0
    else
        return num1 / den1 * (K - (1.0 - num2 / den2) * E)
    end
end
