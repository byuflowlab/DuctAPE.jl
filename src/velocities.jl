######################################################################
#                                                                    #
#                         "MESH" GENERATION                          #
#                                                                    #
######################################################################
"""

calculate "mesh" geometry without creating a mesh object
returns zeros if ring is on (or approximately on) the axis of rotation (zero radius)

**Arguments:**
- `controlpoint::Vector{Float}` [z r] coordinates of point being influenced
- `node::Vector{Float}` : [z r] coordinates of singularity ring
"""
function calculate_xrm(controlpoint, node)
    if isapprox(node[2], 0.0)
        return 0.0, 0.0, 0.0, 0.0
    else
        # normalized axial distance
        xi = (controlpoint[1] - node[1]) / node[2]

        # normalized radial distance
        rho = controlpoint[2] / node[2]

        # elliptic integral parameter
        m = (4.0 * rho) / (xi^2 + (rho + 1)^2)

        # influence point radial position
        rj = node[2]

        return xi, rho, m, rj
    end
end

######################################################################
#                                                                    #
#                        Elliptic Functions                          #
#                                                                    #
######################################################################
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

######################################################################
#                                                                    #
#                         INDUCED VELOCITIES                         #
#                                                                    #
######################################################################
#---------------------------------#
#    Unit Induced Velocities      #
#---------------------------------#

##### ----- Vortex ----- #####
"""
axial velocity induced by axisymmetric vortex ring. uses equivalent smoke ring induced velocity for self-induction, and returns zero if vortex ring is on axis of rotation (zero radius).

**Arguments:**
- `xi::Float` : normalized z-coordinate, (z-zo)/ro
- `rho::Float` : normalized r-coordinate, r/ro
- `m::Float` : Elliptic Integral parameter, 4rho/sqrt(xi^2+(rho+1)^2)
- `r_influence::Float` : radial location of vortex ring, ro
- `influence_length::Float` : length of panel used in calculating self-induction

**Returns:**
- `vz::Float` : axially induced velocity of vortex ring
"""
function vortex_ring_vz(xi, rho, m, r_influence, influence_length)

    # check panel locations
    if abs(r_influence)<=eps()
        # if influence on the axis, the influence is set to zero
        return 0.0
    elseif (xi^2 + (rho - 1.0)^2 <= eps())
        # set self-induced case is "smoke ring" self influence in axial direction only.
        return smoke_ring_vz(r_influence, influence_length)
    else
        #get the first denominator
        den1 = 2.0 * pi * r_influence * sqrt(xi^2 + (rho + 1.0)^2)

        #get numerator and denominator of second fraction
        num2 = 2.0 * (rho - 1.0)
        den2 = xi^2 + (rho - 1.0)^2

        #get values for elliptic integrals
        K, E = get_elliptics(m)

        # negative in Lewis' version is due to Lewis' convention that the vortex is postive clockwise
        # return -1.0 / den1 * (K - (1.0 + num2 / den2) * E)
        # No negative in your derivation with vortex postive according to right hand coordinate system
        # TODO: this may affect coupling with rotor model, need to check signs there.
        return 1.0 / den1 * (K - (1.0 + num2 / den2) * E)
    end
end

"""
equivalent "smoke" ring self-induced velocity
"""
function smoke_ring_vz(r_influence, influence_length)
    # return -1.0 / (4.0 * pi * r_influence) * (log(8.0 * pi * r_influence / influence_length) - 0.25)
    # Lamb has negative out front due to vortex in opposite direction to you
    return 1.0 / (4.0 * pi * r_influence) * (log(8.0 * pi * r_influence / influence_length) - 0.25)
end


"""
radial velocity induced by axisymmetric vortex ring. returns zero if vortex ring is on axis of rotation (zero radius), the point of influence is on the axis, or if self-inducing velocity.

**Arguments:**
- `xi::Float` : normalized z-coordinate, (z-zo)/ro
- `rho::Float` : normalized r-coordinate, r/ro
- `m::Float` : Elliptic Integral parameter, 4rho/sqrt(xi^2+(rho+1)^2)
- `r_influence::Float` : radial location of vortex ring, ro

**Returns:**
- `vr::Float` : radially induced velocity of vortex ring
"""
function vortex_ring_vr(xi, rho, m, r_influence)

    # return 0.0 for self-induced, influence on axis, or target on axis cases
    if (xi^2 + (rho - 1.0)^2 <= eps()) || abs(r_influence)<=eps() || isapprox(rho, 0.0)
        return 0.0
    else
        #get numerator and denominator of first fraction
        num1 = xi / rho
        den1 = 2.0 * pi * r_influence * sqrt(xi^2 + (rho + 1.0)^2)

        #get numerator and denominator of second fraction
        num2 = 2.0 * rho
        den2 = xi^2 + (rho - 1.0)^2

        #get values for elliptic integrals
        K, E = get_elliptics(m)

        # positive is what lewis had using Gamma in opposite direction to you
        # return num1 / den1 * (K - (1.0 + num2 / den2) * E)
        # negative sign is what you got in your derivation
        return -num1 / den1 * (K - (1.0 + num2 / den2) * E)
    end
end

##### ----- Source ----- #####

"""
axial velocity induced by axisymmetric source ring. returns zero if source ring is on axis of rotation (zero radius), the point of influence is on the axis, or if self-inducing velocity.

**Arguments:**
- `xi::Float` : normalized z-coordinate, (z-zo)/ro
- `rho::Float` : normalized r-coordinate, r/ro
- `m::Float` : Elliptic Integral parameter, 4rho/sqrt(xi^2+(rho+1)^2)
- `r_influence::Float` : radial location of vortex ring, ro

**Returns:**
- `vz::Float` : axially induced velocity of source ring
"""
function source_ring_vz(xi, rho, m, r_influence)

    # return zero for the self-induced off body case
    if (xi^2 + (rho - 1.0)^2 <= eps()) || isapprox(r_influence, 0.0) || isapprox(rho, 0.0)
        return 0.0
    else

        #get values for elliptic integrals
        K, E = get_elliptics(m)

        #get the first denominator
        den1 = 2.0 * pi * r_influence * sqrt(xi^2 + (rho + 1.0)^2)

        #get numerator and denominator of second fraction
        num2 = 2 * xi * E
        den2 = xi^2 + (rho - 1)^2

        return 1.0 / den1 * (num2 / den2)
    end
end

"""
radial velocity induced by axisymmetric source ring. returns zero if source ring is on axis of rotation (zero radius), the point of influence is on the axis, or if self-inducing velocity.

**Arguments:**
- `xi::Float` : normalized z-coordinate, (z-zo)/ro
- `rho::Float` : normalized r-coordinate, r/ro
- `m::Float` : Elliptic Integral parameter, 4rho/sqrt(xi^2+(rho+1)^2)
- `r_influence::Float` : radial location of vortex ring, ro

**Returns:**
- `vr::Float` : radially induced velocity of source ring
"""
function source_ring_vr(xi, rho, m, r_influence)

    # return zero for the self-induced off-body case
    if (xi^2 + (rho - 1.0)^2 <= eps()) || isapprox(r_influence, 0.0) || isapprox(rho, 0.0)
        return 0.0
    else

        #get values for elliptic integrals
        K, E = get_elliptics(m)

        #get numerator and denominator of first fraction
        num1 = 1.0 / rho
        den1 = 2.0 * pi * r_influence * sqrt(xi^2 + (rho + 1.0)^2)

        #get numerator and denominator of second fraction
        num2 = 2 * rho * (rho - 1.0)
        den2 = xi^2 + (rho - 1)^2

        return num1 / den1 * (K - (1.0 - num2 / den2) * E)
    end
end

#---------------------------------#
#   "Panel" Induced Velocities    #
#---------------------------------#

##### ----- Vortex ----- #####
"""
in place calculation of axial and radial components of induced velocity due to an axisymmetric vortex ring

**Arguments:**
- `vel::Vector{Float}` : [vz, vr] vector of induced velocity components
- `controlpoint::Vector{Float}` [z r] coordinates of point being influenced
- `node1::Vector{Float}` : [z r] coordinates of first node
- `node2::Vector{Float}` : [z r] coordinates of second node
- `influence_length::Float` : length over which vortex ring influence is applied on the surface.
- `gamma::Vecto{Float}` : panel node strengths, default = [1.0;1.0] (unit strength)
"""
function vortex_induced_velocity!(
        vel, controlpoint, node1, node2, influence_length, gamma=[1.0; 1.0]
)

    #TODO: replace with integration stuff
vel .= nominal_vortex_panel_integration(
    node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0], debug=false
)
    #TODO: how to determine self-induced panel? check if controlpoint is on node line?
    # get full induced velocities of vortex ring
    vel[1] +=
        vortex_ring_vz(xi, rho, m, rj, influence_length) * #length input needed here in case of self-induced case
        gamma *
        influence_length
    vel[2] += vortex_ring_vr(xi, rho, m, rj) * gamma * influence_length

    return vel
end

"""
out of place calculation of axial and radial components of induced velocity due to an axisymmetric vortex ring

**Arguments:**
- `controlpoint::Vector{Float}` [z r] coordinates of point being influenced
- `node::Vector{Float}` : [z r] coordinates of vortex ring
- `influence_length::Float` : length over which vortex ring influence is applied on the surface.
- `gamma::Float` : vortex constant circulation value, default = 1.0 (unit vortex)

**Returns:**
- `vel::Vector{Float}` : [vz, vr] vector of induced velocity components
"""
function vortex_induced_velocity(
    controlpoint::AbstractVector{T3},
    node::AbstractVector{T1},
    influence_length::T2,
    gamma::T4=1.0,
) where {T1,T2,T3,T4}

    # Initialize
    T = promote_type(T1, T2, T3, T4)
    vel = zeros(T, 2)

    # get velocities
    vortex_induced_velocity!(
        vel, controlpoint, node, influence_length, gamma
    )

    return vel
end

"""
out of place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric vortex rings

Used for constructing the influence matrices from the body to the rotor/wake in order to calculate the body induced velocities on the rotor/wake.

**Arguments:**
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `node::Matrix{Float}` : [z r] coordinates of vortex rings
- `influence_length::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
- `gamma::Vector{Float}` : vortex constant circulation values

**Returns:**
- `AIC::Array{Float}` : N-controlpoint x N-node x [vz, vr] array of induced velocity components
"""
function influencefromvortices(
    controlpoints::AbstractMatrix{T1},
    nodes::AbstractMatrix{T2},
    influence_lengths::AbstractVector{T3},
    strengths::AbstractVector{T4},
) where {T1,T2,T3,T4}

# Initialize
    T = promote_type(T1, T2, T3, T4)
    AIC = zeros(T, size(controlpoints, 1), size(nodes, 1), 2)

    influencefromvortices!(
        AIC, controlpoints, nodes, influence_lengths, strengths
    )

    return AIC
end

"""
in place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric vortex rings

Used for constructing the influence matrices from the body to the rotor/wake in order to calculate the body induced velocities on the rotor/wake.

**Arguments:**
- `AIC::Array{Float}` : N-controlpoint x N-node x [vz, vr] array of induced velocity components
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `node::Matrix{Float}` : [z r] coordinates of vortex rings
- `influence_length::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
- `gamma::Vector{Float}` : vortex constant circulation values
"""
function influencefromvortices!(
    AIC, controlpoints, nodes, influence_lengths, strengths
)

# loop through control points
    for (i, cpi) in enumerate(eachrow(controlpoints))
        # loop through panels doing the influencing
        for (j, (gamma, nj, lj)) in
            enumerate(zip(strengths, eachrow(nodes), influence_lengths))

            # get unit induced velocity from the panel onto the control point
            vortex_induced_velocity!(view(AIC, i, j, :), cpi, nj, lj, gamma)
        end
    end

    return AIC
end

"""
out of place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric vortex rings but only returning sum of velocities on each control point

Used in calculating velocities on body surfaces in preparation to obtain tangential components and eventually pressure distributions.
Note: there is probably a more efficient way to achieve this functionality.

**Arguments:**
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `node::Matrix{Float}` : [z r] coordinates of vortex rings
- `influence_length::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
- `gamma::Vector{Float}` : vortex constant circulation values

**Returns:**
- `V::Array{Float}` : N-controlpoint x [vz, vr] array of summed induced velocity components due to vortex rings
"""
function vfromvortices(
    controlpoints::AbstractMatrix{T1},
    nodes::AbstractMatrix{T2},
    influence_lengths::AbstractVector{T3},
    strengths::AbstractArray{T4},
) where {T1,T2,T3,T4}

# Initialize
    T = promote_type(T1, T2, T3, T4)
    V = zeros(T, size(controlpoints, 1), 2)

    vfromvortices!(V, controlpoints, nodes, influence_lengths, strengths)

    return V
end

"""
In place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric vortex rings but only returning sum of velocities on each control point

Used in calculating velocities on body surfaces in preparation to obtain tangential components and eventually pressure distributions.
Note: there is probably a more efficient way to achieve this functionality.

**Arguments:**
- `V::Array{Float}` : N-controlpoint x [vz, vr] array of summed induced velocity components due to vortex rings
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `node::Matrix{Float}` : [z r] coordinates of vortex rings
- `influence_length::Vector{Float}` : lengths over which vortex ring influence is applied on the surface.
- `gamma::Vector{Float}` : vortex constant circulation values
"""
function vfromvortices!(
        V,
    controlpoints,
    nodes,
    influence_lengths,
    strengths,
)

# Loop through control points
    for (i, (cpi, vel)) in enumerate(zip(eachrow(controlpoints), eachrow(V)))
        # loop through panels doing the influencing
        for (j, (gamma, nj, lj)) in
            enumerate(zip(strengths, eachrow(nodes), influence_lengths))

            # get unit induced velocity from the panel onto the control point
            vortex_induced_velocity!(vel, cpi, nj, lj, gamma)
        end
    end

    return nothing
end

##### ----- Source ----- #####
"""
in place calculation of axial and radial components of induced velocity due to an axisymmetric source ring

**Arguments:**
- `vel::Vector{Float}` : [vz, vr] vector of induced velocity components
- `controlpoint::Vector{Float}` [z r] coordinates of point being influenced
- `node::Vector{Float}` : [z r] coordinates of source ring
- `influence_length::Float` : length over which source ring influence is applied on the surface.
- `sigma::Float` : source constant circulation value, default = 1.0 (unit source)
"""
function source_induced_velocity!(
    vel, controlpoint, node, influence_length, sigma=1.0
)

    # get relative geometry
    xi, rho, m, rj = calculate_xrm(controlpoint, node)

    # get unit induced velocities of source ring
    vel[1] += source_ring_vz(xi, rho, m, rj) * sigma * influence_length
    vel[2] += source_ring_vr(xi, rho, m, rj) * sigma * influence_length

    return nothing
end

"""
out of place calculation of axial and radial components of induced velocity due to an axisymmetric source ring

**Arguments:**
- `controlpoint::Vector{Float}` [z r] coordinates of point being influenced
- `node::Vector{Float}` : [z r] coordinates of source ring
- `influence_length::Float` : length over which source ring influence is applied on the surface.
- `sigma::Float` : source constant circulation value, default = 1.0 (unit source)

**Returns:**
- `vel::Vector{Float}` : [vz, vr] vector of induced velocity components
"""
function source_induced_velocity(
    controlpoint::AbstractVector{T3},
    node::AbstractVector{T1},
    influence_length::T2,
    sigma::T4=1.0,
) where {T1,T2,T3,T4}

    # Initialize
    T = promote_type(T1, T2, T3, T4)
    vel = zeros(T, 2)

    # get velocities
    source_induced_velocity!(
        vel, controlpoint, node, influence_length, sigma
    )

    return vel
end

"""
out of place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric source rings

Used for constructing the influence matrices from the rotor to the body/wake in order to calculate the rotor profile drag induced velocities on the body/wake.

**Arguments:**
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `node::Matrix{Float}` : [z r] coordinates of source rings
- `influence_length::Vector{Float}` : lengths over which source ring influence is applied on the surface.
- `gamma::Vector{Float}` : source constant circulation values

**Returns:**
- `AIC::Array{Float}` : N-controlpoint x N-node x [vz, vr] array of induced velocity components
"""
function influencefromsources(
    controlpoints::AbstractMatrix{T1},
    nodes::AbstractArray{T2},
    influence_lengths::AbstractVector{T3},
    strengths::AbstractArray{T4},
) where {T1,T2,T3,T4}

# Initialize
    T = promote_type(T1, T2, T3, T4)
    AIC = zeros(T, size(controlpoints, 1), size(nodes, 1), 2)

    influencefromsources!(
        AIC, controlpoints, nodes, influence_lengths, strengths
    )

    return AIC
end

"""
in place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric source rings

Used for constructing the influence matrices from the rotor to the body/wake in order to calculate the rotor profile drag induced velocities on the body/wake.

**Arguments:**
- `AIC::Array{Float}` : N-controlpoint x N-node x [vz, vr] array of induced velocity components
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `node::Matrix{Float}` : [z r] coordinates of source rings
- `influence_length::Vector{Float}` : lengths over which source ring influence is applied on the surface.
- `gamma::Vector{Float}` : source constant circulation values
"""
function influencefromsources!(
    AIC, controlpoints, nodes, influence_length, strengths
)
    for (i, cpi) in enumerate(eachrow(controlpoints))
        # loop through panels doing the influencing
        for (j, (sigma, nj, lj)) in
            enumerate(zip(strengths, eachrow(nodes), influence_length))

            # get unit induced velocity from the panel onto the control point
            source_induced_velocity!(view(AIC, i, j, :), cpi, nj, lj, sigma)
        end
    end

    return nothing
end

"""
out of place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric source rings but only returning sum of velocities on each control point

Used in calculating velocities on body surfaces in preparation to obtain tangential components and eventually pressure distributions.
Note: there is probably a more efficient way to achieve this functionality.

**Arguments:**
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `node::Matrix{Float}` : [z r] coordinates of source rings
- `influence_length::Vector{Float}` : lengths over which source ring influence is applied on the surface.
- `sigma::Vector{Float}` : source constant circulation values

**Returns:**
- `V::Array{Float}` : N-controlpoint x [vz, vr] array of summed induced velocity components due to source rings
"""
function vfromsources(
    controlpoints::AbstractMatrix{T1},
    nodes::AbstractArray{T2},
    influence_lengths::AbstractVector{T3},
    strengths::AbstractArray{T4},
) where {T1,T2,T3,T4}

# Initialize
    T = promote_type(T1, T2, T3, T4)
    V = zeros(T, size(controlpoints, 1), 2)

    vfromsources!(V, controlpoints, nodes, influence_lengths, strengths)

    return V
end

"""
in place calculation of axial and radial components of induced velocity for a set of control points due to a set of axisymmetric source rings but only returning sum of velocities on each control point

Used in calculating velocities on body surfaces in preparation to obtain tangential components and eventually pressure distributions.
Note: there is probably a more efficient way to achieve this functionality.

**Arguments:**
- `V::Array{Float}` : N-controlpoint x [vz, vr] array of summed induced velocity components due to source rings
- `controlpoint::Matrix{Float}` [z r] coordinates of points being influenced
- `node::Matrix{Float}` : [z r] coordinates of source rings
- `influence_length::Vector{Float}` : lengths over which source ring influence is applied on the surface.
- `sigma::Vector{Float}` : source constant circulation values
"""
function vfromsources!(
    V, controlpoints, nodes, influence_length, strengths
)
    for (i, (cpi, vel)) in enumerate(zip(eachrow(controlpoints), eachrow(V)))
        # loop through panels doing the influencing
        for (j, (sigma, nj, lj)) in
            enumerate(zip(strengths, eachrow(nodes), influence_length))

            # get unit induced velocity from the panel onto the control point
            source_induced_velocity!(vel, cpi, nj, lj, sigma)
        end
    end

    return nothing
end

