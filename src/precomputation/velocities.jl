#=
Functions for obtaining velocities induced by unit strength rings and arbitrary strength linear panels
=#

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
#  Unit Ring Induced Velocities   #
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
    if abs(r_influence) <= eps()
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
    return 1.0 / (4.0 * pi * r_influence) *
           (log(8.0 * pi * r_influence / influence_length) - 0.25)
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
    if (xi^2 + (rho - 1.0)^2 <= eps()) || abs(r_influence) <= eps() || isapprox(rho, 0.0)
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
    if (xi^2 + (rho - 1.0)^2 <= eps()) || abs(r_influence) < eps() || abs(rho) < eps()
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