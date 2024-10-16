#---------------------------------#
#         Setup functions         #
#---------------------------------#

"""
    arc_lengths_from_panel_lengths(duct_panel_lengths, bl_ids)

Cumulative sum of panel lengths for the given section of surface associated with the upper or lower boundary layer.

# Arguments:
- `duct_panel_lengths::Vector{Float}` : vector of panel lengths (called influence_length in body_vortex_panels) associated with the duct (casing + nacelle).

# Returns:
- `s::Vector{Float}` : cumulative sum of panel lengths between control points in the given index range, starting from zero.
"""
function arc_lengths_from_panel_lengths(duct_panel_lengths)
    return cumsum(
        [
            0.0
            [
                0.5 * (duct_panel_lengths[i] + duct_panel_lengths[i - 1]) for
                i in 2:length(duct_panel_lengths)
            ]
        ],
    )
end

"""
   split_at_stagnation_point(duct_panel_lengths, cp_duct)

Find the stagnation point as the point with local cp closest to 1.0.

# Arguments:
- `duct_panel_lengths::Vector{Float}` : Vector of panel lengths for the duct from casing trailing edge clockwise to nacelle trailing edge.
- `cp_duct::Vector{Float}` : Vector of surface pressure coefficients on duct from casing trailing edge clockwise to nacelle trailing edge.

# Returns:
- `s::Vector{Float}` : cumulative sum of distance along panels between control points starting at 0.
- `s_stagnation::Float` : surface length at which the stagnation point occurs, calculated using extrapolation technique described below
- `lower_length::Float` : length from stagnation point to casing trailing edge control point
- `upper_length::Float` : length from stagnation point to nacelle trailing edge control point


To determine point of intersection between extrapolated lines, we do the following:

Determine the point-slope form of the lines extrapolated from either side of the panel with a pressure coefficient closest to 1.0, and set up equations to find a point somewhere along each line

intersection    Known points on    unknown    slopes of
   point         each line         factors    each line
    P1 =       [s[i-2]; cp[i-2]] + x1    *   [ds1; dcp1]
    P2 =       [s[i+1]; cp[i+1]] + x2    *   [ds2; dcp2]

where if we set P1=P2, in other words, to be the point the lines intersect, we can solve for the unknown factors that yield the point of intersection.

We therefore set P1=P2 and assemble a system of linear equations, then plug one of the solved facotrs back into its associated equation to obtain the s location of intersection of the extrapolated lines and take this to be the stagnation point, regardless of what the pressure coefficient ends up being.
"""
function split_at_stagnation_point(duct_panel_lengths, n_panels_casing)
    # function split_at_stagnation_point(duct_panel_lengths, cp_duct)

    # get surface length along entire surface
    s_upper = arc_lengths_from_panel_lengths(duct_panel_lengths[n_panels_casing:end])
    s_lower = arc_lengths_from_panel_lengths(duct_panel_lengths[n_panels_casing:-1:1])


    # # find index of cp closest to 1.0
    # _, spi = findmin(abs.(1.0 .- cp_duct))
    # # println("cp~1 at ", spi, " of ", length(duct_panel_lengths))

    # # - Extrapolate to determine smooth stagnation point along surface length - #
    # # set up system of equations to determine intersection of extrapolated lines on either side of the index closest to cp = 1.0
    # if (length(s) - spi < 2) || (spi < 2)
    #     # if stagnation point is somehow at the trailing edge just return the location
    #     # this should never happen though.
    #     # If it does, you may get non-smooth behavior,
    #     # but at least it won't throw an error.
    #     s_stagnation = s[spi]
    # else
    #     A = [
    #         s[spi - 1]-s[spi - 2] s[spi + 1]-s[spi + 2]
    #         cp_duct[spi - 1]-cp_duct[spi - 2] cp_duct[spi + 1]-cp_duct[spi + 2]
    #     ]
    #     b = [s[spi + 1] - s[spi - 2]; cp_duct[spi + 1] - cp_duct[spi - 2]]

    #     if (det(A) < eps()) || (length(s) - spi < 2)
    #         # singular matrix, lines are parallel
    #         # note: this should not happen either,
    #         # and you may get non-smooth behavior,
    #         # but at least it won't throw an error.
    #         s_stagnation = s[spi]
    #     else
    #         # solve linear system for unknown factors
    #         x = ImplicitAD.linear_solve(A, b)
    #         # get stagnation point from one of the equations
    #         s_stagnation = s[spi - 2] + x[1] * (s[spi - 1] - s[spi - 2])
    #     end
    # end



    return s_upper, s_lower
end

"""
    bl_step_fun(n, m, p)

Function used in determining step sizes for boundary layer calculation. f(n) = m*n^p

Given a number of steps, n ∈ [1:N], provides the cumulative step lengths according to the power, p, and the multiplicative factor, m; where p determined from the `set_boundary_layer_steps` functions.
"""
function bl_step_fun(n, m, p)
    return m .* n .^ p
end

"""
    set_boundary_layer_steps(N::Int, first_step_size, total_length)

Sets boundary layer steps based on desired number of steps (must be an Integer), an initial step size, and the total cumulative length of the steps.

# Arguments:
- `N::Int` : Number of steps to take
- `first_step_size::Float` : size of first step (which is `m` in `bl_step_fun`)
- `total_length::Float` : total surface length to divide up.

# Returns:
- `steps::Vector{Float}` : steps along surface length satisfying the equation: f(n) = m*n^p with the condition that `m` is the first step size and f(N) = `total_length`
"""
function set_boundary_layer_steps(N::Int, first_step_size, total_length)

    # println("N steps: ", N)
    # println("first step size: ", first_step_size)
    # println("total length: ", ForwardDiff.value(total_length))
    # println()

    # residual
    function res(y, x, p)
        total_length = x[1]
        first_step_size = x[2]
        N = x[3]
        return total_length - first_step_size * N^y
    end

    # solver (uses Roots.jl, for which the docs indicate implicit differentiation is required)
    function solvewrap(x, p)
        reswrap(y) = res(y, x, p)
        return Roots.find_zero(reswrap, 1.0)
    end

    # solve for power coefficient
    p = ImplicitAD.implicit(solvewrap, res, [total_length; first_step_size; N])

    # return steps
    return bl_step_fun(1:N, first_step_size, p)
end
