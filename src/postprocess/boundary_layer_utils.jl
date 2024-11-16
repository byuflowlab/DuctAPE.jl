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
            duct_panel_lengths[1] # partial panel length
            [
                0.5 * (duct_panel_lengths[i] + duct_panel_lengths[i - 1]) for
                i in 3:length(duct_panel_lengths)
            ]
        ],
    )
end

"""
   split_at_stagnation_point(duct_panel_lengths, n_panels_casing)

Split the duct body surface at the leading edge of the duct.

# Arguments:
- `duct_panel_lengths::Vector{Float}` : Vector of panel lengths for the duct from casing trailing edge clockwise to nacelle trailing edge.
- `n_panels_casing::Int` : number of panels comprising the casing side of the duct

# Returns:
- `s_upper::Vector{Float}` : cumulative sum of upper side (nacelle) panel lengths
- `s_lower::Vector{Float}` : cumulative sum of lower side (casing) panel lengths
"""
function split_at_stagnation_point(duct_panel_lengths, duct_panel_tangents, Vtot_duct)

    # initialize stagnation point finder
    stag_ids = [1, 1]
    dp = zeros(eltype(Vtot_duct), 2)
    dp[1] = dot(duct_panel_tangents[:, 1], Vtot_duct[:, 1])

    # loop through panels and stop when dot product of panel vector and velocity vector changes sign
    for i in 2:length(duct_panel_lengths)
        dp[2] = dot(duct_panel_tangents[:, i], Vtot_duct[:, i])

        if sign(dp[1]) != sign(dp[2])
            stag_ids[:] .= [i - 1; i]
            break
        end
        dp[1] = dp[2]
    end

    if dp[1] == dp[2]
        # we're likely in a hover-ish case and there is no stagnation point.

        s_upper = nothing
        s_lower = arc_lengths_from_panel_lengths(duct_panel_lengths[end:-1:1])
        split_ratio = 1.0
        stag_ids .= length(duct_panel_lengths)
    else

        # interpolate the lengths between control points
        sum_length = 0.5 * sum(duct_panel_lengths[stag_ids])
        stag_interp = FLOWMath.linear(dp, [0.0, sum_length], 0.0)

        partial_panel_lengths = [stag_interp, sum_length - stag_interp]

        s_upper = arc_lengths_from_panel_lengths(
            [abs(partial_panel_lengths[2]); duct_panel_lengths[stag_ids[2]:end]]
        )

        s_lower = arc_lengths_from_panel_lengths(
            [abs(partial_panel_lengths[1]); duct_panel_lengths[stag_ids[1]:-1:1]]
        )
        split_ratio = stag_interp / sum_length
    end

    return s_upper, s_lower, stag_ids, split_ratio
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

function set_boundary_layer_steps(N::Int, total_length)
    return range(0.0, total_length; length=N)
end
