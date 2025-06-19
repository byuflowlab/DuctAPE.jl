#---------------------------------#
#         Setup functions         #
#---------------------------------#

"""
    arc_lengths_from_panel_lengths(duct_panel_lengths)

Cumulative sum of panel lengths for the given section of surface associated with the upper or lower boundary layer.

# Arguments
- `duct_panel_lengths::Vector{Float}` : vector of panel lengths (called influence_length in body_vortex_panels) associated with the duct (casing + nacelle).

# Returns
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
    split_at_stagnation_point(duct_panel_lengths, duct_panel_tangents, Vtot_duct, Vtan_duct, first_step_size)

Split the duct body surface at the leading edge of the duct by locating the stagnation point.

# Arguments:
- `duct_panel_lengths::Vector{Float}` : Vector of panel lengths for the duct from casing trailing edge clockwise to nacelle trailing edge.
- `duct_panel_tangents::Matrix{Float}` : Tangent vectors of each duct panel (dimension 2 × number_of_panels).
- `Vtot_duct::Matrix{Float}` : Total velocity vectors at each duct panel (dimension 2 × number_of_panels).
- `Vtan_duct::Vector{Float}` : Tangential velocity magnitude at each duct panel.
- `first_step_size::Float` : Reference step size used for detecting stagnation location.

# Returns:
- `s_upper::Union{Vector{Float}, Nothing}` : cumulative panel lengths of the upper (nacelle) side starting at stagnation point, or `nothing` if no stagnation point found.
- `s_lower::Vector{Float}` : cumulative panel lengths of the lower (casing) side starting at stagnation point.
- `stag_ids::Vector{Int}` : indices bounding the stagnation point on the duct panel cumulative length vector.
- `stag_point::Float` : arc length coordinate of the stagnation point along the duct surface.
- `split_ratio::Float` : ratio of stagnation point location relative to total duct length.
- `dots::Vector{Float}` : dot product values of panel tangents with total velocity vectors (used for stagnation detection).
"""
function split_at_stagnation_point(
    duct_panel_lengths, duct_panel_tangents, Vtot_duct, Vtan_duct, first_step_size
)
    offset = findfirst(x -> x > 1.1 * first_step_size, cumsum(duct_panel_lengths[end:-1:1]))

    s_tot = arc_lengths_from_panel_lengths(duct_panel_lengths[1:(end - offset)])

    dots = [
        dot(duct_panel_tangents[:, i], Vtot_duct[:, i]) for
        i in 1:(length(duct_panel_lengths) - offset)
    ]

    if allequal(sign.(dots))
        #=
            We're likely in a hover-ish case and there is no stagnation point.
            Want to find the minimum Vtan in this case, since either the stagnation point is the trailing edge, or it has lifted off the surface and the most recent point of stagnation could be used.
        =#

        _, minvtid = findmin(Vtan_duct)
        vtsp = smooth_Akima(s_tot, Vtan_duct)

        if minvtid >= length(s_tot)
            s_upper = nothing
            s_lower = arc_lengths_from_panel_lengths(duct_panel_lengths[end:-1:1])
            split_ratio = 1.0
            stag_point = s_tot[1]
            stag_ids = ones(Int, 2) .* length(duct_panel_lengths)

            return s_upper, s_lower, stag_ids, stag_point, split_ratio, dots
        else

            # print("in past bug area that is hard to recreate.  ")
            # print("min Vt index: ", minvtid)
            # println("  length s_tot: ", length(s_tot))
            # println("  length Vtan_duct: ", length(Vtan_duct))

            # - Check Bracket actually brackets - #
            bracket = (s_tot[max(minvtid - 1, 1)], s_tot[min(minvtid + 1, length(s_tot))])
            bracketvals = FLOWMath.derivative.(Ref(vtsp), bracket)

            # if not, get a bracket
            bidl = 0
            bidr = 1
            bid = 1
            while sign(bracketvals[1]) == sign(bracketvals[2])
                # not a bracketing interval
                lb = max(minvtid - bidl, 1)
                ub = min(minvtid + bidr, length(s_tot))
                bracket = (s_tot[lb], s_tot[ub])

                bracketvals = FLOWMath.derivative.(Ref(vtsp), bracket)

                if lb == 1 && ub == length(s_tot)
                    break
                else
                    if bid % 2 == 0
                        bidr += 2
                        bid += 1
                    else
                        bidt = bidl
                        bidl = bidr
                        bidr = bidt
                        bid += 1
                    end
                end
            end

            # printdebug("bracket:", bracket)
            # printdebug("bracketvals:", FLOWMath.derivative.(Ref(vtsp), bracket))

            # if still no bracket found, then try one last attempt without bracketing method
            if sign(bracketvals[1]) == sign(bracketvals[2])

                # println("using non-brent search")

                stag_point = Roots.find_zero(
                    x -> FLOWMath.derivative(vtsp, x), stot[minvtid]; atol=eps()
                )
                # use bracketing method if you can
            else

                # println("using brent search")

                stag_point = Roots.find_zero(
                    x -> FLOWMath.derivative(vtsp, x), bracket, Roots.Brent(); atol=eps()
                )
            end
        end

        # elseif mindotid == length(duct_panel_lengths) - offset

        #     s_upper = nothing
        #     s_lower = arc_lengths_from_panel_lengths(duct_panel_lengths[(end - 1):-1:1])
        #     split_ratio = 1.0
        #     stag_ids = ones(Int,2) .* (length(duct_panel_lengths) - offset)
    else
        mindotid = findfirst(x -> sign(x) != sign(dots[1]), dots)
        dotsp = smooth_Akima(s_tot, dots)

        # println("mindotid: ", mindotid)
        # println("mindot: ", dots[mindotid])
        # println("bracket: ", (s_tot[max(mindotid - 1, 1)], s_tot[mindotid]))
        # println(
        #     "bracket vals: ", dotsp([(s_tot[max(mindotid - 1, 1)], s_tot[mindotid])...])
        # )
        # println("dots: ", dots)

        stag_point = Roots.find_zero(
            dotsp, (s_tot[max(mindotid - 1, 1)], s_tot[mindotid]), Roots.Brent(); atol=eps()
        )
    end

    stag_ids = searchsortedfirst(s_tot, stag_point) .+ [-1, 0]
    split_ratio = stag_point / s_tot[end]

    partial_panel_lengths = [
        stag_point - s_tot[stag_ids[1]]
        s_tot[stag_ids[2]] - stag_point
    ]

    s_upper = arc_lengths_from_panel_lengths(
        [abs(partial_panel_lengths[2]); duct_panel_lengths[stag_ids[2]:end]]
    )

    s_lower = arc_lengths_from_panel_lengths(
        [abs(partial_panel_lengths[1]); duct_panel_lengths[stag_ids[1]:-1:1]]
    )

    return s_upper, s_lower, stag_ids, stag_point, split_ratio, dots
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

# Arguments
- `N::Int` : Number of steps to take
- `first_step_size::Float` : size of first step (which is `m` in `bl_step_fun`)
- `total_length::Float` : total surface length to divide up.

# Returns
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
