"""
"""
function update_body_geometry(
    duct_coordinates, hub_coordinates, xwake, nhub_inlet, nduct_inlet; finterp=FLOWMath.akima
)

    # - separate inner and outer duct coordinates - #
    dle, dleidx = findmin(view(duct_coordinates, :, 1))
    duct_inner_coordinates = view(duct_coordinates, dleidx:-1:1, :)
    duct_outer_coordinates = view(duct_coordinates, dleidx:size(duct_coordinates, 1), :)

    # - interpolate inner duct geometry to provided grid locations - #
    xduct_inner = view(duct_inner_coordinates, :, 1)
    rduct_inner = view(duct_inner_coordinates, :, 2)
    duct_trailing_index = min(length(xwake), searchsortedfirst(xwake, xduct_inner[end]))
    new_xduct_grid = xwake[1:duct_trailing_index]
    new_rduct_grid = finterp(xduct_inner, rduct_inner, new_xduct_grid)

    # - interpolate outer duct geometry to provided grid locations - #
    xduct_outer = view(duct_outer_coordinates, :, 1)
    rduct_outer = view(duct_outer_coordinates, :, 2)

    new_rduct_grid_outer = finterp(xduct_outer, rduct_outer, new_xduct_grid)

    # - interpolate hub geometry to provided grid locations - #
    xhub = view(hub_coordinates, :, 1)
    rhub = view(hub_coordinates, :, 2)
    hub_trailing_index = min(length(xwake), searchsortedfirst(xwake, xhub[end]))
    new_xhub_grid = xwake[1:hub_trailing_index]
    new_rhub_grid = finterp(xhub, rhub, new_xhub_grid)

    # - update hub geometry between leading edge and rotor - #
    if isnothing(nduct_inlet)
        # find the index of the first rotor on the hub
        ridx = searchsortedfirst(xhub, xwake[1])

        # use existing hub geometry between leading edge and rotor
        new_xhub = view(xhub, 1:ridx, 1)
        new_rhub = view(rhub, 1:ridx, 2)
    else
        # interpolate hub geometry between leading edge and rotor
        # new_xhub = range(xhub[1, 1], new_xhub_grid[1], nhub_inlet + 1)

        scale = new_xhub_grid[1] - xhub[1]
        transform = xhub[1]
        new_xhub = scaled_cosine_spacing(nhub_inlet + 1, 2 * scale, transform; mypi=pi / 2)
        new_rhub = finterp(xhub, rhub, new_xhub)
    end

    # - update inner duct geometry between leading edge and rotor - #
    if isnothing(nduct_inlet)
        # find the index of the first rotor on the duct inside surface
        ridx = length(xduct_inner) - searchsortedfirst(xduct_inner, xwake[1]) + 1

        # use existing inner duct geometry between leading edge and rotor
        new_xduct_inner = view(duct_coordinates, dleidx:-1:ridx, 1)
        new_rduct_inner = view(duct_coordinates, dleidx:-1:ridx, 2)
    else
        # interpolate inner duct geometry between leading edge and rotor
        # new_xduct_inner = range(
        #     duct_coordinates[dleidx, 1], new_xduct_grid[1], nduct_inlet + 1
        # )

        scale = new_xduct_grid[1] - duct_coordinates[dleidx, 1]
        transform = duct_coordinates[dleidx, 1]
        new_xduct_inner = scaled_cosine_spacing(
            nduct_inlet + 1, 2 * scale, transform; mypi=pi / 2
        )

        new_rduct_inner = finterp(xduct_inner, rduct_inner, new_xduct_inner)
    end

    # update outer duct geometry
    if isnothing(nduct_inlet)
        new_xduct_outer = xduct_outer
        new_rduct_outer = rduct_outer
    else
        # new_xduct_outer = range(xduct_outer[1, 1], xduct_outer[end, 1], nduct_outer + 1)

        # in order to make sure that the trailing edge panels are (nearly) the same length it's easiest just to use the same x-spacing for the outer side of the duct as well.
        new_rduct_outer = finterp(xduct_outer, rduct_outer, new_xduct_inner)
    end

    # assemble new duct coordinates
    updated_duct_coordinates = vcat(
        hcat(reverse(new_xduct_grid), reverse(new_rduct_grid)),
        hcat(reverse(new_xduct_inner)[2:end], reverse(new_rduct_inner)[2:end]),
        hcat(new_xduct_inner[2:(end - 1)], new_rduct_outer[2:(end - 1)]),
        hcat(new_xduct_grid, new_rduct_grid_outer),
    )

    # assemble new hub coordinates
    updated_hub_coordinates = vcat(
        hcat(new_xhub[1:(end - 1)], new_rhub[1:(end - 1)]),
        hcat(new_xhub_grid, new_rhub_grid),
    )

    return updated_duct_coordinates, updated_hub_coordinates
end

"""
transforms duct radial coordinates such that the leading rotor radius touches the duct wall.
Also finds the various rotor hub and tip radii based on the hub and duct geometry
Note that this function is called AFTER the repanling function is called, such that the xrotor locations should line up directly with the duct and hub coordinates.
"""
function place_duct(duct_coordinates, hub_coordinates, Rtip, tip_gaps, xrotors)

    # - Get hub and tip wall indices - #
    ihub = zeros(Int, length(xrotors))
    iduct = zeros(Int, length(xrotors))
    for i in 1:length(xrotors)
        #indices
        _, ihub[i] = findmin(x -> abs(x - xrotors[i]), view(hub_coordinates, :, 1))
        _, iduct[i] = findmin(x -> abs(x - xrotors[i]), view(duct_coordinates, :, 1))
    end

    # get current radial position of duct wall at leading rotor location
    rduct = duct_coordinates[iduct[1], 2]

    # - transform duct r-coordinates up by Rtip+tip gap of first rotor - #
    # need to account for current radial position of duct at leading rotor location as well.
    duct_coordinates[:, 2] .+= Rtip .+ tip_gaps[1] .- rduct[1]

    # - Get hub and tip radial positions - #
    Rhubs = hub_coordinates[ihub, 2]
    #need to shift the tips down by the distance of the tip gaps to get the actual tip radii
    #note that for stators, the tip gap should be zero anyway.
    Rtips = duct_coordinates[iduct, 2] .- tip_gaps

    return duct_coordinates, Rtips, Rhubs
end

"""
    generate_body_panels(duct_coordinates, hub_coordinates, discretization; kwargs...)

Define the paneling (see [`FLOWFoil.AxisymmetricPanel`](@ref)) for the duct and hub.

# Arguments:
- `duct_coordinates::Array{Float64,2}` : duct coordinates, starting from trailing edge, going clockwise
- `hub_coordinates::Array{Float64,2}` : hub coordinates, starting from leading edge

# Keyword Arguments:
- `method` : Axisymmetric method as defined by FLOWFoil, defaults to `AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true])`

"""
function generate_body_panels(duct_coordinates, hub_coordinates)
    # use FLOWFoil to generate panels
    if duct_coordinates != nothing
        duct_method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false])
        duct_panels = ff.generate_panels(duct_method, duct_coordinates)
    end

    if hub_coordinates != nothing
        hub_method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true])
        hub_panels = ff.generate_panels(hub_method, hub_coordinates)
    end

    if duct_coordinates != nothing && hub_coordinates != nothing
        return [duct_panels; hub_panels]
    elseif duct_coordinates == nothing && hub_coordinates != nothing
        return [hub_panels]
    elseif duct_coordinates != nothing && hub_coordinates == nothing
        return [duct_panels]
    else
        return []
    end
end

"""
"""
function scaled_cosine_spacing(N, scale, transform; mypi=pi)
    return transform .+ scale * [0.5 * (1 - cos(mypi * (i - 1) / (N - 1))) for i in 1:N]
end

# - Function for adding in xlocations - #
# from https://stackoverflow.com/questions/25678112/insert-item-into-a-sorted-list-with-julia-with-and-without-duplicates
function insert_and_dedup!(v, x)
    for i in 1:length(x)
        # find ranges and replace with discrete values (thus deleting duplicates if present)
        v = (splice!(v, searchsorted(v, x[i]), x[i]); v)
    end
end
