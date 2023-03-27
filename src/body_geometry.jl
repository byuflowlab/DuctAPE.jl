"""
"""
function update_body_geometry(
    duct_coordinates,
    hub_coordinates,
    xwake,
    nhub,
    nduct_inner,
    nduct_outer;
    finterp=FLOWMath.akima,
)

    # separate inner and outer duct coordinates
    _, leidx = findmin(view(duct_coordinates, :, 1))
    duct_inner_coordinates = view(duct_coordinates, leidx:-1:1, :)
    duct_outer_coordinates = view(duct_coordinates, leidx:size(duct_coordinates, 1), :)

    # existing outer duct geometry
    xduct_outer = view(duct_outer_coordinates, :, 1)
    rduct_outer = view(duct_outer_coordinates, :, 2)

    # update outer duct geometry
    if isnothing(nduct_inner)
        new_xduct_outer = xduct_outer
        new_rduct_outer = rduct_outer
    else
        new_xduct_outer = range(xduct_outer[1, 1], xduct_outer[end, 1], nduct_outer + 1)
        new_rduct_outer = finterp(xduct_outer, rduct_outer, new_xduct_outer)
    end

    # interpolate inner duct geometry to provided grid locations
    xduct_inner = view(duct_inner_coordinates, :, 1)
    rduct_inner = view(duct_inner_coordinates, :, 2)
    duct_trailing_index = min(length(xwake), searchsortedfirst(xwake, xduct_inner[end]))
    new_xduct_grid = xwake[1:duct_trailing_index]
    new_rduct_grid = finterp(xduct_inner, rduct_inner, new_xduct_grid)

    #TODO: something isn't right around here. see plots for non-symmetric geometry
    #probably ask taylor what he's trying to do here
    #also consider if this is even a good approach to the overall problem
    # find the index of the first rotor on the duct inside surface
    ridx = length(xduct_inner) - searchsortedfirst(xduct_inner, xwake[1]) + 1

    # update inner duct geometry between leading edge and rotor
    if isnothing(nduct_inner)
        # use existing inner duct geometry between leading edge and rotor
        new_xduct_inner = view(duct_coordinates, leidx:-1:ridx, 1)
        new_rduct_inner = view(duct_coordinates, leidx:-1:ridx, 2)
    else
        # interpolate inner duct geometry between leading edge and rotor
        new_xduct_inner = range(
            duct_coordinates[leidx, 1], duct_coordinates[ridx, 1], nduct_inner + 1
        )
        new_rduct_inner = finterp(xduct_inner, rduct_inner, new_xduct_inner)
    end

    # interpolate hub geometry to provided grid locations
    xhub = view(hub_coordinates, :, 1)
    rhub = view(hub_coordinates, :, 2)
    hub_trailing_index = min(length(xwake), searchsortedfirst(xwake, xhub[end]))
    new_xhub_grid = xwake[1:hub_trailing_index]
    new_rhub_grid = finterp(xhub, rhub, new_xhub_grid)

    # find the index of the first rotor on the hub
    ridx = searchsortedfirst(xhub, xwake[1])

    # update hub geometry between leading edge and rotor
    if isnothing(nduct_inner)
        # use existing hub geometry between leading edge and rotor
        new_xhub = view(xhub, 1:ridx, 1)
        new_rhub = view(rhub, 1:ridx, 2)
    else
        # interpolate hub geometry between leading edge and rotor
        new_xhub = range(xhub[1, 1], xhub[ridx, 1], nhub + 1)
        new_rhub = finterp(xhub, rhub, new_xhub)
    end

    # assemble new duct coordinates
    updated_duct_coordinates = vcat(
        hcat(reverse(new_xduct_grid), reverse(new_rduct_grid)),
        hcat(reverse(new_xduct_inner)[2:end], reverse(new_rduct_inner)[2:end]),
        hcat(new_xduct_outer[2:end], new_rduct_outer[2:end]),
    )

    # assemble new hub coordinates
    updated_hub_coordinates = vcat(
        hcat(new_xhub[1:(end - 1)], new_rhub[1:(end - 1)]),
        hcat(new_xhub_grid, new_rhub_grid),
    )

    return updated_duct_coordinates, updated_hub_coordinates
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
function generate_body_geometry(
    duct_coordinates,
    hub_coordinates;
    method=ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true]),
)
    body_panels = ff.generate_panels(method, (duct_coordinates, hub_coordinates))

    return body_panels
end
