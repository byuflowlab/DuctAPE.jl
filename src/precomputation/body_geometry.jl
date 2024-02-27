#"""
#"""
#function reinterpolate_bodies(
#    duct_coordinates,
#    centerbody_coordinates,
#    zwake,
#    ncenterbody_inlet,
#    nduct_inlet;
#    finterp=FLOWMath.akima,
#)

##TODO: initialize repaneled coordinates

#return  reinterpolate_bodies!(
#    duct_coordinates,
#    centerbody_coordinates,
#    zwake,
#    ncenterbody_inlet,
#    nduct_inlet;
#    finterp=FLOWMath.akima,
#)

#end

"""
"""
function reinterpolate_bodies!(
    rp_duct_coordinates,
    rp_centerbody_coordinates,
    duct_coordinates,
    centerbody_coordinates,
    zwake,
    ncenterbody_inlet,
    nduct_inlet;
    finterp=FLOWMath.akima,
)

    # - separate inner and outer duct coordinates - #
    dle, dleidx = findmin(view(duct_coordinates, :, 1))
    duct_inner_coordinates = view(duct_coordinates, dleidx:-1:1, :)
    duct_outer_coordinates = view(duct_coordinates, dleidx:size(duct_coordinates, 1), :)

    # - interpolate inner duct geometry to provided grid locations - #
    zduct_inner = view(duct_inner_coordinates, :, 1)
    rduct_inner = view(duct_inner_coordinates, :, 2)
    duct_trailing_index = min(length(zwake), searchsortedfirst(zwake, zduct_inner[end]))
    new_zduct_grid = zwake[1:duct_trailing_index]
    new_rduct_grid = finterp(zduct_inner, rduct_inner, new_zduct_grid)

    # - interpolate outer duct geometry to provided grid locations - #
    zduct_outer = view(duct_outer_coordinates, :, 1)
    rduct_outer = view(duct_outer_coordinates, :, 2)

    new_rduct_grid_outer = finterp(zduct_outer, rduct_outer, new_zduct_grid)

    # - interpolate hub geometry to provided grid locations - #
    zhub = view(centerbody_coordinates, :, 1)
    rhub = view(centerbody_coordinates, :, 2)
    centerbody_trailing_index = min(length(zwake), searchsortedfirst(zwake, zhub[end]))
    new_zcenterbody_grid = zwake[1:centerbody_trailing_index]
    new_rcenterbody_grid = finterp(zhub, rhub, new_zcenterbody_grid)

    # - update hub geometry between leading edge and rotor - #
    if isnothing(nduct_inlet)
        # find the index of the first rotor on the hub
        ridx = searchsortedfirst(zhub, zwake[1])

        # use existing hub geometry between leading edge and rotor
        new_zhub = view(zhub, 1:ridx, 1)
        new_rhub = view(rhub, 1:ridx, 2)
    else
        # interpolate hub geometry between leading edge and rotor
        # new_zhub = range(zhub[1, 1], new_zcenterbody_grid[1], ncenterbody_inlet + 1)

        scale = new_zcenterbody_grid[1] - zhub[1]
        transform = zhub[1]
        new_zhub = scaled_cosine_spacing(
            ncenterbody_inlet + 1, 2 * scale, transform; mypi=pi / 2
        )
        new_rhub = finterp(zhub, rhub, new_zhub)
    end

    # - update inner duct geometry between leading edge and rotor - #
    if isnothing(nduct_inlet)
        # find the index of the first rotor on the duct inside surface
        ridx = length(zduct_inner) - searchsortedfirst(zduct_inner, zwake[1]) + 1

        # use existing inner duct geometry between leading edge and rotor
        new_zduct_inner = view(duct_coordinates, dleidx:-1:ridx, 1)
        new_rduct_inner = view(duct_coordinates, dleidx:-1:ridx, 2)
    else
        # interpolate inner duct geometry between leading edge and rotor
        # new_zduct_inner = range(
        #     duct_coordinates[dleidx, 1], new_zduct_grid[1], nduct_inlet + 1
        # )

        scale = new_zduct_grid[1] - duct_coordinates[dleidx, 1]
        transform = duct_coordinates[dleidx, 1]
        new_zduct_inner = scaled_cosine_spacing(
            nduct_inlet + 1, 2 * scale, transform; mypi=pi / 2
        )
        new_rduct_inner = finterp(zduct_inner, rduct_inner, new_zduct_inner)
    end

    # update outer duct geometry
    if isnothing(nduct_inlet)
        new_zduct_outer = zduct_outer
        new_rduct_outer = rduct_outer
    else
        # new_zduct_outer = range(zduct_outer[1, 1], zduct_outer[end, 1], nduct_outer + 1)

        # in order to make sure that the trailing edge panels are (nearly) the same length it's easiest just to use the same x-spacing for the outer side of the duct as well.
        new_rduct_outer = finterp(zduct_outer, rduct_outer, new_zduct_inner)
    end

    # assemble new duct coordinates
    rp_duct_coordinates .= hcat(
        [reverse(new_zduct_grid)'; reverse(new_rduct_grid)'],
        [reverse(new_zduct_inner)[2:end]'; reverse(new_rduct_inner)[2:end]'],
        [new_zduct_inner[2:(end - 1)]'; new_rduct_outer[2:(end - 1)]'],
        [new_zduct_grid'; new_rduct_grid_outer'],
    )

    # assemble new hub coordinates
    rp_centerbody_coordinates .= hcat(
        [new_zhub[1:(end - 1)]'; new_rhub[1:(end - 1)]'],
        [new_zcenterbody_grid'; new_rcenterbody_grid'],
    )

    return rp_duct_coordinates, rp_centerbody_coordinates
end

"""
transforms duct radial coordinates such that the leading rotor radius touches the duct wall.
Also finds the various rotor hub and tip radii based on the hub and duct geometry
Note that this function is called AFTER the repanling function is called, such that the rotorzloc locations should line up directly with the duct and hub coordinates.
"""
function place_duct!(duct_coordinates, Rtip, rotorzloc, tip_gap)

    # get current radial position of duct wall at leading rotor location
    _, iduct = findmin(x -> abs(x - rotorzloc), view(duct_coordinates, 1, :))
    rduct = duct_coordinates[2, iduct]

    # - transform duct r-coordinates up by Rtip+tip gap of first rotor - #
    # need to account for current radial position of duct at leading rotor location as well.
    duct_coordinates[2, :] .+= Rtip .+ tip_gap .- rduct

    return duct_coordinates
end
