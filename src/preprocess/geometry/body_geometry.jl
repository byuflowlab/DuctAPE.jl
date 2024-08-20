"""
    reinterpolate_bodies!(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        duct_coordinates,
        centerbody_coordinates,
        zwake,
        ncenterbody_inlet,
        nduct_inlet;
        finterp=FLOWMath.akima,
    )

Reinterpolate duct and centerbody coordinates in order to make them compatible with the calculated wake sheet panel axial positions.

# Arguments
- `rp_duct_coordinates::Matrix{Float}` : the re-paneled duct coordinates
- `rp_centerbody_coordinates::Matrix{Float}` : the re-paneled centerbody coordinates
- `duct_coordinates::Matrix{Float}` : the input duct coordinates
- `centerbody_coordinates::Matrix{Float}` : the input centerbody coordinates
- `zwake::Matrix{Float}` : the wake sheet panel node axial positions
- `ncenterbody_inlet::Matrix{Float}` : the number of panels to use for the centerbody inlet
- `nduct_inlet::Matrix{Float}` : the number of panels to use for the duct inlet

# Keyword Arguments
- `finterp::Function=FLOWMath.akima` : interpolation method
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
    new_zduct_grid = @view(zwake[1:duct_trailing_index])
    new_rduct_grid = finterp(zduct_inner, rduct_inner, new_zduct_grid)

    # - interpolate outer duct geometry to provided grid locations - #
    zduct_outer = view(duct_outer_coordinates, :, 1)
    rduct_outer = view(duct_outer_coordinates, :, 2)

    new_rduct_grid_outer = finterp(zduct_outer, rduct_outer, new_zduct_grid)

    # - interpolate centerbody geometry to provided grid locations - #
    zcenterbody = view(centerbody_coordinates, :, 1)
    rcenterbody = view(centerbody_coordinates, :, 2)
    centerbody_trailing_index = min(
        length(zwake), searchsortedfirst(zwake, zcenterbody[end])
    )
    new_zcenterbody_grid = @view(zwake[1:centerbody_trailing_index])
    new_rcenterbody_grid = finterp(zcenterbody, rcenterbody, new_zcenterbody_grid)

    scale = new_zcenterbody_grid[1] - zcenterbody[1]
    transform = zcenterbody[1]
    new_zcenterbody = scaled_cosine_spacing(
        ncenterbody_inlet + 1, 2 * scale, transform; mypi=pi / 2
    )
    new_rcenterbody = finterp(zcenterbody, rcenterbody, new_zcenterbody)

    scale = new_zduct_grid[1] - duct_coordinates[dleidx, 1]
    transform = duct_coordinates[dleidx, 1]
    new_zduct_inner = scaled_cosine_spacing(
        nduct_inlet + 1, 2 * scale, transform; mypi=pi / 2
    )
    new_rduct_inner = finterp(zduct_inner, rduct_inner, new_zduct_inner)

    # update outer duct geometry
    # in order to make sure that the trailing edge panels are (nearly) the same length it's easiest just to use the same x-spacing for the outer side of the duct as well.
    new_rduct_outer = finterp(zduct_outer, rduct_outer, new_zduct_inner)

    # assemble new duct coordinates
    rp_duct_coordinates .= hcat(
        [reverse(new_zduct_grid)'; reverse(new_rduct_grid)'],
        [reverse(new_zduct_inner)[2:end]'; reverse(new_rduct_inner)[2:end]'],
        [new_zduct_inner[2:(end - 1)]'; new_rduct_outer[2:(end - 1)]'],
        [new_zduct_grid'; new_rduct_grid_outer'],
    )

    # assemble new centerbody coordinates
    rp_centerbody_coordinates .= hcat(
        [new_zcenterbody[1:(end - 1)]'; new_rcenterbody[1:(end - 1)]'],
        [new_zcenterbody_grid'; new_rcenterbody_grid'],
    )

    # check that the splining didn't put any of the center body radial coordinates in the negative.
    for rpcb in eachcol(rp_centerbody_coordinates)
        if rpcb[2] < 0.0 && rpcb[2] > -2.0*eps()
            rpcb[2] = 0.0
        end
    end

    return rp_duct_coordinates, rp_centerbody_coordinates
end

"""
    place_duct!(rp_duct_coordinates, Rtip, rotorzloc, tip_gap)

Transform the duct radial coordinates such that the leading rotor radius touches the duct wall.

Note that this function is called AFTER the repanling function is called, such that the rotorzloc locations should line up directly with the duct and centerbody coordinates.

# Arguments
- `rp_duct_coordinates::Matrix{Float}` : the re-paneled duct coordinates
- `Rtip::Vector{Float}` : Tip radii for the rotor(s)
- `rotorzloc::Vector{Float}` : axial position(s) of the rotor(s)
- `tip_gap::Vector{Float}` : tip gap for the fore-most rotor (MUST BE ZERO for now)
"""
function place_duct!(rp_duct_coordinates, Rtip, rotorzloc, tip_gap)

    # get current radial position of duct wall at leading rotor location
    _, iduct = findmin(x -> abs(x - rotorzloc), rp_duct_coordinates[1, :])
    rduct = rp_duct_coordinates[2, iduct]

    # - transform duct r-coordinates up by Rtip+tip gap of first rotor - #
    # need to account for current radial position of duct at leading rotor location as well.
    rp_duct_coordinates[2, :] .+= Rtip .+ tip_gap .- rduct

    return rp_duct_coordinates
end
