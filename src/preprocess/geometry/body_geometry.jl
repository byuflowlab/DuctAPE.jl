"""
    reinterpolate_bodies!(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        duct_coordinates,
        centerbody_coordinates,
        zwake,
        duct_le_coordinates,
        ncenterbody_inlet,
        nduct_inlet,
        finterp=FLOWMath.akima,
    )

Reinterpolate duct and centerbody coordinates in order to make them compatible with the calculated wake sheet panel axial positions.

# Arguments
- `rp_duct_coordinates::Matrix{Float}` : the re-paneled duct coordinates
- `rp_centerbody_coordinates::Matrix{Float}` : the re-paneled centerbody coordinates
- `duct_coordinates::Matrix{Float}` : the input duct coordinates
- `centerbody_coordinates::Matrix{Float}` : the input centerbody coordinates
- `zwake::Vector{Float}` : the wake sheet panel node axial positions
- `duct_le_coordinates::Matrix{Float}` : [z r] coordinates of duct leading edge
- `ncenterbody_inlet::Int` : the number of panels to use for the centerbody inlet
- `nduct_inlet,::Int` : the number of panels to use for the duct inlet

# Keyword Arguments
- `finterp::Function=FLOWMath.akima` : interpolation method
"""
function reinterpolate_bodies!(
    rp_duct_coordinates,
    rp_centerbody_coordinates,
    duct_coordinates,
    centerbody_coordinates,
    zwake,
    duct_le_coordinates,
    ncenterbody_inlet,
    nduct_inlet;
    finterp=FLOWMath.akima,
)

    # - separate inner and outer duct coordinates - #
    dle, dleidx = findmin(view(duct_coordinates, :, 1))
    duct_inlet_length = zwake[1] - duct_le_coordinates[1]

    if duct_le_coordinates[2] > duct_coordinates[dleidx, 2]
        duct_inner_coordinates = [
            duct_le_coordinates
            view(duct_coordinates, dleidx:-1:1, :)
        ]
        duct_outer_coordinates = [
            duct_le_coordinates
            view(duct_coordinates, (dleidx + 1):size(duct_coordinates, 1), :)
        ]
    elseif duct_le_coordinates[2] < duct_coordinates[dleidx, 2]
        duct_inner_coordinates = [
            duct_le_coordinates
            view(duct_coordinates, (dleidx - 1):-1:1, :)
        ]
        duct_outer_coordinates = [
            duct_le_coordinates
            view(duct_coordinates, dleidx:size(duct_coordinates, 1), :)
        ]
    else
        # coordinate is the zero
        duct_inner_coordinates = view(duct_coordinates, dleidx:-1:1, :)
        duct_outer_coordinates = view(duct_coordinates, dleidx:size(duct_coordinates, 1), :)
    end

    # - interpolate inner duct geometry to provided grid locations - #

    # rename for convenience
    casing_z = view(duct_inner_coordinates, :, 1)
    casing_r = view(duct_inner_coordinates, :, 2)

    # index of casing trailing edge in wake discretization
    casing_te_id = min(length(zwake), searchsortedfirst(zwake, casing_z[end]))

    # new points for casing aft of rotor
    casing_in_wake_z = @view(zwake[1:casing_te_id])
    casing_in_wake_r = finterp(casing_z, casing_r, casing_in_wake_z)

    # repaneled casing inlet nodes
    casing_inlet_z = scaled_cosine_spacing(
        nduct_inlet + 1, 2 * duct_inlet_length, duct_le_coordinates[1]; mypi=pi / 2.0
    )
    casing_inlet_r = finterp(casing_z, casing_r, casing_inlet_z)

    # - interpolate outer duct geometry to provided grid locations - #

    # rename for convenience
    nacelle_z = view(duct_outer_coordinates, :, 1)
    nacelle_r = view(duct_outer_coordinates, :, 2)

    nacelle_in_wake_z = collect(
        range(zwake[1], nacelle_z[end]; length=length(casing_in_wake_z))
    )
    nacelle_in_wake_r = finterp(nacelle_z, nacelle_r, nacelle_in_wake_z)

    # repaneled nacelle inlet nodes
    nacelle_inlet_z = scaled_cosine_spacing(
        nduct_inlet + 1, 2 * duct_inlet_length, duct_le_coordinates[1]; mypi=pi / 2.0
    )
    nacelle_inlet_r = finterp(nacelle_z, nacelle_r, nacelle_inlet_z)

    # - interpolate centerbody geometry to provided grid locations - #

    # rename for convenience
    centerbody_z = view(centerbody_coordinates, :, 1)
    centerbody_r = view(centerbody_coordinates, :, 2)

    centerbody_te_id = min(length(zwake), searchsortedfirst(zwake, centerbody_z[end]))
    centerbody_in_wake_z = @view(zwake[1:centerbody_te_id])
    centerbody_in_wake_r = finterp(centerbody_z, centerbody_r, centerbody_in_wake_z)

    centerbody_inlet_length = centerbody_in_wake_z[1] - centerbody_z[1]
    centerbody_inlet_z = scaled_cosine_spacing(
        ncenterbody_inlet + 1, 2 * centerbody_inlet_length, centerbody_z[1]; mypi=pi / 2.0
    )
    centerbody_inlet_r = finterp(centerbody_z, centerbody_r, centerbody_inlet_z)

    # assemble new duct coordinates
    rp_duct_coordinates .= hcat(
        [reverse(casing_in_wake_z)'; reverse(casing_in_wake_r)'],
        [reverse(casing_inlet_z)[2:end]'; reverse(casing_inlet_r)[2:end]'],
        [nacelle_inlet_z[2:end]'; nacelle_inlet_r[2:end]'],
        [nacelle_in_wake_z[2:end]'; nacelle_in_wake_r[2:end]'],
    )

    # assemble new centerbody coordinates
    rp_centerbody_coordinates .= hcat(
        [centerbody_inlet_z[1:(end - 1)]'; centerbody_inlet_r[1:(end - 1)]'],
        [centerbody_in_wake_z'; centerbody_in_wake_r'],
    )

    # check that the splining didn't put any of the center body radial coordinates in the negative.
    for rpcb in eachcol(rp_centerbody_coordinates)
        if rpcb[2] < 0.0 && rpcb[2] > -2.0 * eps()
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
