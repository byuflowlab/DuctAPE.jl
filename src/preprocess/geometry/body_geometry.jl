"""
    reinterpolate_bodies!(
        rp_duct_coordinates,
        rp_center_body_coordinates,
        duct_coordinates,
        center_body_coordinates,
        zwake,
        duct_le_coordinates,
        ncenter_body_inlet,
        nduct_inlet,
        finterp=FLOWMath.akima,
    )

Reinterpolate duct and center_body coordinates in order to make them compatible with the calculated wake sheet panel axial positions.

# Arguments
- `rp_duct_coordinates::Matrix{Float}` : the re-paneled duct coordinates
- `rp_center_body_coordinates::Matrix{Float}` : the re-paneled center_body coordinates
- `duct_coordinates::Matrix{Float}` : the input duct coordinates
- `center_body_coordinates::Matrix{Float}` : the input center_body coordinates
- `zwake::Vector{Float}` : the wake sheet panel node axial positions
- `duct_le_coordinates::Matrix{Float}` : [z r] coordinates of duct leading edge
- `ncenter_body_inlet::Int` : the number of panels to use for the center_body inlet
- `nduct_inlet,::Int` : the number of panels to use for the duct inlet

# Keyword Arguments
- `finterp::Function=FLOWMath.akima` : interpolation method
"""
function reinterpolate_bodies!(
    rp_duct_coordinates,
    rp_center_body_coordinates,
    duct_coordinates,
    center_body_coordinates,
    zwake,
    duct_le_coordinates,
    ncenter_body_inlet,
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

    # - interpolate center_body geometry to provided grid locations - #

    # rename for convenience
    center_body_z = view(center_body_coordinates, :, 1)
    center_body_r = view(center_body_coordinates, :, 2)

    center_body_te_id = min(length(zwake), searchsortedfirst(zwake, center_body_z[end]))
    center_body_in_wake_z = @view(zwake[1:center_body_te_id])
    center_body_in_wake_r = finterp(center_body_z, center_body_r, center_body_in_wake_z)

    center_body_inlet_length = center_body_in_wake_z[1] - center_body_z[1]
    center_body_inlet_z = scaled_cosine_spacing(
        ncenter_body_inlet + 1, 2 * center_body_inlet_length, center_body_z[1]; mypi=pi / 2.0
    )
    center_body_inlet_r = finterp(center_body_z, center_body_r, center_body_inlet_z)

    # assemble new duct coordinates
    rp_duct_coordinates .= hcat(
        [reverse(casing_in_wake_z)'; reverse(casing_in_wake_r)'],
        [reverse(casing_inlet_z)[2:end]'; reverse(casing_inlet_r)[2:end]'],
        [nacelle_inlet_z[2:end]'; nacelle_inlet_r[2:end]'],
        [nacelle_in_wake_z[2:end]'; nacelle_in_wake_r[2:end]'],
    )

    # assemble new center_body coordinates
    rp_center_body_coordinates .= hcat(
        [center_body_inlet_z[1:(end - 1)]'; center_body_inlet_r[1:(end - 1)]'],
        [center_body_in_wake_z'; center_body_in_wake_r'],
    )

    # check that the splining didn't put any of the center body radial coordinates in the negative.
    for rpcb in eachcol(rp_center_body_coordinates)
        if rpcb[2] < 0.0 && rpcb[2] > -2.0 * eps()
            rpcb[2] = 0.0
        end
    end

    return rp_duct_coordinates, rp_center_body_coordinates
end

"""
    place_duct!(rp_duct_coordinates, Rtip, rotor_axial_position, tip_gap)

Transform the duct radial coordinates such that the leading rotor radius touches the duct wall.

Note that this function is called AFTER the repanling function is called, such that the rotor_axial_position locations should line up directly with the duct and center_body coordinates.

# Arguments
- `rp_duct_coordinates::Matrix{Float}` : the re-paneled duct coordinates
- `Rtip::Vector{Float}` : Tip radii for the rotor(s)
- `rotor_axial_position::Vector{Float}` : axial position(s) of the rotor(s)
- `tip_gap::Vector{Float}` : tip gap for the fore-most rotor (MUST BE ZERO for now)
"""
function place_duct!(rp_duct_coordinates, Rtip, rotor_axial_position, tip_gap)

    # get current radial position of duct wall at leading rotor location
    _, iduct = findmin(x -> abs(x - rotor_axial_position), rp_duct_coordinates[1, :])
    rduct = rp_duct_coordinates[2, iduct]

    # - transform duct r-coordinates up by Rtip+tip gap of first rotor - #
    # need to account for current radial position of duct at leading rotor location as well.
    rp_duct_coordinates[2, :] .+= Rtip .+ tip_gap .- rduct

    return rp_duct_coordinates
end
