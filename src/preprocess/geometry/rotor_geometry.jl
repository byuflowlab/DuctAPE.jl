"""
    interpolate_blade_elements(
        rsp, Rtips, Rhubs, rotor_panel_centers, nbe; finterp=FLOWMath.linear
    )

Interpolate blade elements based on Rotor inputs and number of desired blade elements (from number of wake sheet in PanelingConstants input)

# Arguments
- `rsp::Rotor` : A Rotor object
- `Rtips::Vector{Float}' : Vector of rotor tip radii
- `Rhubs::Vector{Float}' : Vector of rotor hub radii
- `rotor_panel_centers::Vector{Float}' : Vector of rotor panel centers
- `nbe::Int` : number of blade elements per rotor

# Keyword Arguments
- `finterp::Function=FLOWMath.linear` : interpolation method (note, using Akima splines as is done for the body geometry can lead to negative chord in some cases)

# Returns
- `blade_element_cache::NamedTuple` : A named tuple containing the cacheable blade element information excluding the airfoil data.
- `airfoils::NamedTuple` : A named tuple containing vectors of inner and outer airfoil polar data for each blade element, used in interpolating the input data at blade element locations.
"""
function interpolate_blade_elements(
    rsp, Rtips, Rhubs, rotor_panel_centers, nbe; finterp=FLOWMath.linear
)
    nrotor = length(rsp.B)
    Rtip = Rtips
    Rhub = Rhubs
    B = rsp.B
    is_stator = rsp.is_stator
    chords = similar(rsp.chords, nbe, nrotor) .= 0
    twists = similar(rsp.twists, nbe, nrotor) .= 0
    stagger = similar(rsp.twists, nbe, nrotor) .= 0
    solidity = similar(rsp.chords, nbe, nrotor) .= 0
    # outer_airfoil = similar(rsp.airfoils[1], nbe, nrotor)
    # inner_airfoil = similar(rsp.airfoils[1], nbe, nrotor)
    outer_airfoil = Array{eltype(rsp.airfoils[1][1])}(undef, nbe, nrotor)
    inner_airfoil = Array{eltype(rsp.airfoils[1][1])}(undef, nbe, nrotor)
    inner_fraction = similar(rsp.r, nbe, nrotor) .= 0

    for irotor in 1:nrotor
        rpcs = rotor_panel_centers[(nbe * (irotor - 1) + 1):(nbe * irotor)]

        # dimensionalize the blade element radial positions
        rblade = finterp([0.0; 1.0], [0.0; Rtip[irotor]], rsp.r[:, irotor])

        # update chord lengths
        chords[:, irotor] .= finterp(rblade, rsp.chords[:, irotor], rpcs)

        # update twists
        twists[:, irotor] .= finterp(rblade, rsp.twists[:, irotor], rpcs)

        # update stagger
        stagger[:, irotor] .= get_stagger(twists[:, irotor])

        # update solidity
        solidity[:, irotor] .= get_local_solidity(B[irotor], chords[:, irotor], rpcs)

        for ir in 1:nbe
            # outer airfoil
            io = min(length(rblade), searchsortedfirst(rblade, rpcs[ir]))
            outer_airfoil[ir, irotor] = rsp.airfoils[irotor][io]

            # inner airfoil
            ii = max(1, io - 1)
            inner_airfoil[ir, irotor] = rsp.airfoils[irotor][ii]

            # fraction of inner airfoil's polars to use
            if rblade[io] == rblade[ii]
                inner_fraction[ir, irotor] = 1.0
            else
                inner_fraction[ir, irotor] =
                    (rpcs[ir] - rblade[ii]) / (rblade[io] - rblade[ii])
            end

            # Check incorrect extrapolation
            if inner_fraction[ir, irotor] > 1.0
                inner_fraction[ir, irotor] = 1.0
            end
        end
    end

    return (;
        Rtip,
        Rhub,
        rotor_panel_centers=reshape(rotor_panel_centers, (nbe, nrotor)),
        B,
        is_stator,
        chords,
        twists,
        stagger,
        solidity,
        inner_fraction,
        outer_airfoil,
        inner_airfoil,
    )
end

"""
    populate_airfoil_cache(cache, airfoil)

Copy the contents of airfoil into cache of the same type and size.
"""
function populate_airfoil_cache!(cache, airfoil)
    for pn in propertynames(airfoil)
        getproperty(cache, pn) .= getproperty(airfoil, pn)
    end
    return nothing
end

"""
    interpolate_blade_elements!(
        blade_element_cache, rsp, rotor_panel_centers, nbe; finterp=FLOWMath.linear
    )

In-place version of `interpolate_blade_elements`.

# Returns
- `airfoils::NamedTuple` : A named tuple containing vectors of inner and outer airfoil polar data for each blade element, used in interpolating the input data at blade element locations.
"""
function interpolate_blade_elements!(
    blade_element_cache, rsp, rotor_panel_centers, nbe; finterp=FLOWMath.linear
)
    nrotor = length(rsp.B)
    Rtip = blade_element_cache.Rtip .= rsp.Rtip
    Rhub = blade_element_cache.Rhub .= rsp.Rhub
    blade_element_cache.B .= rsp.B
    blade_element_cache.is_stator .= rsp.is_stator

    for irotor in 1:nrotor
        rpcs = @view(rotor_panel_centers[(nbe * (irotor - 1) + 1):(nbe * irotor)])

        # dimensionalize the blade element radial positions
        rblade = linear_transform((0.0, 1.0), (0.0, Rtip[irotor]), @view(rsp.r[:, irotor]))

        # update chord lengths
        blade_element_cache.chords[:, irotor] .= finterp(
            rblade, @view(rsp.chords[:, irotor]), rpcs
        )

        # update twists
        blade_element_cache.twists[:, irotor] .= finterp(
            rblade, @view(rsp.twists[:, irotor]), rpcs
        )

        # update stagger
        blade_element_cache.stagger[:, irotor] .= get_stagger(
            @view(blade_element_cache.twists[:, irotor])
        )

        # update solidity
        blade_element_cache.solidity[:, irotor] .= get_local_solidity(
            blade_element_cache.B[irotor],
            @view(blade_element_cache.chords[:, irotor]),
            rpcs,
        )

        for ir in 1:nbe
            # outer airfoil
            io = min(length(rblade), searchsortedfirst(rblade, rpcs[ir]))
            populate_airfoil_cache!(
                blade_element_cache.outer_airfoil[ir, irotor], rsp.airfoils[irotor][io]
            )

            # inner airfoil
            ii = max(1, io - 1)
            populate_airfoil_cache!(
                blade_element_cache.inner_airfoil[ir, irotor], rsp.airfoils[irotor][ii]
            )

            # fraction of inner airfoil's polars to use
            if rblade[io] == rblade[ii]
                blade_element_cache.inner_fraction[ir, irotor] = 1.0
            else
                blade_element_cache.inner_fraction[ir, irotor] =
                    (rpcs[ir] - rblade[ii]) / (rblade[io] - rblade[ii])
            end

            # Check incorrect extrapolation
            if blade_element_cache.inner_fraction[ir, irotor] > 1.0
                blade_element_cache.inner_fraction[ir, irotor] = 1.0
            end
        end
    end

    blade_element_cache.rotor_panel_centers .= reshape(rotor_panel_centers, (nbe, nrotor))

    return blade_element_cache
end

"""
    get_blade_ends_from_body_geometry(
        rp_duct_coordinates, rp_center_body_coordinates, tip_gaps, rotor_axial_position
    )

Obtain rotor hub and tip radii based on duct and center_body geometry.

# Arguments
- `var::type` :
- `rp_duct_coordinates::Matrix{Float}` : re-paneled duct coordinates
- `rp_center_body_coordinates::Matrix{Float}` : re-paneled center_body coordinates
- `tip_gaps::Vector{Float}` : gaps between blade tips and duct surface (MUST BE ZEROS for now)
- `rotor_axial_position::Vector{Float}` : rotor lifting line axial positions.

# Returns
- `Rtips::Vector{Float}` : rotor tip radii
- `Rhubs::Vector{Float}` : rotor hub radii
"""
function get_blade_ends_from_body_geometry(
    rp_duct_coordinates, rp_center_body_coordinates, tip_gaps, rotor_axial_position
)
    TF = promote_type(
        eltype(rp_duct_coordinates),
        eltype(rp_center_body_coordinates),
        eltype(tip_gaps),
        eltype(rotor_axial_position),
    )

    Rtip = zeros(TF, length(rotor_axial_position))
    Rhub = zeros(TF, length(rotor_axial_position))

    return get_blade_ends_from_body_geometry!(
        Rtip, Rhub, rp_duct_coordinates, rp_center_body_coordinates, tip_gaps, rotor_axial_position
    )
end

"""
    get_blade_ends_from_body_geometry!(
        Rtip,
        Rhub,
        rp_duct_coordinates,
        rp_center_body_coordinates,
        tip_gaps,
        rotor_axial_position;
        silence_warnings=true,
    )

In-place version of `get_blade_ends_from_body_geometry`.
"""
function get_blade_ends_from_body_geometry!(
    Rtip,
    Rhub,
    rp_duct_coordinates,
    rp_center_body_coordinates,
    tip_gaps,
    rotor_axial_position;
    silence_warnings=true,
)

    # - Get hub and tip wall indices - #
    ihub = zeros(Int, length(rotor_axial_position))
    iduct = zeros(Int, length(rotor_axial_position))
    for i in eachindex(rotor_axial_position)
        #indices
        _, ihub[i] = findmin(
            x -> abs(x - rotor_axial_position[i]), view(rp_center_body_coordinates, 1, :)
        )
        _, iduct[i] = findmin(
            x -> abs(x - rotor_axial_position[i]),
            view(rp_duct_coordinates, 1, 1:ceil(Int, size(rp_duct_coordinates, 2) / 2)),
        )
    end

    # - Add warnings about over writing Rhub and Rtip
    if !silence_warnings
        for (irotor, (R, r)) in enumerate(zip(Rtip, Rhub))
            if r !== rp_center_body_coordinates[2, ihub[irotor]]
                @info "Overwriting Rhub for rotor $(ForwardDiff.value(irotor)) to place it at the center_body wall.  Moving from $(ForwardDiff.value(r)) to $(ForwardDiff.value(rp_center_body_coordinates[2, ihub[irotor]]))"
            end
            if R !== rp_duct_coordinates[2, iduct[irotor]] .- tip_gaps[irotor]
                @info "Overwriting Rtip for rotor $(ForwardDiff.value(irotor)) to place it at the correct tip gap relative to the casing wall. Moving from $(ForwardDiff.value(R)) to $(ForwardDiff.value(rp_duct_coordinates[2, iduct[irotor]] .- tip_gaps[irotor]))"
            end
        end
    end

    # - Get hub and tip radial positions - #
    Rhub .= rp_center_body_coordinates[2, ihub]

    #need to shift the tips down by the distance of the tip gaps to get the actual tip radii
    #note that for stators, the tip gap should be zero anyway.
    Rtip .= rp_duct_coordinates[2, iduct] .- tip_gaps

    # check that the rotor radii aren't messed up
    for (irotor, (R, r)) in enumerate(zip(Rtip, Rhub))
        @assert R > r "Rotor #$(irotor) Tip Radius is set to be less than its Hub Radius. Consider setting the `autoshiftduct` option to true and/or check the input geometry."
    end

    return Rtip, Rhub
end

"""
    get_local_solidity(B, chord, r)

Calculate local solidity from local chord, radial position, and number of blades.

# Arguments
- `B::Float` : number of blades on rotor (usually an integer, but not necessarily).
- `chord::Vector{Float}` : chord lengths at each radial station.
- `r::Vector{Float}` : dimensional radial positions.

# Returns
- `solidity::Vector{Float}` : local solidity at each radial station
"""
function get_local_solidity(B, chord, r)
    return B .* chord ./ (2.0 * pi * r)
end

"""
    get_stagger(twists)

Convert twist angle to stagger angle
"""
function get_stagger(twists)
    return 0.5 * pi .- twists
end

"""
    generate_rotor_panels(rotor_axial_position, wake_grid, rotor_indices_in_wake, num_wake_sheets)

Generate rotor panel objects.

# Arguments
- `rotor_axial_position::Vector{Float}` : rotor lifting line axial position
- `wake_grid::Array{Float,3}` : wake elliptic grid axial and radial locations
- `rotor_indices_in_wake::Vector{Int}` : indices of where along wake the rotors are placed
- `num_wake_sheets::Int` : number of wake sheets

# Returns
- `rotor_source_panels::NamedTuple` : A named tuple containing the rotor source panel variables.
"""
function generate_rotor_panels(rotor_axial_position, wake_grid, rotor_indices_in_wake, num_wake_sheets)
    TF = promote_type(eltype(rotor_axial_position), eltype(wake_grid))

    xr = [zeros(TF, num_wake_sheets, 2) for i in 1:length(rotor_axial_position)]

    for irotor in eachindex(rotor_axial_position)
        @views xr[irotor] = [
            fill(rotor_axial_position[irotor], num_wake_sheets)'
            wake_grid[2, rotor_indices_in_wake[irotor], 1:num_wake_sheets]'
        ]
    end

    return generate_panels(xr; isbody=false, isrotor=true)
end

"""
    generate_rotor_panels!(
        rotor_source_panels, rotor_axial_position, wake_grid, rotor_indices_in_wake, num_wake_sheets
    )

In-place version of `generate_rotor_panels`.
"""
function generate_rotor_panels!(
    rotor_source_panels, rotor_axial_position, wake_grid, rotor_indices_in_wake, num_wake_sheets
)
    TF = promote_type(eltype(rotor_axial_position), eltype(wake_grid))

    xr = [zeros(TF, num_wake_sheets, 2) for i in 1:length(rotor_axial_position)]

    for irotor in eachindex(rotor_axial_position)
        @views xr[irotor] = [
            fill(rotor_axial_position[irotor], num_wake_sheets)'
            wake_grid[2, Int(rotor_indices_in_wake[irotor]), 1:num_wake_sheets]'
        ]
    end

    return generate_panels!(rotor_source_panels, xr; isbody=false, isrotor=true)
end
