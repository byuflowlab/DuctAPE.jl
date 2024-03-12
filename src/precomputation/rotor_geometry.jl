# """
#     BladeElements{TF, TAF}

# # Fields:
# - `B::Int`: number of blades on the rotor
# - `Omega::TF` : rotor rotation rate (rad/s)
# - `rotorzloc::TF`: x-position of the rotor (dimensionalized)
# - `r::Vector{TF}`: radial coordinate of each blade element (dimensionalized)
# - `chords::Vector{TF}` : chord length of each blade element
# - `twists::Vector{TF}` : twist of each blade element (in radians)
# - `solidity::Vector{TF}` : solidity of each blade element
# - `outer_airfoils::Vector{TAF}` : outer bounding airfoil polar
# - `inner_airfoils::Vector{TAF}` : inner bounding airfoil polar
# - `inner_fraction::Vector{TAF}`: fraction of inner bounding airfoil polar to use
# """
# struct BladeElements{TF,TAF}
#     B::Int
#     Omega::TF
#     rotorzloc::TF
#     r::Vector{TF}
#     chords::Vector{TF}
#     twists::Vector{TF}
#     solidity::Vector{TF}
#     outer_airfoils::Vector{TAF}
#     inner_airfoils::Vector{TAF}
#     inner_fraction::Vector{TF}
# end

"""
"""
function get_blade_ends_from_body_geometry(
    duct_coordinates, centerbody_coordinates, tip_gaps, rotorzloc
)
    TF = promote_type(
        eltype(duct_coordinates),
        eltype(centerbody_coordinates),
        eltype(tip_gaps),
        eltype(rotorzloc),
    )

    Rtip = zeros(TF, length(rotorzloc))
    Rhub = zeros(TF, length(rotorzloc))

    return get_blade_ends_from_body_geometry!(
        Rtip, Rhub, duct_coordinates, centerbody_coordinates, tip_gaps, rotorzloc
    )
end

function get_blade_ends_from_body_geometry!(
    Rtip,
    Rhub,
    duct_coordinates,
    centerbody_coordinates,
    tip_gaps,
    rotorzloc;
    silence_warnings=true,
)

    # - Get hub and tip wall indices - #
    ihub = zeros(Int, length(rotorzloc))
    iduct = zeros(Int, length(rotorzloc))
    for i in 1:length(rotorzloc)
        #indices
        _, ihub[i] = findmin(x -> abs(x - rotorzloc[i]), view(centerbody_coordinates, 1, :))
        _, iduct[i] = findmin(x -> abs(x - rotorzloc[i]), view(duct_coordinates, 1, :))
    end

    # - Add warnings about over writing Rhub and Rtip
    if !silence_warnings
        for (irotor, (R, r)) in enumerate(zip(Rtip, Rhub))
            if r !== centerbody_coordinates[2, ihub[irotor]]
                @warn "Overwriting Rhub for rotor $(irotor) to place it at the centerbody wall.  Moving from $(r) to $(centerbody_coordinates[2, ihub[irotor]])"
            end
            if R !== duct_coordinates[2, iduct[irotor]] .- tip_gaps[irotor]
                @warn "Overwriting Rtip for rotor $(irotor) to place it at the correct tip gap relative to the casing wall. Moving from $(R) to $(duct_coordinates[2, iduct[irotor]] .- tip_gaps[irotor])"
            end
        end
    end

    # - Get hub and tip radial positions - #
    Rhub .= centerbody_coordinates[2, ihub]

    #need to shift the tips down by the distance of the tip gaps to get the actual tip radii
    #note that for stators, the tip gap should be zero anyway.
    Rtip .= duct_coordinates[2, iduct] .- tip_gaps

    # check that the rotor radii aren't messed up
    for (irotor, (R, r)) in enumerate(zip(Rtip, Rhub))
        @assert R > r "Rotor #$(irotor) Tip Radius is set to be less than its Hub Radius. Consider setting the `autoshiftduct` option to true and/or check the input geometry."
    end

    return Rtip, Rhub
end

"""
    generate_blade_elements(B, Omega, rotorzloc, rblade, chords, twists, solidity, airfoils,
        duct_coordinates, centerbody_coordinates, r)

Use the duct and hub geometry to dimensionalize the non-dimensional radial positions in
`rblade`. Then return a [`BladeElements`](@ref) struct.

# Arguments:
- `B::Int64` : number of blades on the rotor
- `Omega::Float` : rotor rotation rate (rad/s)
- `rotorzloc::Float` : x-position of the the rotor (dimensional)
- `rblade::Vector{Float}` : non-dimensional radial positions of each blade element (zero at hub, one at tip)
- `chords::Vector{Float}` : values of the chord lengths of the blade elements
- `twists::Vector{Float}` : values of the twist from the plane of rotation (in radians) of the blade elements
- `airfoils::Vector{Airfoil Function Type}` : Airfoil polar data objects associated with each blade element
- `body_geometry::BodyGeometry` : BodyGeometry object for duct and hub
- `nr::Int64` : desired number of radial stations on the blade
- `r::Vector{Float}`: prescribed dimensional radial stations along the blade (optional)

Future Work: add tip-gap capabilities
"""
function generate_blade_elements(
    B, Omega, rotorzloc, rnondim, chords, twists, airfoils, Rtip, Rhub, rbe; fliplift=false
)

    # get floating point type
    TF = promote_type(
        eltype(B),
        eltype(Omega),
        eltype(rotorzloc),
        eltype(rnondim),
        eltype(chords),
        eltype(twists),
        eltype(Rtip),
        eltype(Rhub),
        eltype(rbe),
    )

    # dimensionalize the blade element radial positions
    rblade = FLOWMath.linear([0.0; 1.0], [0.0; Rtip], rnondim)

    # update chord lengths
    chords = FLOWMath.akima(rblade, chords, rbe)

    # update twists
    twists = FLOWMath.akima(rblade, twists, rbe)

    # update stagger
    stagger = get_stagger(twists)

    # update solidity
    solidity = get_local_solidity(B, chords, rbe)

    # get bounding airfoil polars
    outer_airfoil = similar(airfoils, length(rbe))
    inner_airfoil = similar(airfoils, length(rbe))
    inner_fraction = similar(airfoils, TF, length(rbe))
    for i in 1:(length(rbe))
        # panel radial location
        # ravg = (rwake[i] + rwake[i + 1]) / 2
        # outer airfoil
        io = min(length(rblade), searchsortedfirst(rblade, rbe[i]))
        outer_airfoil[i] = airfoils[io]
        # inner airfoil
        ii = max(1, io - 1)
        inner_airfoil[i] = airfoils[ii]
        # fraction of inner airfoil's polars to use
        if rblade[io] == rblade[ii]
            inner_fraction[i] = 1.0
        else
            inner_fraction[i] = (rbe[i] - rblade[ii]) / (rblade[io] - rblade[ii])
        end

        # Check incorrect extrapolation
        if inner_fraction[i] > 1.0
            inner_fraction[i] = 1.0
        end
    end

    # return blade elements
    return (;
        B,
        Omega,
        rotorzloc,
        rbe,
        chords,
        twists,
        stagger,
        solidity,
        outer_airfoil,
        inner_airfoil,
        inner_fraction,
        Rtip,
        Rhub,
        fliplift,
        # cl=zeros(TF, length(chords)),
        # cd=zeros(TF, length(chords)),
    )
end

function get_local_solidity(B, chord, r)
    return B .* chord ./ (2.0 * pi * r)
end

function get_stagger(twists)
    return 0.5 * pi .- twists
end

# generates rotor panels
function generate_rotor_panels(rotorzloc, rwake)
    x = fill(rotorzloc, length(rwake))
    xr = @views [x'; rwake']

    return generate_panels(xr; isbody=false, isrotor=true)
end

function generate_rotor_panels(
    rotorzloc, wake_grid, rotor_indices_in_wake, nwake_sheets
)
    TF = promote_type(eltype(rotorzloc), eltype(wake_grid))

    xr = [zeros(TF, nwake_sheets, 2) for i in 1:length(rotorzloc)]

    for irotor in 1:length(rotorzloc)
        @views xr[irotor] = [
            fill(rotorzloc[irotor], nwake_sheets)'
            wake_grid[2, rotor_indices_in_wake[irotor], 1:nwake_sheets]'
        ]
    end

    return generate_panels(xr; isbody=false, isrotor=true)
end

"""
Needs to be tested
"""
function generate_rotor_panels!(
    rotor_source_panels, rotorzloc, wake_grid, rotor_indices_in_wake, nwake_sheets
)
    TF = promote_type(eltype(rotorzloc), eltype(wake_grid))

    xr = [zeros(TF, nwake_sheets, 2) for i in 1:length(rotorzloc)]

    for irotor in 1:length(rotorzloc)
        @views xr[irotor] = [
            fill(rotorzloc[irotor], nwake_sheets)'
            wake_grid[2, rotor_indices_in_wake[irotor], 1:nwake_sheets]'
        ]
    end

    return generate_panels!(rotor_source_panels, xr; isbody=false, isrotor=true)
end
