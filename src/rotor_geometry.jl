#=

Functions Related to Rotor Geometry

Authors: Judd Mehr,

=#

"""
    BladeElements

**Fields:**
- `num_radial_stations::Int64` : number of radial stations on the blade
- `radial_positions::Vector{Float}` : radial positions (dimensional) of the blade elements on the blade
- `chords::Vector{Float}` : values of the chord lengths of the blade elements
- `twists::Vector{Float}` : values of the twist (in radians) of the blade elements
- `solidities::Vector{Float}` : values of the local solidities of the blade elements
- `airfoils::Vector{Airfoil Function Type}` : Airfoil polar data objects associated with each blade element
- `num_blades::Int64` : number of blades on the rotor
- `rotor_x_position::Float` : x-position (dimensional) of the the rotor
- `omega::Float` : Rotation rate, in radians per second, of rotor
"""
struct BladeElements{TF,TAF}
    num_radial_stations::Vector{Int64}
    radial_positions::Vector{TF}
    chords::Vector{TF}
    twists::Vector{TF}
    solidities::Vector{TF}
    airfoils::Vector{TAF}
    num_blades::Vector{Int64}
    rotor_x_position::Vector{TF}
    omega::Vector{TF}
end

"""
    generate_blade_elements(
    rotor_x_position,
    radial_positions,
    chords,
    twists,
    airfoils,
    num_radial_stations,
    num_blades,
    body_geometry;
    kwargs
)

takes in non-dimensional radial positions, and dimensional everything else.
uses body geometry to determine rotor hub and tip radius
twist in degrees

**Arguments:**
- `rotor_x_position::Float` : x-position (dimensional) of the the rotor
- `radial_positions::Vector{Float}` : radial positions (non-dimensional, zero at hub, one at tip) of the blade elements on the blade
- `chords::Vector{Float}` : values of the chord lengths of the blade elements
- `twists::Vector{Float}` : values of the twist (in degrees) of the blade elements
- `solidities::Vector{Float}` : values of the local solidities of the blade elements
- `airfoils::Vector{Airfoil Function Type}` : Airfoil polar data objects associated with each blade element
- `num_radial_stations::Int64` : number of radial stations on the blade, if more than the length of the inputs, things will be interpolated to refine the blade.
- `num_blades::Int64` : number of blades on the rotor
- `omega::Float` : Rotation rate, in radians per second, of rotor
- `body_geometry::BodyGeometry` : BodyGeometry object for duct and hub

**Keyword Arguments:**
- `method::FLOWFoil.AxisymmetricProblem` : default = AxisymmetricProblem(Vortex(Constant()), Neumann(), [false, true]),

Future Work: add tip-gap capabilities
"""
function generate_blade_elements(
    rotor_x_position,
    radial_positions,
    chords,
    twists,
    airfoils,
    num_radial_stations,
    num_blades,
    omega,
    body_geometry;
    updated_radial_positions=nothing,
)

    # - Find Tip Radius based on body geometry - #

    # Make sure that the rotor is inside the duct
    @assert rotor_x_position > body_geometry.duct_range[1] "Rotor is in front of duct leading edge."
    @assert rotor_x_position < body_geometry.duct_range[2] "Rotor is behind duct trailing edge."
    @assert rotor_x_position > body_geometry.hub_range[1] "Rotor is in front of hub leading edge."
    @assert rotor_x_position < body_geometry.hub_range[2] "Rotor is behind hub trailing edge."

    # Sample the splines to get the radius values
    Rtip = body_geometry.duct_inner_spline(rotor_x_position)
    Rhub = body_geometry.hub_spline(rotor_x_position)

    # - Dimensionalize the blade element radial positions - #
    dim_radial_positions = lintran([Rhub; Rtip], [0.0; 1.0], radial_positions)

    #---------------------------------#
    #       Refine Distribuions       #
    #---------------------------------#
    # Get fine radial positions
    if updated_radial_positions == nothing
        updated_radial_positions = collect(range(Rhub, Rtip; length=num_radial_stations))
    end

    # Chords
    fine_chords = fm.akima(dim_radial_positions, chords, updated_radial_positions)

    # Twists (convert to radians here)
    fine_twists = fm.akima(
        dim_radial_positions, twists .* pi / 180.0, updated_radial_positions
    )

    # - Calculate Solidity - #
    solidities = (2.0 .* pi .* updated_radial_positions) ./ num_blades

    # - Define Airfoil Distribution - #
    # Initialize
    afdist = Array{typeof(airfoils[1])}(undef, num_radial_stations)

    # get average values of radial stations to better place airfoils (center defined airfoils in refined blade)
    mean_rad_stash =
        (dim_radial_positions[1:(end - 1)] .+ dim_radial_positions[2:end]) ./ 2.0

    # loop through refined radial stations and apply appropriate airfoil
    for i in 1:(num_radial_stations)
        ridx = findfirst(x -> x > updated_radial_positions[i], mean_rad_stash)
        if ridx != nothing
            afdist[i] = airfoils[ridx]
        else
            afdist[i] = airfoils[end]
        end
    end

    return BladeElements(
        [num_radial_stations],
        updated_radial_positions,
        fine_chords,
        fine_twists,
        solidities,
        afdist,
        [num_blades],
        [rotor_x_position],
        [omega],
    )
end

function generate_blade_elements!(
    blade_elements,
    rotor_x_position,
    radial_positions,
    chords,
    twists,
    airfoils,
    num_radial_stations,
    num_blades,
    omega,
    body_geometry;
    updated_radial_positions=nothing,
)

    # - Find Tip Radius based on body geometry - #

    # Make sure that the rotor is inside the duct
    @assert rotor_x_position > body_geometry.duct_range[1] "Rotor is in front of duct leading edge."
    @assert rotor_x_position < body_geometry.duct_range[2] "Rotor is behind duct trailing edge."
    @assert rotor_x_position > body_geometry.hub_range[1] "Rotor is in front of hub leading edge."
    @assert rotor_x_position < body_geometry.hub_range[2] "Rotor is behind hub trailing edge."

    # Sample the splines to get the radius values
    Rtip = body_geometry.duct_inner_spline(rotor_x_position)
    Rhub = body_geometry.hub_spline(rotor_x_position)

    # - Dimensionalize the blade element radial positions - #
    dim_radial_positions = lintran([Rhub; Rtip], [0.0; 1.0], radial_positions)

    #---------------------------------#
    #       Refine Distribuions       #
    #---------------------------------#
    # Get fine radial positions
    if updated_radial_positions == nothing
        blade_elements.radial_positions .= collect(
            range(Rhub, Rtip; length=num_radial_stations)
        )
    else
        blade_elements.radial_positions .= updated_radial_positions
        num_radial_stations = length(updated_radial_positions)
    end

    # Chords
    blade_elements.chords .= fm.akima(
        dim_radial_positions, chords, blade_elements.radial_positions
    )

    # Twists (convert to radians here)
    blade_elements.twists .= fm.akima(
        dim_radial_positions, twists .* pi / 180.0, blade_elements.radial_positions
    )

    # - Calculate Solidity - #
    blade_elements.solidities .=
        (2.0 .* pi .* blade_elements.radial_positions) ./ num_blades

    # - Define Airfoil Distribution - #

    # get average values of radial stations to better place airfoils (center defined airfoils in refined blade)
    mean_rad_stash =
        (dim_radial_positions[1:(end - 1)] .+ dim_radial_positions[2:end]) ./ 2.0

    # loop through refined radial stations and apply appropriate airfoil
    for i in 1:(num_radial_stations)
        ridx = findfirst(x -> x > blade_elements.radial_positions[i], mean_rad_stash)
        if ridx != nothing
            blade_elements.airfoils[i] = airfoils[ridx]
        else
            blade_elements.airfoils[i] = airfoils[end]
        end
    end

    # - Update the rest - #
    blade_elements.rotor_x_position[1] = rotor_x_position
    blade_elements.omega[1] = omega
    blade_elements.num_blades[1] = num_blades

    return nothing
end

function initialize_blade_elements(rotor_parameters, body_geometry, num_rotors)
    return [
        generate_blade_elements(
            rotor_parameters.rotor_x_position,
            rotor_parameters.radial_positions,
            rotor_parameters.chords,
            rotor_parameters.twists,
            rotor_parameters.airfoils,
            rotor_parameters.num_radial_stations,
            rotor_parameters.num_blades,
            rotor_parameters.omega,
            body_geometry;
            updated_radial_positions=nothing,
        ) for i in 1:num_rotors
    ]
end

"""
"""
function initialize_rotor_panels(
    blade_elements,
    num_rotors;
    method=ff.AxisymmetricProblem(Source(Constant()), Neumann(), [true]),
)
    return [
        ff.generate_panels(
            method,
            [blade_elements.rotor_x_position[1] .*
             ones(blade_elements.num_radial_stations[1]) blade_elements.radial_positions],
        ) for i in 1:num_rotors
    ]
end

"""
"""
function initialize_dummy_rotor_panels(
    blade_elements,
    num_rotors;
    method=ff.AxisymmetricProblem(Source(Constant()), Neumann(), [true]),
)
    TF = eltype(blade_elements.chords)

    # - Set up the panel edges such that the blade elements are at the centers - #
    edges = zeros(TF, blade_elements.num_radial_stations[1] + 1)

    rp = view(blade_elements.radial_positions, :)

    edges[2:(end - 1)] = (rp[1:(end - 1)] .+ rp[2:end]) / 2.0

    edges[1] = 2.0 * rp[1] - edges[2]
    edges[end] = 2.0 * rp[end] - edges[end - 1]

    return [
        ff.generate_panels(
            method,
            [blade_elements.rotor_x_position[1] .*
             ones(blade_elements.num_radial_stations[1] + 1) edges],
        ) for i in 1:num_rotors
    ]
end

"""
"""
function update_dummy_rotor_panels!(
    blade_elements,
    dummy_panels,
    r_grid_points;
    method=ff.AxisymmetricProblem(Source(Constant()), Neumann(), [true]),
)

    # - Set up the panel edges such that the blade elements are at the centers - #
    TF = eltype(dummy_panels.panel_center)
    edges = zeros(TF, blade_elements.num_radial_stations[1] + 1)

    rp = view(r_grid_points, :)

    edges[2:(end - 1)] = (rp[1:(end - 1)] .+ rp[2:end]) / 2.0

    edges[1] = 2.0 * rp[1] - edges[2]
    edges[end] = 2.0 * rp[end] - edges[end - 1]

    ff.generate_panels!(
        method,
        dummy_panels,
        [blade_elements.rotor_x_position[1] .*
         ones(blade_elements.num_radial_stations[1] + 1) edges],
    )

    return nothing
end
