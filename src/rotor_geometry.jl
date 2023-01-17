#=

Functions Related to Rotor Geometry

=#

"""
all dimensional, twist in radians
"""
struct BladeElements
    radial_positions
    chords
    twists
    airfoils
end

"""
takes in non-dimensional radial positions, and dimensional everything else.
uses body geometry to determine rotor hub and tip radius
twist in degrees

some day add tip-gap capabilities
"""
function generate_blade_elements(
    rotor_x_position,
    radial_positions,
    chords,
    twists,
    airfoils,
    body_geometry;
    method=ff.AxisymmetricProblem(Vortex(Constant()), Neumann(), [false, true]),
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

    # - Generate Rotor Panels - #

    rotor_panels = ff.generate_panels(
        method, [rotor_x_position .* ones(length(radial_positions)) dim_radial_positions]
    )

    return BladeElements(dim_radial_positions, chords, twists * pi / 180.0, airfoils),
    rotor_panels
end
