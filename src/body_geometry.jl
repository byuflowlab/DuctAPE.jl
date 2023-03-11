#=

Functions related to generating and modifying body geometry

Authors: Judd Mehr,

=#

"""
    BodyGeometry

**Fields:**
- `duct_inner_spline::FLOWMath.Akima` : Akima spline object describing the duct geometry inner surface from the leading to trailing edge.
- `duct_outer_spline::FLOWMath.Akima` : Akima spline object describing the duct geometry outer surface from the leading to trailing edge.
- `hub_spline::FLOWMath.Akima` : Akima spline object describing the hub geometry from the leading to trailing edge.
- `duct_range::Vector{Int}` : Vector containing the leading and trailing x-locations for the duct geometry.
- `hub_range::Vector{Int}` : Vector containing the leading and trailing x-locations for the hub geometry.
"""
struct BodyGeometry{TF,TX,TY,TC}
    duct_inner_spline::FLOWFoil.Akima{TX,TY,TC}
    duct_outer_spline::FLOWFoil.Akima{TX,TY,TC}
    hub_spline::FLOWFoil.Akima{TX,TY,TC}
    duct_range::Vector{TF}
    hub_range::Vector{TF}
end

"""
    generate_body_geometry(duct_coordinates, hub_coordinates; kwargs)

Generates BodyGeometry and duct and hub panel objects from the input coordinates.

Note: assumes duct coordinates are provided starting and ending at te and progressing clockwise, assumes hub coordinates start at leading edge and end at trailing edge.

**Arguments:**
- `duct_coordinates::Array{Float64,2}` : [x y] coordinates of duct
- `hub_coordinates::Array{Float64,2}` : [x y] coordinates of hub

**Keyword Arguments:**
- `method::FLOWFoil.AxisymmetricProblem` : default = AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true]),

**Returns:**
- `body_geometry::BodyGeometry` : BodyGeometry object for duct and hub
- `body_panels::Vector{FLOWFoil.AxisymmetricPanel}` : Vector of flowfoil panel objects for the duct and hub
"""
function generate_body_geometry(
    duct_coordinates,
    hub_coordinates;
    method=ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true]),
)

    # - Use FLOWFoil Functions to generate array of panel objects for duct and hub. - #
    body_panels = ff.generate_panels(method, (duct_coordinates, hub_coordinates))

    # - Generate Akima Splines - #
    duct_inner_spline, duct_outer_spline, hub_spline = generate_splines(
        duct_coordinates, hub_coordinates
    )

    # - Return DuctTAPE specific geometry and FLOWFoil panel objects. - #
    return BodyGeometry(
        duct_inner_spline,
        duct_outer_spline,
        hub_spline,
        [minimum(duct_coordinates[:, 1]); maximum(duct_coordinates[:, 1])],
        [minimum(hub_coordinates[:, 1]); maximum(hub_coordinates[:, 1])],
    ),
    body_panels
end

"""
    generate_splines(duct_coordinates, hub_coordinates)

Assembles Akima spline objects for the given duct and hub coordinates.

Note: assumes duct coordinates are provided starting and ending at te and progressing clockwise, assumes hub coordinates start at leading edge and end at trailing edge.

**Arguments:**
- `duct_coordinates::Array{Float64,2}` : [x y] coordinates of duct
- `hub_coordinates::Array{Float64,2}` : [x y] coordinates of hub

**Returns:**
- `duct_inner_spline::FLOWMath.Akima` : Akima spline object describing the duct geometry inner surface from the leading to trailing edge.
- `duct_outer_spline::FLOWMath.Akima` : Akima spline object describing the duct geometry outer surface from the leading to trailing edge.
- `hub_spline::FLOWMath.Akima` : Akima spline object describing the hub geometry from the leading to trailing edge.
"""
function generate_splines(duct_coordinates, hub_coordinates)

    # - Break Duct coordinates up - #
    _, leidx = findmin(duct_coordinates[:, 1])
    ndc = length(duct_coordinates[:, 1])

    # - Reverse first half of coordinates for Akima spline inputs - #
    duct_inner_coordinates = reverse(view(duct_coordinates, 1:leidx, :); dims=1)
    duct_outer_coordinates = view(duct_coordinates, leidx:ndc, :)

    # - Create Akima Splines - #
    duct_inner_spline = fm.Akima(duct_inner_coordinates[:, 1], duct_inner_coordinates[:, 2])
    duct_outer_spline = fm.Akima(duct_outer_coordinates[:, 1], duct_outer_coordinates[:, 2])
    hub_spline = fm.Akima(hub_coordinates[:, 1], hub_coordinates[:, 2])

    # - Return - #
    return duct_inner_spline, duct_outer_spline, hub_spline
end
