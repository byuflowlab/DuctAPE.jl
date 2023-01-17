#=

Functions related to generating and modifying body geometry

Authors: Judd Mehr,

=#

"""
"""
struct BodyGeometry
    duct_inner_spline
    duct_outer_spline
    hub_spline
    duct_range
    hub_range
end

"""
take in duct and hub coordinates and output geometry and panel objects
ASSUMES DUCT COORDINATES ARE PROVIDED STARTING AND ENDING AT TE AND PROGRESSING CLOCKWISE, ASSUMES HUB COORDINATES START AT LE AND END AT TE
"""
function generate_body_geometry(
    duct_coordinates,
    hub_coordinates;
    method=ff.AxisymmetricProblem(Vortex(Constant()), Neumann(), [false, true]),
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
ASSUMES DUCT COORDINATES ARE PROVIDED STARTING AND ENDING AT TE AND PROGRESSING CLOCKWISE, ASSUMES HUB COORDINATES START AT LE AND END AT TE
"""
function generate_splines(duct_coordinates, hub_coordinates)

    # - Break Duct coordinates up - #
    _, leidx = findmin(duct_coordinates[:, 1])
    ndc = length(duct_coordinates[:, 1])

    # - Reverse first half of coordinates for Akima spline inputs - #
    duct_inner_coordinates = reverse(view(duct_coordinates, 1:leidx, :),dims=1)
    duct_outer_coordinates = view(duct_coordinates, leidx:ndc, :)

    # - Create Akima Splines - #
    duct_inner_spline = fm.Akima(duct_inner_coordinates[:, 1], duct_inner_coordinates[:, 2])
    duct_outer_spline = fm.Akima(duct_outer_coordinates[:, 1], duct_outer_coordinates[:, 2])
    hub_spline = fm.Akima(hub_coordinates[:, 1], hub_coordinates[:, 2])

    # - Return - #
    return duct_inner_spline, duct_outer_spline, hub_spline
end
