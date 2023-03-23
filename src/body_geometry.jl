"""
    BodyGeometry

Contains splines of the duct and hub geometry.

# Fields:
- `duct_inner_spline::FLOWMath.Akima`: Spline of duct inner surface from LE to TE
- `duct_outer_spline::FLOWMath.Akima`: Spline of duct outer surface from LE to TE
- `hub_spline::FLOWMath.Akima` : Spline of hub geometry from the LE to TE
- `duct_range::Vector{Int}`: Indices of duct leading and trailing x-locations.
- `hub_range::Vector{Int}`: Indices of hub leading and trailing x-locations.
"""
struct BodyGeometry{TF,TX,TY,TC}
    duct_inner_spline::FLOWFoil.Akima{TX,TY,TC}
    duct_outer_spline::FLOWFoil.Akima{TX,TY,TC}
    hub_spline::FLOWFoil.Akima{TX,TY,TC}
    duct_range::Vector{TF}
    hub_range::Vector{TF}
end

"""
    generate_body_geometry(duct_coordinates, hub_coordinates; kwargs...)

Define the geometry (see [`BodyGeometry`](@ref)) and paneling (see [`FLOWFoil.AxisymmetricPanel`](@ref)) 
of the duct and center body using their coordinates.  

# Arguments:
- `duct_coordinates::Array{Float64,2}` : duct coordinates, starting from trailing edge, going clockwise
- `hub_coordinates::Array{Float64,2}` : hub coordinates, starting from leading edge

# Keyword Arguments:
- `method` : Axisymmetric method as defined by FLOWFoil, defaults to `AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true])`

# Returns:
- `body_geometry`: Splines of the center body and duct geometry, of type [`BodyGeometry`](@ref)
- `body_panels`: Vector containing duct and center body paneling, each of type [`FLOWFoil.AxisymmetricPanel`](@ref)
"""
function generate_body_geometry(duct_coordinates, hub_coordinates;
    method=ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true]))

    # generate splines
    duct_inner_spline, duct_outer_spline, hub_spline = generate_splines(duct_coordinates, hub_coordinates)
    duct_range = [minimum(duct_coordinates[:, 1]); maximum(duct_coordinates[:, 1])]
    hub_range = [minimum(hub_coordinates[:, 1]); maximum(hub_coordinates[:, 1])]

    # define body geometry
    body_geometry = BodyGeometry(duct_inner_spline, duct_outer_spline, hub_spline,
        duct_range, hub_range),

    # define duct and center body paneling using FLOWFoil
    body_panels = ff.generate_panels(method, (duct_coordinates, hub_coordinates))

    return body_geometry, body_panels
end

"""
    generate_splines(duct_coordinates, hub_coordinates)

Assembles Akima spline objects for the given duct and hub coordinates.

# Arguments:
- `duct_coordinates::Array{Float64,2}` : duct coordinates, starting from trailing edge, going clockwise
- `hub_coordinates::Array{Float64,2}` : hub coordinates, starting from leading edge

# Returns:
- `duct_inner_spline::FLOWMath.Akima`: Spline of duct inner surface from LE to TE
- `duct_outer_spline::FLOWMath.Akima`: Spline of duct outer surface from LE to TE
- `hub_spline::FLOWMath.Akima` : Spline of hub geometry from the LE to TE
"""
function generate_splines(duct_coordinates, hub_coordinates)

    # split duct into inner and outer coordinates
    _, leidx = findmin(duct_coordinates[:, 1])
    ndc = length(duct_coordinates[:, 1])
    duct_inner_coordinates = reverse(view(duct_coordinates, 1:leidx, :); dims=1)
    duct_outer_coordinates = view(duct_coordinates, leidx:ndc, :)

    # construct splines of duct and hub geometry
    duct_inner_spline = fm.Akima(duct_inner_coordinates[:, 1], duct_inner_coordinates[:, 2])
    duct_outer_spline = fm.Akima(duct_outer_coordinates[:, 1], duct_outer_coordinates[:, 2])
    hub_spline = fm.Akima(hub_coordinates[:, 1], hub_coordinates[:, 2])

    # return result
    return duct_inner_spline, duct_outer_spline, hub_spline
end
