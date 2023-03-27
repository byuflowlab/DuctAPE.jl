# """
#     BladeElements{TF, TAF}

# # Fields:
# - `B::Int`: number of blades on the rotor
# - `omega::TF` : rotor rotation rate (rad/s)
# - `xrotor::TF`: x-position of the rotor (dimensionalized)
# - `r::Vector{TF}`: radial coordinate of each blade element (dimensionalized)
# - `chords::Vector{TF}` : chord length of each blade element
# - `twists::Vector{TF}` : twist of each blade element (in radians)
# - `solidities::Vector{TF}` : solidity of each blade element
# - `outer_airfoils::Vector{TAF}` : outer bounding airfoil polar
# - `inner_airfoils::Vector{TAF}` : inner bounding airfoil polar
# - `inner_fraction::Vector{TAF}`: fraction of inner bounding airfoil polar to use
# """
# struct BladeElements{TF,TAF}
#     B::Int
#     omega::TF
#     xrotor::TF
#     r::Vector{TF}
#     chords::Vector{TF}
#     twists::Vector{TF}
#     solidities::Vector{TF}
#     outer_airfoils::Vector{TAF}
#     inner_airfoils::Vector{TAF}
#     inner_fraction::Vector{TF}
# end

"""
    generate_blade_elements(B, omega, xrotor, rblade, chords, twists, solidities, airfoils,
        duct_coordinates, hub_coordinates, r)

Use the duct and hub geometry to dimensionalize the non-dimensional radial positions in
`rblade`. Then return a [`BladeElements`](@ref) struct.

# Arguments:
- `B::Int64` : number of blades on the rotor
- `omega::Float` : rotor rotation rate (rad/s)
- `xrotor::Float` : x-position of the the rotor (dimensional)
- `rblade::Vector{Float}` : non-dimensional radial positions of each blade element (zero at hub, one at tip)
- `chords::Vector{Float}` : values of the chord lengths of the blade elements
- `twists::Vector{Float}` : values of the twist (in degrees) of the blade elements
- `airfoils::Vector{Airfoil Function Type}` : Airfoil polar data objects associated with each blade element
- `body_geometry::BodyGeometry` : BodyGeometry object for duct and hub
- `nr::Int64` : desired number of radial stations on the blade
- `r::Vector{Float}`: prescribed dimensional radial stations along the blade (optional)

Future Work: add tip-gap capabilities
"""
function generate_blade_elements(
    B,
    omega,
    xrotor,
    rblade,
    chords,
    twists,
    solidities,
    airfoils,
    duct_coordinates,
    hub_coordinates,
    xwake,
    rwake,
)

    # get hub and tip wall radial positions
    _, ihub = findmin(x -> abs(x - xrotor), view(hub_coordinates, :, 1))
    _, iduct = findmin(x -> abs(x - xrotor), view(duct_coordinates, :, 1))
    Rhub = hub_coordinates[ihub, 2]
    Rtip = duct_coordinates[iduct, 2]

    # dimensionalize the blade element radial positions
    rblade = fm.linear([Rhub; Rtip], [0.0; 1.0], rblade)

    # update chord lengths
    chords = fm.akima(rblade, chords, rwake)

    # update twists (convert to radians here)
    twists = fm.akima(rblade, twists .* pi / 180.0, rwake)

    # update solidities
    solidities = 2 * pi * rwake / B

    # get bounding airfoil polars
    outer_airfoil = similar(airfoils, length(rwake) - 1)
    inner_airfoil = similar(airfoils, length(rwake) - 1)
    inner_fraction = similar(airfoils, Float64, length(rwake) - 1)
    for i in 1:(length(rwake) - 1)
        # panel radial location
        ravg = (rwake[i] + rwake[i + 1]) / 2
        # outer airfoil
        io = min(length(rblade), searchsortedfirst(rblade, ravg))
        outer_airfoil[io] = airfoils[io]
        # inner airfoil
        ii = max(1, io - 1)
        inner_airfoil[ii] = airfoils[ii]
        # fraction of inner airfoil's polars to use
        inner_fraction[i] = (ravg - rblade[ii]) / (rblade[io] - rblade[ii])
    end

    # find index of the rotor's position in the wake grid
    _, wake_index = findmin(x -> abs(x - xrotor), xwake)

    # return blade elements
    return (;
        B,
        omega,
        xrotor,
        r=rwake,
        chords,
        twists,
        solidities,
        outer_airfoil,
        inner_airfoil,
        inner_fraction,
        wake_index,
    )
end

# generates rotor panels
function generate_rotor_panels(blade_elements, method)
    r = blade_elements.r
    x = fill(blade_elements.xrotor, length(r))
    xr = [x r]

    return ff.generate_panels(method, xr)
end

#TODO: pretty sure this can be deleted.
# # generates dummy rotor panels
# function generate_dummy_rotor_panels(blade_elements, method)
#     r = blade_elements.r
#     x = fill(blade_elements.xrotor, length(r))

#     # place blade elements in center
#     rm = similar(r, length(r) + 1)
#     rm[2:(end - 1)] .= (r[1:(end - 1)] .+ r[2:end]) ./ 2.0
#     rm[1] = 2 * r[1] - rm[2]
#     rm[end] = 2 * r[end] - rm[end - 1]

#     xr = [x, rm]

#     return ff.generate_panels(method, xr)
# end
