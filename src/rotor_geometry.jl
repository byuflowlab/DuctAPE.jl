# """
#     BladeElements{TF, TAF}

# # Fields:
# - `B::Int`: number of blades on the rotor
# - `Omega::TF` : rotor rotation rate (rad/s)
# - `xrotor::TF`: x-position of the rotor (dimensionalized)
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
#     xrotor::TF
#     r::Vector{TF}
#     chords::Vector{TF}
#     twists::Vector{TF}
#     solidity::Vector{TF}
#     outer_airfoils::Vector{TAF}
#     inner_airfoils::Vector{TAF}
#     inner_fraction::Vector{TF}
# end

"""
    generate_blade_elements(B, Omega, xrotor, rblade, chords, twists, solidity, airfoils,
        duct_coordinates, hub_coordinates, r)

Use the duct and hub geometry to dimensionalize the non-dimensional radial positions in
`rblade`. Then return a [`BladeElements`](@ref) struct.

# Arguments:
- `B::Int64` : number of blades on the rotor
- `Omega::Float` : rotor rotation rate (rad/s)
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
    B, Omega, xrotor, rnondim, chords, twists, airfoils, Rtip, Rhub, rbe
)

    # get floating point type
    TF = promote_type(
        eltype(chords),
        eltype(twists),
        eltype(Omega),
        eltype(rbe),
        eltype(Rtip),
    )

    # dimensionalize the blade element radial positions
    rblade = fm.linear([0.0; 1.0], [0.0; Rtip], rnondim)

    # update chord lengths
    chords = fm.akima(rblade, chords, rbe)

    # update twists
    twists = fm.akima(rblade, twists, rbe)

    # update solidity
    solidity = B * chords ./ (2.0 * pi * rbe)

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

    # find index of the rotor's position in the wake grid
    # _, wake_index = findmin(x -> abs(x - xrotor), xwake)

    # return blade elements
    return (;
        B,
        Omega,
        xrotor,
        rbe,
        chords,
        twists,
        solidity,
        outer_airfoil,
        inner_airfoil,
        inner_fraction,
        Rtip,
        Rhub,
        # wake_index,
    )
end

# generates rotor panels
function generate_rotor_panels(
    xrotor, rwake, method=ff.AxisymmetricProblem(Source(Constant()), Dirichlet(), [true])
)
    x = fill(xrotor, length(rwake))
    xr = [x rwake]

    return ff.generate_panels(method, xr)
end
