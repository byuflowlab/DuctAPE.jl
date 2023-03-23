"""
    BladeElements{TF, TAF}

# Fields:
- `B::Int`: number of blades on the rotor
- `omega::TF` : rotor rotation rate (rad/s)
- `xrotor::TF`: x-position of the rotor (dimensionalized)
- `r::Vector{TF}`: radial coordinate of each blade element (dimensionalized)
- `chords::Vector{TF}` : chord length of each blade element
- `twists::Vector{TF}` : twist of each blade element (in radians)
- `solidities::Vector{TF}` : solidity of each blade element
- `airfoils::Vector{TAF}` : Airfoil polars for each blade element
"""
struct BladeElements{TF,TAF}
    B::Int
    omega::TF
    xrotor::TF
    r::Vector{TF}
    chords::Vector{TF}
    twists::Vector{TF}
    solidities::Vector{TF}
    airfoils::Vector{TAF}
end

"""
    generate_blade_elements(B, omega, xrotor, rblade, chords, twists, solidities, airfoils,
        body_geometry, nr, r=range(Rhub, Rtip; length=nr))

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
function generate_blade_elements(B, omega, xrotor, rblade, chords, twists, airfoils, 
    body_geometry, nr, r=nothing)

    # Make sure that the rotor is inside the duct
    @assert xrotor > body_geometry.duct_range[1] "Rotor is in front of duct leading edge."
    @assert xrotor < body_geometry.duct_range[2] "Rotor is behind duct trailing edge."
    @assert xrotor > body_geometry.hub_range[1] "Rotor is in front of hub leading edge."
    @assert xrotor < body_geometry.hub_range[2] "Rotor is behind hub trailing edge."

    # get duct and hub wall radial positions
    Rtip = body_geometry.duct_inner_spline(xrotor)
    Rhub = body_geometry.hub_spline(xrotor)

    # dimensionalize the blade element radial positions
    rblade = fm.linear([Rhub; Rtip], [0.0; 1.0], rblade)

    # update radial coordinates
    if isnothing(r)
        r = collect(range(Rhub, Rtip; length=nr))
    end

    # update chord lengths
    chords = fm.akima(rblade, chords, r)

    # update twists (convert to radians here)
    twists = fm.akima(rblade, twists .* pi/180.0, r)

    # update solidities
    solidities = 2*pi*r/B

    # update airfoil polars (using midpoint between stations), note that this will create 
    # a discontinuity if the radial positions are updated during an optimization
    ravg = (dim_radial_positions[1:(end - 1)] .+ dim_radial_positions[2:end]) ./ 2.0
    afdist = similar(airfoils, nr)
    for i in 1:nr
        ridx = findfirst(x -> x > updated_radial_positions[i], ravg)
        afdist[i] = isnothing(ridx) ? airfoils[end] : airfoils[ridx]
    end

    # return blade elements
    return BladeElements(B, omega, xrotor, r, chords, twists, solidities, afdist)
end


# generates rotor panels
function generate_rotor_panels(blade_elements, method)
    
    r = blade_elements.r
    x = fill(blade_elements.xrotor, length(r))
    xr = [x r]

    return ff.generate_panels(method, xr)
end

# generates dummy rotor panels
function generate_dummy_rotor_panels(blade_elements, method)

    r = blade_elements.r
    x = fill(blade_elements.xrotor, length(r))

    # place blade elements in center
    rm = similar(r, length(r)+1)
    rm[2:end-1] .= (r[1:end-1] .+ r[2:end]) ./ 2.0
    rm[1] = 2*r[1] - rm[2]
    rm[end] =  2*r[end] - rm[end-1]

    xr = [x, rm]

    return ff.generate_panels(method, xr)
end
