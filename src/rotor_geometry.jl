#=

Functions Related to Rotor Geometry

=#

"""
all dimensional, twist in radians
"""
struct BladeElements
    num_radial_stations
    radial_positions
    chords
    twists
    solidities
    airfoils
    num_blades
    rotor_x_position
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
    num_radial_stations,
    num_blades,
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

    #---------------------------------#
    #       Refine Distribuions       #
    #---------------------------------#
    # Get fine radial positions
    fine_radial_positions = range(Rhub, Rtip; length=num_radial_stations)

    # Chords
    fine_chords = fm.akima(dim_radial_positions, chords, fine_radial_positions)

    # Twists (convert to radians here)
    fine_twists = fm.akima(
        dim_radial_positions, twists .* pi / 180.0, fine_radial_positions
    )

    # - Calculate Solidity - #
    solidities = (2.0 .* pi .* fine_radial_positions) ./ num_blades

    # - Define Airfoil Distribution - #
    # Initialize
    afdist = Array{typeof(airfoils[1])}(undef, num_radial_stations)

    # get average values of radial stations to better place airfoils (center defined airfoils in refined blade)
    mean_rad_stash =
        (dim_radial_positions[1:(end - 1)] .+ dim_radial_positions[2:end]) ./ 2.0

    # loop through refined radial stations and apply appropriate airfoil
    for i in 1:(num_radial_stations)
        ridx = findfirst(x -> x > fine_radial_positions[i], mean_rad_stash)
        if ridx != nothing
            afdist[i] = airfoils[ridx]
        else
            afdist[i] = airfoils[end]
        end
    end

    # - Generate Rotor Panels - #
    rotor_panels = ff.generate_panels(
        method, [rotor_x_position .* ones(num_radial_stations) fine_radial_positions]
    )

    return BladeElements(
        num_radial_stations,
        fine_radial_positions,
        fine_chords,
        fine_twists,
        solidities,
        afdist,
        num_blades,
        rotor_x_position,
    ),
    rotor_panels
end

##TODO: need to update this function
#"""
#    reinterpolate_rotors(wakegrid, rotor, rotoridx)
#Since the wake grid relaxes and is not aligned with aft rotor radial stations, this function reinterpolates rotor data based on updated radial stations.
#(The `rotor` inputs is the only one updated by this function.)
#**Arguments:**
# - `wakegrid::DuctTAPE.WakeGridGeometry` : wake grid geometry object
# - `rotor::DuctTAPE.RotorGeometry` : the rotor geometry to update
# - `rotoridx::Int` : index in the x direction for where the rotor lies on the wake grid
#"""
#function reinterpolate_rotors!(rotors, wakegrid)

#    #check if less than 2 rotors
#    if length(rotors) <= 1
#        return rotors
#    else
#        #if more, go through and update rotors
#        for r in 2:length(rotors)
#            #rename things for convenience
#            nr = length(rotor[r].radialstations)
#            new_rad_stash = wakegrid.x_grid_points[rotoridx[r], :]

#            ## -- Calculate New Section Properties -- ##

#            # update chords
#            new_chords = FLOWMath.Akima(rotor[r].radialstations, rotor[r].chords)
#            for i in 1:nr
#                rotor[r].chords[i] = new_chords(new_rad_stash[i])
#            end

#            # update twists
#            new_twists = FLOWMath.Akima(rotor[r].radialstations, rotor[r].twists)
#            for i in 1:nr
#                rotor[r].twists[i] = new_twists(new_rad_stash[i])
#            end

#            # update solidity
#            new_solidities = FLOWMath.Akima(rotor[r].radialstations, rotor[r].solidities)
#            for i in 1:nr
#                rotor[r].solidities[i] = new_solidities(new_rad_stash[i])
#            end

#            # update radial stations at the end after using both old and new for splines
#            for i in 1:nr
#                rotor[r].radialstations[i] = new_rad_stash[i]
#            end
#        end
#    end

#    return nothing
#end
