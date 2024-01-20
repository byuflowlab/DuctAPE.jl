#=
Attempt at a "1D" model in order to provide a reasonable starting point for design/optimization

Goal for option 1: input open rotor design, output duct inlet and outlet areas
Goal for option 2: input some performance metric (say, thrust), and some sizing guide (say rotor diameter), and some operating conditions (say, Vinf and Omega), and output reasonable duct inlet and outlet as well as some rotor geometry (chord, twist, camber?)

=#

import FLOWMath as fm

######################################################################
#                                                                    #
#                         1D Model Option 1                          #
#                                                                    #
######################################################################

"""

# Arguments:
- `thrust::Float` : desired thrust force.
- `exit_area::Float` : exit area of nozzle
- `Vinf::Float` : freestream velocity
- `Rtip::Float` : rotor tip radius
- `Rhub::Float` : rotor hub radius
- `radii:Vector{Float}` : radial positions along the blade (dimensional)
- `num_blades::Int` : number of blades
- `lift_polars::Vector{Matrix{Float}}` : matrix [α cℓ] of angles of attack vs lift coefficients for desired airfoils
-

# Keyword Arguments:
- `rho::Float` : air density (ρ), assumed constant = 1.225 kg/m3
- `ambient_static_pressure::Float` : ambient static pressure, default = 101325.0 Pa
- `ambient_static_temperature::Float` : ambient static temperature, default = 288.15 K
- `flow_coeff::Float` : tip flow coefficient (ϕ), default = 0.4
- `adiabatic_stage_efficiency::Float` : adiabatic stage efficiency, default = 0.85
- `c_p::Float` : specific heat at constant pressure, default is 1.005 kJ/kg-K
- `lift_coefficients::Vector{Float}` : desired operating lift coefficient for each blade section, default = [0.6] (constant along blade)
- `verbose::Bool` : if true, will print out warning if solidity exceeds 0.3
- `debug::Bool` : if true, function returns intermediate values in a named tuple

# Returns:
- `chords::Vector{Float}` : vector of chord distribution
- `twists::Vector{Float}` : vector of twist distribution
- `debug::NamedTuple` : NamedTuple of intermediate values
"""
function opt_prelim(
    Vinf,
    thrust,
    exit_area,
    Rtip,
    Rhub,
    radii,
    num_blades,
    lift_polars;
    rho=1.225,
    flow_coeff=0.4,
    ambient_static_pressure=101325.0,
    ambient_static_temperature=288.15,
    adiabatic_stage_efficiency=0.85,
    c_p=1.005,
    lift_coefficients=[0.6],
    verbose=true,
)

    # Compute Exit Velocity m/s
    Vexit = compute_exit_velocity(Vinf, thrust, exit_area; rho=rho)

    # Compute Fan Face Axial Velocity m/s
    fan_face_axial_velocity = compute_fan_face_axial_velocity(Vexit, exit_area, Rtip, Rhub)

    # Compute Tip Rotational Velocity m/s
    tip_rotational_velocity = compute_tip_rotational_velocity(
        fan_face_axial_velocity; flow_coeff=flow_coeff
    )

    # Compute Exhaust Total Pressue Pa = N/m2 = kg/m-s2
    exhaust_total_pressue = compute_exhaust_total_pressure(
        Vexit; rho=rho, ambient_static_pressure=ambient_static_pressure
    )

    # Compute the Fan Pressure Ratio (non-dimensional)
    pressure_ratio = compute_fan_pressure_ratio(
        exhaust_total_pressue,
        Vinf;
        rho=rho,
        ambient_static_pressure=ambient_static_pressure,
    )

    # Compute_fan_temperature_ratio (non-dimensional)
    temperature_ratio = compute_fan_temperature_ratio(
        pressure_ratio; gamma=1.4, adiabatic_stage_efficiency=adiabatic_stage_efficiency
    )

    # Compute Total Specific Enthalpy Rise kJ/kg
    enthalpy_rise = compute_total_specific_enthalpy_rise(
        temperature_ratio,
        Vinf;
        c_p=c_p,
        ambient_static_temperature=ambient_static_temperature,
    )

    # Compute the Required Work Coefficient
    work_coeff = compute_work_coefficient(enthalpy_rise * 1000.0, tip_rotational_velocity)

    # Compute Change in Swirl Velocity at Tip
    change_in_swirl = compute_change_in_swirl(work_coeff, tip_rotational_velocity)

    # Compute Inflow Angles
    inflow_angles = compute_inflow_angles(
        fan_face_axial_velocity, tip_rotational_velocity, Rtip, radii
    )

    # Compute Swirl Velocity Distribution
    swirl_velocities = compute_swirl_velocity_distribution(Rtip, radii, change_in_swirl)

    # Look up Section Angles of Attack
    angles_of_attack = look_up_angle_of_attack(
        lift_polars, lift_coefficients, length(radii)
    )

    # Calculate Stagger Angles
    stagger_angles = calculate_stagger_angles(inflow_angles, angles_of_attack)

    # Calculate Circulation Distribution
    circulations = calculate_circulation(change_in_swirl, radii, num_blades)

    # Calculate Chord
    chords = calculate_chord(
        circulations, lift_coefficients, fan_face_axial_velocity, change_in_swirl
    )

    # Check Solidity
    if verbose
        check_solidity(chords, radii, num_blades)
    end

    return chords,
    pi / 2.0 .- stagger_angles,
    (;
        circulations,
        stagger_angles,
        angles_of_attack,
        swirl_velocities,
        inflow_angles,
        omega=compute_omega(tip_rotational_velocity, Rtip),
        change_in_swirl,
        work_coeff,
        enthalpy_rise,
        temperature_ratio,
        pressure_ratio,
        exhaust_total_pressue,
        tip_rotational_velocity,
        fan_face_axial_velocity,
        fan_annulus_area=compute_annulus_area(Rtip, Rhub),
        Vexit,
    )
end

"""

Apply conservation of mass and momentum to obtain exit velocity.
"""
function compute_exit_velocity(Vinf, thrust, exit_area; rho=1.225)
    return Vinf + 0.5 * sqrt(2.0 * Vinf + 4.0 * thrust / (rho * exit_area))
end

"""
"""
function compute_annulus_area(Rtip, Rhub)
    return pi * (Rtip^2 - Rhub^2)
end

"""
"""
function compute_fan_face_axial_velocity(Vexit, exit_area, Rtip, Rhub)

    # Get fan annulus area
    fan_area = compute_annulus_area(Rtip, Rhub)

    # convervation of mass
    return Vexit * exit_area / fan_area
end

"""
"""
function compute_tip_rotational_velocity(Cz; flow_coeff=0.4)
    return Cz / flow_coeff
end

"""
"""
function compute_exhaust_total_pressure(Vexit; rho=1.225, ambient_static_pressure=101325.0)
    return 0.5 * rho * Vexit^2 + ambient_static_pressure
end

"""
"""
function compute_fan_pressure_ratio(pte, Vinf; rho=1.225, ambient_static_pressure=101325.0)
    ptinf = 0.5 * rho * Vinf^2 + ambient_static_pressure

    return pte / ptinf
end

"""
"""
function compute_fan_temperature_ratio(
    pressure_ratio; gamma=1.4, adiabatic_stage_efficiency=0.85
)
    return (pressure_ratio^((gamma - 1) / gamma) - 1.0) / adiabatic_stage_efficiency + 1.0
end

"""
"""
function compute_total_specific_enthalpy_rise(
    temperature_ratio, Vinf; c_p=1.005, ambient_static_temperature=288.15
)
    return c_p * (temperature_ratio - 1.0) * ambient_static_temperature
end

"""
"""
function compute_work_coefficient(enthalpy_rise, tip_rotational_velocity)
    return enthalpy_rise / tip_rotational_velocity^2
end

"""
"""
function compute_change_in_swirl(work_coeff, tip_rotational_velocity)
    return work_coeff * tip_rotational_velocity
end

"""
"""
function compute_omega(tip_rotational_velocity, Rtip)
    return (tip_rotational_velocity / Rtip)
end

"""

radians
"""
function compute_inflow_angles(
    fan_face_axial_velocity, tip_rotational_velocity, Rtip, radii
)
    omega = compute_omega(tip_rotational_velocity, Rtip)

    return atan.(omega * radii, fan_face_axial_velocity)
end

"""
"""
function compute_swirl_velocity_distribution(Rtip, radii, change_in_swirl)
    return Rtip * change_in_swirl ./ radii
end

"""
"""
function look_up_angle_of_attack(lift_polars, lift_coefficients, ns)

    # Initialize Output
    angles_of_attack = zeros(eltype(lift_coefficients), length(lift_coefficients))
    # Loop through sections
    for (is, cl) in enumerate(lift_coefficients)

        # get entry just before
        i1 = searchsortedlast(lift_polars[is][:, 2], cl)

        # get entry just after
        i2 = searchsortedfirst(lift_polars[is][:, 2], cl)

        # interpolate
        angles_of_attack[is] = fm.linear(
            [lift_polars[is][i1, 2]; lift_polars[is][i2, 2]],
            [lift_polars[is][i1, 1]; lift_polars[is][i2, 1]],
            cl,
        )
    end

    if length(angles_of_attack) == 1
        return fill(angles_of_attack[1], ns)
    else
        return angles_of_attack
    end
end

"""
"""
function calculate_stagger_angles(inflow_angles, angles_of_attack)
    return inflow_angles .- angles_of_attack
end

"""
"""
function calculate_circulation(change_in_swirl, radii, num_blades)
    return 2.0 * pi * change_in_swirl .* radii / num_blades
end

"""
"""
function calculate_chord(
    circulations, lift_coefficients, fan_face_axial_velocity, change_in_swirl
)
    average_relative_velocity = sqrt(fan_face_axial_velocity^2 + (change_in_swirl / 2.0)^2)

    return 2.0 * circulations ./ (lift_coefficients * average_relative_velocity)
end

function check_solidity(chords, radii, num_blades)
    solidity = chords * num_blades ./ (2.0 * pi * radii)
    warnid = findall(x -> x > 0.3, solidity)
    if !isempty(warnid)
        @warn "Solidity greater than 0.3 at sections $(warnid)"
    end
    return nothing
end
######################################################################
#                                                                    #
#                              SCRATCH                               #
#                                                                    #
######################################################################

# """
# """
# function compute_exit_flow_angles()
# end

#
#
#
#
#
#

#=
#-----------------------------------#
# Potentially useful relationships  #
#-----------------------------------#
Thrust = mdot * (exit_velocity - freestream_velocity) (conservation of momentum)
=#

#=

Known (or can be approximated): Rtip, Rhub, chord, twist, airfoil parameters, reynolds number, rotation rate, number of blades, open rotor performance data (CT, c_p, CQ, eta), freestream conditions (Vinf, rhoinf, muinf, asound).

Want to Find: inlet area, outlet area for reasonable ducted performance (at least as much thrust and efficiency as without the duct if possible)

=#

"""
not sure what's happening here, but seems like we might need something along these lines
"""
function pressure_from_open_ct(CT, disk_area)
    T = CT * rho * n^2 * disk_area^4 # is this valid
    dp = T / disk_area

    return dp
end

"""
function to call the other functions and find what reasonable inlet and outlet areas might be
"""
function find_io_areas(
    Rtip, # rotor tip radius
    Rhub, # rotor hub radius
    radial_stations, # radial stations between Rhub and Rtip
    chord, # chord distribution
    twist, # twist distribution (better for this to be stagger?)
    ambient_static_pressure=101325.0, # Pa standard at sea level
)

    # solver variable initial guesses, just set to disk annulus area for now
    inlet_area = pi * (Rtip^2 - Rhub^2)
    outlet_area = inlet_area

    #TODO: choose a reasonable solver approach (NLSolve?)

    return res.zero
end

######################################################################
#                                                                    #
#                         1D Model Option 2                          #
#                                                                    #
######################################################################
#######################################################################
#                                                                    #
######################################################################

