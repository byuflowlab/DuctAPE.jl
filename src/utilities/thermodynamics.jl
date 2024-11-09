"""
    sa1(altitude; hardness=50)

Standard atmosphere temperature and pressure in SI units blended between the first linear portion and the constant portion.
"""
function sa1(altitude; hardness=50)
    # Sigmoid blend of fits for <11000 and <25000
    T1(altitude) = 15.04 - 0.00649 * altitude + 273.1 #units: K
    P1(altitude) = 101.29 * (T1(altitude) / 288.08)^5.256 #units: kPa
    T2(altitude) = -56.46 + 273.1 #units: K
    P2(altitude) = 22.65 * exp(1.73 - 0.000157 * altitude) #units: kPa

    return FLOWMath.sigmoid_blend(T1(altitude), T2(altitude), altitude, 11000, hardness),
    FLOWMath.sigmoid_blend(P1(altitude), P2(altitude), altitude, 11000, hardness)
end

"""
    sa2(altitude; hardness=50)

Standard atmosphere temperature and pressure in SI units blended between the the constant portion and the second linear portion.
"""
function sa2(altitude; hardness=50)
    # Sigmoid blend of <25000 and >25000
    T2(altitude) = -56.46 + 273.1 #units: K
    P2(altitude) = 22.65 * exp(1.73 - 0.000157 * altitude) #units: kPa
    T3(altitude) = -131.21 + 0.00299 * altitude + 273.1 #units: K
    P3(altitude) = 2.488 / ((T3(altitude) / 216.6)^11.388) #units: kPa

    return FLOWMath.sigmoid_blend(T2(altitude), T3(altitude), altitude, 25000, hardness),
    FLOWMath.sigmoid_blend(P2(altitude), P3(altitude), altitude, 25000, hardness)
end

"""
    ideal_gas_rho(P, T)

Ideal gas law for calculating density in SI units.
"""
function ideal_gas_rho(P, T)
    return @. P / (0.2869 * T)
end

"""
    sutherlands_law(
        static_temperure, mu_sea_level=1.789e-5, T_sea_level=288.15, S=110.4
    )

Sutherland's law in SI units for calculating air viscosity relative to sea level.
"""
function sutherlands_law(
    static_temperure, mu_sea_level=1.789e-5, T_sea_level=288.15, S=110.4
)
    return mu_sea_level * (static_temperure / T_sea_level)^(3.0 / 2.0) * (T_sea_level + S) /
           (static_temperure + S)
end

"""
    standard_atmosphere(altitude; hardness=25)
    standard_atmosphere(::Imperial, altitude; hardness=25)

Smoothed fits to the Standard Atmosphere model.

Assumes calorically imperfect gas.

# Arguments
- `altitude::Float` : Altitude in meters for SI units, or feet for Imperial units

# Keyword Arguments:
- `hardness::float` : hardness factor for sigmoid blend

# Returns
- `static_temperature::Float` : Static temperature
- `static_pressure::Float` : Static pressure
- `static_density::Float` : Static density
- `static_dynamic_viscosity::Float` : Static dynamic Viscosity
"""
function standard_atmosphere(altitude; hardness=25)

    #Get temperature (T) and pressure (P) from table fits
    if altitude < (11000 + 25000) / 2.0
        T, P = sa1(altitude; hardness=hardness)
    else
        T, P = sa2(altitude; hardness=hardness)
    end

    # return T, P, rho, mu
    return T, P * 1000, ideal_gas_rho(P, T), sutherlands_law(T)
end

function standard_atmosphere(imperial_units, altitude; hardness=25)

    # convert from feet to meters
    altitude *= 0.3048

    #Get temperature (T) and pressure (P) from table fits in SI units
    if altitude < (11000 + 25000) / 2.0
        T, P = sa1(altitude; hardness=hardness)
    else
        T, P = sa2(altitude; hardness=hardness)
    end

    # - Convert to Imperial Units - #

    # convert from kg/m^3 to slugs/ft^3
    rho = ideal_gas_rho(P, T) * 0.00194032

    # convert from Pa-s to slugs/ft-s
    mu = sutherlands_law(T) * 0.0208854342

    # convert from Kelvin to Fahrenheit
    T -= 273.15
    T *= 9.0 / 5.0
    T += 32.0

    # convert from kilo pascals to slugs/ft^2
    P *= 20.885434273039

    return T, P, rho, mu
end

"""
    speed_of_sound(static_pressure, static_density; gamma=1.4)

Speed of sound from isentropic relations
"""
function speed_of_sound(static_pressure, static_density; gamma=1.4)
    return sqrt(gamma * static_pressure / static_density)
end

"""
    calculate_mach(edge_velocity, speed_of_sound)

Mach number from velocity and speed of sound
"""
function calculate_mach(edge_velocity, speed_of_sound)
    return edge_velocity / speed_of_sound
end

"""
    total_pressure(static_pressure, M; gamma=1.4)

Total pressure from isentropic relations
"""
function total_pressure(static_pressure, M; gamma=1.4)
    return static_pressure * (1.0 + (gamma - 1.0) * M^2 / 2.0)^(gamma / (gamma - 1.0))
end

"""
    total_temperature(static_temperature, M; gamma=1.4)

Total temperature from isentropic relations
"""
function total_temperature(static_temperature, M; gamma=1.4)
    return static_temperature * (1.0 + (gamma - 1.0) * M^2 / 2.0)
end

"""
    static_pressure(total_pressure, M; gamma=1.4)

Static pressure from isentropic relations
"""
function static_pressure(total_pressure, M; gamma=1.4)
    return total_pressure / ((1.0 + (gamma - 1.0) * M^2 / 2.0)^(gamma / (gamma - 1.0)))
end

"""
    static_temperature(total_temperature, M; gamma=1.4)

Static temperature from isentropic relations
"""
function static_temperature(total_temperature, M; gamma=1.4)
    return total_temperature / (1.0 + (gamma - 1.0) * M^2 / 2.0)
end

"""
    static_density(static_pressure, speed_of_sound; gamma=1.4)

Static density from isentropic relations
"""
function static_density(static_pressure, speed_of_sound; gamma=1.4)
    return gamma * static_pressure / speed_of_sound^2
end

"""
    convert_temperature_to_kelvin(::Units, T)

Convert from Fahrenheit to Kelvin or Return temperature if already in SI units
"""
function convert_temperature_to_kelvin(::SI, T)
    return T
end

function convert_temperature_to_kelvin(::Imperial, T)

    # convert to celsius
    T -= 32.0
    T *= 5.0 / 9.0

    # convert to kelvin
    T += 273.15

    return T
end

"""
    convert_viscosity(::SI, mu)

Convert viscosity from Imperial units to SI or Return input if already SI units.
"""
function convert_viscosity(::SI, mu)
    return mu
end

function convert_viscosity(::Imperial, mu)
    return mu * 0.0208854342
end

