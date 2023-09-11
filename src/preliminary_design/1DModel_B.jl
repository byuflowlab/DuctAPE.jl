#=
This version takes a given open rotor geometry and operating point and outputs the required exit velocity to maintain a similar performance.
=#

"""

**Arguments:**
- `Rtip::Float` : Rotor tip radius
- `RPM::Float` : rotations per minute
- `J::Float` : Advance ratio
- `ct::Float` : thrust coefficient
- `Vinf::Float` : freestream velocity for ducted operation

**Keyword Arguments:**
- `rho::Float` : air density, default = 1.225 kg/m^3

**Returns:**
- `exit_area::Float` : exit area
- `debug::NamedTuple` : named tuple containing intermediate calculation values
"""
function size_exit(Rtip, RPM, J, ct, Vinf; rho=1.225)
    vstar = calculate_vstar(Rtip, RPM, J)

    thrust = calculate_thrust(ct, RPM, Rtip; rho=1.225)

    dP = calculate_pressure_rise(thrust, Rtip)

    vjetstar = calculate_jet_velocity(vstar, dP; rho=1.225)

    Vface = caluclate_fan_face_velocity(vstar, vjetstar)

    Vjet = calculate_jet_velocity(Vinf, dP; rho=1.225)

    exit_area = calculate_jet_area(Rtip, Vface, Vjet)

    debug = (; vstar, thrust, dP, vjetstar, Vface, Vjet)

    return exit_area, debug
end

"""
"""
function calculate_vstar(Rtip, RPM, J)
    return (2.0 * Rtip) * (RPM / 60) * J
end

"""
"""
function calculate_thrust(ct, RPM, Rtip; rho=1.225)
    return ct * rho * (RPM / 60)^2 * (2.0 * Rtip)^4
end

"""
"""
function calculate_pressure_rise(thrust, Rtip)
    return thrust / (pi * Rtip^2)
end

"""
"""
function calculate_jet_velocity(Vinf, dP; rho=1.225)
    return sqrt(Vinf^2 + 2.0 * dP / rho)
end

"""
"""
function caluclate_fan_face_velocity(vstar, vjet)
    return 0.5 * (vstar + vjet)
end

"""
"""
function calculate_jet_area(Rtip, Vface, Vjet)
    return pi * Rtip^2 * Vface / Vjet
end

