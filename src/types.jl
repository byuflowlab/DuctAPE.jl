#=
Custom Composite Type Definitions
=#

"""
    Freestream{TVI,TVR,TF}

**Fields:**
 - `vinf::Float` : Freestream velocity
 - `rho::Float` : Air density value, default = 1.225 kg/m^3
 - `mu::Float` : Air dynamic viscosity value, default = 1.81e-5 Pa-s
 - `asound::Float` : Speed of sound value, default = 341.0 m/s
"""
struct Freestream{TF,TC}
    vinf::TF
    rho::TC # defaults to 1.225 kg/m^3
    mu::TC # defaults to 1.81e-5 Pa-s
    asound::TC # defaults to 341.0 m/s
end

function Freestream(vinf)
    return Freestream(vinf, 1.225, 1.81e-5, 341.0)
end
