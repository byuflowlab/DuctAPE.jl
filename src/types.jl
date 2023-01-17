#=
Custom Composite Type Definitions
=#

"""
    OperatingConditions{TVI,TVR,TF}

**Fields:**
 - `vinf::Float` : Freestream velocity
 - `rho::Float` : Air density value
 - `mu::Float` : Air dynamic viscosity value
 - `asound::Float` : Speed of sound value
"""
struct OperatingConditions{TF,TC}
    vinf::TF
    rho::TC
    mu::TC
    asound::TC
    RPM::TF
end
