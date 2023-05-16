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
    Vinf::TF
    rho::TC # defaults to 1.225 kg/m^3
    mu::TC # defaults to 1.81e-5 Pa-s
    asound::TC # defaults to 341.0 m/s
end

function Freestream(Vinf)
    return Freestream(Vinf, 1.225, 1.81e-5, 341.0)
end

"""
- `alpha0::Float` : zero lift angle of attack
- `clmax::Float` : maximum cl
- `clmin::Float` : minimum cl
- `dclda::Float` : lift curve slope
- `dclda_stall::Float` :  lift curve slope post-stall
- `dcl_stall::Float` : TODO: explain this
- `cdmin::Float` : minimum cd
- `cldmin::Float` : cl and cdmin
- `dcdcl2::Float` : quadratic curve factor for cd curve
- `cmcon::Float` : pitching moment constant
- `Re_ref::Float` : reference Reynolds number
- `Re_exp::Float` : Reynolds number exponent
- `mcrit::Float` : critical Mach number
"""
struct DFDCairfoil{TF}
    alpha0::TF
    clmax::TF
    clmin::TF
    dclda::TF
    dclda_stall::TF
    dcl_stall::TF
    cdmin::TF
    cldmin::TF
    dcdcl2::TF
    cmcon::TF
    Re_ref::TF
    Re_exp::TF
    mcrit::TF
end
