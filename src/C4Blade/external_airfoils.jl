"""
    ExternalAirfoil

# Fields
- `parameters::Vector` : parameters for geoemtry used in `compute_aero`.
- `compute_aero::Function` : funciton of the form `f(parameters, angle_of_attack, reynolds_number, mach_number, local_solidity, stagger_angle)`; returns cl and cd
"""
@kwdef struct ExternalAirfoil{TP,TH}
    parameters::TP
    compute_aero::TH
end

"""
    function external_eval(
        inflow_magnitude,
        local_reynolds,
        local_solidity,
        local_stagger,
        alpha,
        airfoil_object,
        asound;
        verbose=false,
        is_stator=0,
    )

DFDC-like polar function.

# Arguments
- `inflow_magnitude::Float` : Velocity magnitude of inflow
- `local_reynolds::Float` : local Reynolds number
- `local_solidity::Float` : local Solidity
- `local_stagger::Float` : local Stagger angle (radians)
- `alpha::Float` : local angle of attack (radians)
- `airfoil_object::C4Blade.ExternalAirfoil` : ExternalAirfoil object
- `asound::Float` : speed of sound

# Keyword Arguments
- `verbose::Bool=false::` : print verbose statements
- `is_stator::Int=0` : flag to flip lift values (e.g. for stators)

# Returns
- `cl::Float` : lift coefficient corrected for compressibility, solidity and stagger as required.
- `cd::Float` : drag coefficient corrected for compressibility, solidity and stagger as required.
"""
function external_eval(
    Wmag,
    reynolds,
    solidity,
    stagger,
    alpha,
    inner_airfoil,
    asound;
    verbose=verbose,
    is_stator=0.0,
)
    return inner_airfoil.compute_aero(
        inner_airfoil.parameters,
        iszero(is_stator) ? alpha : -alpha,
        reynolds,
        Wmag / asound,
        solidity,
        stagger,
    )
end
