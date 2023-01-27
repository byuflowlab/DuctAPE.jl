module DuctTAPE

## -- DEPENDENCIES -- ##

using FLOWFoil
const ff = FLOWFoil #rename FLOWFoil for convenience

using FLOWMath
const fm = FLOWMath #rename FLOWMath for convenience

using CCBlade
const ccb = CCBlade #rename ccblade for convenience

using NLsolve
using ImplicitAD

## -- INCLUDES -- ##

include("utils.jl")

# Body Geometry Functions
include("body_geometry.jl")

# Rotor Geometry Functions
include("rotor_geometry.jl")

# Wake Geometry Functions
include("wake_geometry.jl")

# Additional Meshing Functions
include("mesh.jl")

# Additional Influence Coefficient Functions
include("coefficient_matrix.jl")

# Rotor Aerodynamic Functions
include("rotor_aerodynamics.jl")

# Pre-solve initializations
include("initialize.jl")

end
