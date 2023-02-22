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

using LineSearches: BackTracking

using SpecialFunctions

## -- INCLUDES -- ##

include("utils.jl")

include("types.jl")

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

# Wake Aerodynamic Functions
include("wake_aerodynamics.jl")

# Body Aerodynamic Functions
include("body_aerodynamics.jl")

# Pre-solve initializations
include("initialize.jl")

# Solver
include("solve.jl")

# Post Processing
include("post_process.jl")

end
