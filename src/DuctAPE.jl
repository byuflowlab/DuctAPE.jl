module DuctAPE

## -- DEPENDENCIES -- ##

# using FLOWFoil
# const ff = FLOWFoil #rename FLOWFoil for convenience

using FLOWMath
const fm = FLOWMath #rename FLOWMath for convenience

using CCBlade
const ccb = CCBlade #rename ccblade for convenience

using NLsolve
using ImplicitAD

using SpecialFunctions

using LinearAlgebra: factorize, mul!, lu!, ldiv!, issuccess, NoPivot

using QuadGK

## -- INCLUDES -- ##

include("types.jl")

# Utility Functions
include("utils.jl")

# Airfoil Parameterizations
include("airfoil_parameters/naca_65series.jl")

# Cascade Functions
include("cascade.jl")

# Body Geometry Functions
include("body_geometry.jl")
include("panel.jl")

# Rotor Geometry Functions
include("rotor_geometry.jl")

# Wake Geometry Functions
include("wake_geometry.jl")

# Additional Meshing Functions
# include("mesh.jl")

# Additional Influence Coefficient Functions
# include("coefficient_matrix.jl") #TODO: delete or move this?
include("integrals.jl")
include("velocities.jl")
include("influence_coefficient_matrices.jl")
include("induced_velocity_matrices.jl")

# Rotor Aerodynamic Functions
include("rotor_aerodynamics.jl")
include("airfoil_corrections.jl")

# Wake Aerodynamic Functions
include("wake_aerodynamics.jl")

# Body Aerodynamic Functions
include("body_aerodynamics.jl")

# Kutta Condition Residual
include("pressure_residual.jl")

# Pre-solve initializations
include("initialize.jl")

# Solver
include("solve.jl")

# Post Process
include("post_process.jl")

# -- ROTOR ONLY -- ##
include("solve_rotor_only.jl")
include("initialize_rotor_only.jl")
include("post_process_rotor.jl")

# -- Debugging -- ##
include("initialize_manual.jl")
end
