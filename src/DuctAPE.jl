module DuctAPE

#---------------------------------#
#          DEPENDENCIES           #
#---------------------------------#

using FLOWMath # used for various items, mostly interpolation
const fm = FLOWMath #rename FLOWMath for convenience

using CCBlade # used for airfoil types
const ccb = CCBlade #rename ccblade for convenience

using SpecialFunctions # required for elliptic integrals
using QuadGK # required for integration of linear panels
using StaticArrays

using LinearAlgebra: factorize, mul!, lu!, ldiv!, issuccess, NoPivot # used in linear system assembly and solve

using NLsolve
using ImplicitAD

#---------------------------------#
#            INCLUDES             #
#---------------------------------#

##### ----- UTILITIES ----- #####
include("utilities/types.jl")

# general utility functions
include("utilities/utils.jl")

# Cascade Functions
include("utilities/cascade.jl")
include("utilities/airfoil_corrections.jl")

##### ----- PRECOMPUTATION ----- #####
# Pre-solve initializations
include("precomputation/initialize.jl")

# Body Geometry Functions
include("precomputation/body_geometry.jl")
include("precomputation/panel.jl")

# Rotor Geometry Functions
include("precomputation/rotor_geometry.jl")

# Wake Geometry Functions
include("precomputation/wake_geometry.jl")

# Aero Influence Matrices
include("precomputation/integrals.jl")
include("precomputation/velocities.jl")
include("precomputation/body_aic.jl")
include("precomputation/induced_velocity_matrices.jl")

##### ----- SOLVER ----- #####

include("solve/solve.jl")

# Rotor Aerodynamic Functions
include("solve/rotor_aerodynamics.jl")

# Wake Aerodynamic Functions
include("solve/wake_aerodynamics.jl")

# Body Aerodynamic Functions
include("solve/body_aerodynamics.jl")

##### ----- POST-PROCESSING ----- #####

include("postprocess/post_process.jl")

##### ----- SPECIALTY ----- #####

# Airfoil Parameterizations
include("airfoil_parameters/naca_65series.jl")

# 1D Models
include("preliminary_design/1DModel_A.jl")
include("preliminary_design/1DModel_B.jl")

### TODO Decide what below should be removed and when
# Additional Influence Coefficient Functions
# include("coefficient_matrix.jl") #TODO: delete or move this?
# Kutta Condition Residual
# include("pressure_residual.jl")

# # -- ROTOR ONLY -- ##
# include("solve_rotor_only.jl")
# include("initialize_rotor_only.jl")
# include("post_process_rotor.jl")

# # -- Debugging -- ##
# include("initialize_manual.jl")

end