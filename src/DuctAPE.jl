module DuctAPE

#---------------------------------#
#          DEPENDENCIES           #
#---------------------------------#

using FLOWMath # used for various items, mostly interpolation

include("C4Blade/C4Blade.jl") # augmented CCBlade implementation (cascade compatible CCBlade)
const c4b = C4Blade

using SpecialFunctions # required for elliptic integrals
using QuadGK # required for integration of linear panels
using StaticArrays # used in miscellaneous places for code efficiency

using LinearAlgebra

using PreallocationTools

# new solve required pacakges
using NLsolve #for newton solver
using LineSearches
using ImplicitAD
using ForwardDiff

using Printf # used when verbose option is selected

#---------------------------------#
#             EXPORTS             #
#---------------------------------#
export c4b

#---------------------------------#
#            INCLUDES             #
#---------------------------------#
##### ----- TYPES ----- #####
include("types.jl")

##### ----- UTILITIES ----- #####
# general utility functions
# TODO: organize utility
include("utilities/utils.jl")
include("solve/io.jl") #TODO: move to utilities
include("solve/cache.jl")
include("solve/bookkeeping.jl")

##### ----- Analysis ----- #####
include("analysis/analyses.jl")

##### ----- PRECOMPUTATION ----- #####
# Pre-solve initializations
include("precomputation/initialize.jl")
include("precomputation/initialize_iad.jl")

# Body Geometry Functions
include("precomputation/body_geometry.jl")
include("precomputation/panel.jl")

# Rotor Geometry Functions
include("precomputation/rotor_geometry.jl")

# Wake Geometry Functions
include("precomputation/wake_geometry.jl")
include("precomputation/wake_geometry_residual.jl")

# Aero Influence Matrices
include("precomputation/integrals.jl")
include("precomputation/velocities.jl")
include("precomputation/body_aic.jl")
include("precomputation/induced_velocity_matrices.jl")

##### ----- SOLVER ----- #####

include("solve/solve.jl")
include("solve/solve_iad.jl")

# Rotor Aerodynamic Functions
include("solve/rotor_aerodynamics.jl")

# Wake Aerodynamic Functions
include("solve/wake_aerodynamics.jl")

# Body Aerodynamic Functions
include("solve/body_aerodynamics.jl")

##### ----- POST-PROCESSING ----- #####

# include("postprocess/post_process.jl")
include("postprocess/utils.jl")
include("postprocess/post_process_iad.jl")

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
include("precomputation/initialize_rotor_only.jl")
include("solve/solve_rotor_only.jl")
include("postprocess/post_process_rotor.jl")

##### ----- DEBUGGING ----- #####
include("../test/test_utils.jl")

end
