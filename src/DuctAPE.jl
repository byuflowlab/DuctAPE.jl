module DuctAPE

#---------------------------------#
#          DEPENDENCIES           #
#---------------------------------#

using FLOWMath # used for various items, mostly interpolation

include("C4Blade/C4Blade.jl") # augmented CCBlade implementation (cascade compatible CCBlade)
const c4b = C4Blade

using SpecialFunctions # required for elliptic integrals
using QuadGK # integration
using FastGaussQuadrature # integration
using Primes # integration
using StaticArrays # used in miscellaneous places for code efficiency

using LinearAlgebra # linear solve and LU decomposition

using PreallocationTools # caches

# new solve required pacakges
using ImplicitAD # used for all solves
using NLsolve #for newton solver
using LineSearches # used in newton solver
using ForwardDiff # used for jacobian for newton solver

using Printf # used when verbose option is selected

#---------------------------------#
#             EXPORTS             #
#---------------------------------#
export c4b

#---------------------------------#
#            INCLUDES             #
#---------------------------------#

##### ----- UTILITIES ----- #####
# general utility functions
include("utilities/utils.jl")
include("utilities/types.jl")
include("utilities/package_states.jl")
include("utilities/caches.jl")
include("utilities/bookkeeping.jl")

# Airfoil utility functions
include("utilities/airfoils/airfoil_utilities.jl")
include("utilities/airfoils/naca_65series.jl")

##### ----- Analysis ----- #####
include("analysis/analyses.jl")

##### ----- PREPROCESS ----- #####
# Pre-solve initializations
include("preprocess/preprocess.jl")

# Geometry Functions
include("preprocess/geometry/body_geometry.jl")
include("preprocess/geometry/panel_geometry.jl")
include("preprocess/geometry/rotor_geometry.jl")
include("preprocess/geoemtry/wake_geometry.jl")
include("preprocess/geometry/elliptic_grid_residuals.jl")

# Induced Velocity Functions
include("preprocess/velocities/unit_induced_velocities.jl")
include("preprocess/velocities/integrals.jl")
include("preprocess/velocities/induced_velocity_matrices.jl")
include("preprocess/velocities/body_aic.jl")

##### ----- PROCESS ----- #####

include("process/solvers/CSORsolve.jl")
include("process/solvers/solve.jl")

# Aerodynamics Functions
include("process/aerodynamics/body_aerodynamics.jl")
include("process/aerodynamics/rotor_aerodynamics.jl")
include("process/aerodynamics/wake_aerodynamics.jl")

##### ----- POST-PROCESSING ----- #####

# include("postprocess/post_process.jl")
include("postprocess/utils.jl")
include("postprocess/postprocess.jl")

##### ----- DEBUGGING ----- #####
include("../test/test_utils.jl")

end
