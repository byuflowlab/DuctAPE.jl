module DuctAPE

#---------------------------------#
#          DEPENDENCIES           #
#---------------------------------#

# - augmented CCBlade implementation (cascade compatible CCBlade) - #
# NOTE: used for robust state initialziation
include("C4Blade/C4Blade.jl")
const c4b = C4Blade

# - Packages for Calculating Unit Induced Velocities - #
# For Kernels
using SpecialFunctions # required for elliptic integrals

# For Integration
using FastGaussQuadrature
using QuadGK

# - Packages for Code Efficiency - #
using StaticArrays # used in miscellaneous places for code efficiency
using PreallocationTools # caches

# - Packages for Solves - #
using ImplicitAD # used for all solves
using LinearAlgebra # linear solve and LU decomposition

# General Nonlinear solves
using SimpleNonlinearSolve # SimpleDFSane and SimpleNewtonRaphson
# using PolyesterForwardDiff

# Quasi-Newton
using SIAMFANLEquations
using MINPACK

# Fixed-Point Iteration Solvers
using SpeedMapping
using FixedPoint

# For using NLsolve
using NLsolve #Includes Anderson Solver
using LineSearches # used in newton solver
using ForwardDiff # used for jacobian for newton solver

# - Utility Packages - #
using FLOWMath # used for various items, mostly interpolation
using Printf # used when verbose option is selected

#---------------------------------#
#             EXPORTS             #
#---------------------------------#

# - Nested Modules - #
export c4b

# - Types - #

# - Analyses - #

# - Preprocess - #

# - Pocess - #

# - Postprocess - #

#---------------------------------#
#            INCLUDES             #
#---------------------------------#

##### ----- UTILITIES ----- #####
# general utility functions
include("utilities/misc.jl")
include("utilities/options.jl")
include("utilities/inputs.jl")
include("utilities/update_propulsor.jl")
include("utilities/package_states.jl")
include("utilities/bookkeeping.jl")
include("utilities/caching/caches.jl")
include("utilities/caching/allocate_caches.jl")
include("utilities/caching/reshape_caches.jl")
include("utilities/caching/integration_caches.jl")

# Airfoil utility functions
include("utilities/airfoils/airfoil_utilities.jl")
include("utilities/airfoils/naca_65series.jl")

##### ----- Analysis ----- #####
include("analysis/setup.jl")
include("analysis/analyses.jl")

##### ----- PREPROCESS ----- #####
# Pre-solve initializations
include("preprocess/preprocess.jl")
include("preprocess/initialize_states.jl")

# Geometry Functions
include("preprocess/geometry/body_geometry.jl")
include("preprocess/geometry/panel_geometry.jl")
include("preprocess/geometry/rotor_geometry.jl")
include("preprocess/geometry/wake_geometry.jl")
include("preprocess/geometry/elliptic_grid_residuals.jl")

# Induced Velocity Functions
include("preprocess/velocities/unit_induced_velocities.jl")
include("preprocess/velocities/induced_velocity_matrices.jl")
include("preprocess/velocities/body_aic.jl")

# Quadrature
include("preprocess/velocities/quadrature/integrands.jl")
include("preprocess/velocities/quadrature/out_of_place_integrals.jl")
include("preprocess/velocities/quadrature/gausslegendre_integrals.jl")
include("preprocess/velocities/quadrature/romberg_integrals.jl")
include("preprocess/velocities/quadrature/gausskronrod_integrals.jl")

##### ----- PROCESS ----- #####
# Solve and Residual Functions
include("process/process.jl")
include("process/solve.jl")
include("process/residuals/CSORresidual.jl")
include("process/residuals/systemresidual.jl")

# Aerodynamics Functions
include("process/aerodynamics/body_aerodynamics.jl")
include("process/aerodynamics/rotor_aerodynamics.jl")
include("process/aerodynamics/wake_aerodynamics.jl")

##### ----- POST-PROCESSING ----- #####
include("postprocess/postprocess.jl")
include("postprocess/velocities.jl")
include("postprocess/pressures.jl")
include("postprocess/rotor_performance.jl")
include("postprocess/utils.jl")

##### ----- DEBUGGING ----- #####
include("../test/test_utils.jl")
end
