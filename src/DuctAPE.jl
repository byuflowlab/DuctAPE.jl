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

# For boundary layer stuff
using Roots

# - Utility Packages - #
using FLOWMath # used for various items, mostly interpolation
using Printf # used when verbose option is selected
using RecipesBase # for plotting

# - Visualization Packages - #
import Colors.RGB
using LaTeXStrings
using ProgressMeter

#---------------------------------#
#             EXPORTS             #
#---------------------------------#

# - Nested Modules - #
export c4b

# - Types - #

# Inputs
export DuctedRotor, Rotor, OperatingPoint, PanelingConstants, ReferenceParameters

# - Preallocations - #
export allocate_prepost_container_cache,
    allocate_solve_parameter_cache, allocate_solve_container_cache

# Options
export Options, set_options, DFDC_options
export IntegrationOptions, GaussLegendre, GaussKronrod, Romberg
export SLORGridSolverOptions, GridSolverOptions
export ChainSolverOptions,
    CompositeSolverOptions,
    NLsolveOptions,
    NonlinearSolveOptions,
    MinpackOptions,
    SIAMFANLEOptions,
    SpeedMappingOptions,
    FixedPointOptions,
    CSORSolverOptions
export BoundaryLayerOptions

# - Preprocess - #
export setup_analysis

# - Analyses - #
export analyze

# - Visualization - #
export generate_plots
export plotGeometry,
    plotDuctGeometry,
    plotBodyGeometry,
    underlayGeometry,
    plotCP,
    plotVtan,
    plotStagnation,
    plotMomentum,
    plotStreamlines,
    staticPlots,
    animatedPlots

#---------------------------------#
#            INCLUDES             #
#---------------------------------#

##### ----- UTILITIES ----- #####
# general utility functions
include("utilities/misc.jl")
include("utilities/inputs.jl")
include("utilities/options.jl")
include("utilities/package_states.jl")
include("utilities/bookkeeping.jl")
include("utilities/caching/caches.jl")
include("utilities/caching/allocate_caches.jl")
include("utilities/caching/reshape_caches.jl")
include("utilities/caching/integration_caches.jl")
include("utilities/caching/elliptic_grid_parameter_packaging.jl")
include("utilities/thermodynamics.jl")
include("utilities/ode_solvers.jl")

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
# Dispatches
include("process/process.jl")
include("process/solve.jl")

# Internal Solvers
include("process/solvers/modCSORsolver.jl")

# Residual Functions
include("process/residuals/CSORresidual.jl")
include("process/residuals/modCSORresidual.jl")
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
include("postprocess/boundary_layer_utils.jl")
include("postprocess/boundary_layer_green.jl")
include("postprocess/boundary_layer_head.jl")
include("postprocess/viscous_drag.jl")

##### ----- VISUALIZATION ----- #####
# include("visualization/plot_recipe_defaults.jl")
include("visualization/plot_recipes.jl")
include("visualization/calculate_streamlines.jl")
include("visualization/convenience_plots.jl")

##### ----- DEBUGGING ----- #####
include("../test/test_utils.jl")
end
