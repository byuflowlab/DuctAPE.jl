#=
AERO
!  #sections
     1
!   Xisection
   0.0000
!       A0deg        dCLdA        CLmax         CLmin
   0.0000       6.2800       1.5000      -1.0000
!  dCLdAstall     dCLstall      Cmconst         Mcrit
  0.50000      0.20000       0.0000      0.70000
!       CDmin      CLCDmin     dCDdCL^2
  0.12000E-01  0.10000      0.50000E-02
!       REref        REexp
  0.20000E+06  0.35000
ENDAERO

function generate_airfoil_data(lift_params, drag_params, Re; filename=nothing)
    # extract parameters
    (; clmax, clmin, alpha0, dclda, dclda_stall, blend_hardness) = lift_params
    (; cdo, clo, b, Re_ref, f) = drag_params
=#

project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

include(project_dir * "/plots_default.jl")

include(project_dir * "/examples/dfdc_comp/simple_airfoil.jl")

clmax = 1.5
clmin = -1.0
alpha0 = 0.0
dclda = 2.0 * pi
dclda_stall = -0.5
blend_hardness = 5
lift_params = (; clmax, clmin, alpha0, dclda, dclda_stall, blend_hardness)

cdo = 0.012
clo = 0.1
b = 0.005
Re_ref = 2e5
f = 0.35
drag_params = (; cdo, clo, b, Re_ref, f)

generate_airfoil_data(
    lift_params,
    drag_params,
    5e5;
    filename=project_dir * "/examples/dfdc_comp/dfdc_af_test1.dat",
)
