#=
Check variations in wake contraction and expansion to make sure outputs are as expected
=#

project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

include(project_dir * "/plots_default.jl")

using DuctAPE
const dt = DuctAPE

# CCBlade used for it's airfoils function objects here.
using CCBlade
const ccb = CCBlade
include(project_dir * "/dev_debug_archive/rotor_only/run_ccblade.jl")

using FLOWMath
const fm = FLOWMath

#---------------------------------#
#             Geometry            #
#---------------------------------#

# Blade Tip Radius, in meters
Rtip = 10 / 2.0 * 0.0254  # inches to meters

# Blade Hub radius, in meters
Rhub = 0.10 * Rtip

# number of blades
B = 2

# Blade section non-dimensional radial positions, chords lengths, and local twists angles in degrees
propgeom = [
    0.15 0.130 32.76
    0.20 0.149 37.19
    0.25 0.173 33.54
    0.30 0.189 29.25
    0.35 0.197 25.64
    0.40 0.201 22.54
    0.45 0.200 20.27
    0.50 0.194 18.46
    0.55 0.186 17.05
    0.60 0.174 15.97
    0.65 0.160 14.87
    0.70 0.145 14.09
    0.75 0.128 13.39
    0.80 0.112 12.84
    0.85 0.096 12.25
    0.90 0.081 11.37
    0.95 0.061 10.19
    # 1.00 0.041 8.99
]

# extract non-dimensional radial positions
r = propgeom[:, 1]
# Dimensionalize chords
chords = propgeom[:, 2] * Rtip
# convert twists to radians
twists = propgeom[:, 3] * pi / 180

airfoils = fill(ccb.AlphaAF("test/data/rotorzloc_af_test.dat"), length(r))

nwake_sheets = 15

rotorzloc = 0.25

wake_length = 2.0
wake_sections = [0.5; 1.0; 1.5]
wake_x_refine = 0.01

#---------------------------------#
#       Operating Conditions      #
#---------------------------------#

# rotor rotation rate in rad/s
Omega = 5400 * pi / 30  # convert from RPM to rad/s

# freestream conditions
rho = 1.225 #kg/m^3
mu = 1.81e-5 # kg/(m⋅s)
asound = 341.0 #m/s

#---------------------------------#
#          Define Inputs          #
#---------------------------------#

# Rotor Parameters
rotor_parameters = [(;
    rotorzloc, nwake_sheets, r, chords, twists, airfoils, Rtip, Rhub, B, Omega
)]

# Freestream Parameters
Vinf = 5.0
freestream = (; rho, mu, asound, Vinf)

#---------------------------------#
#       Initilaize Inputs         #
#---------------------------------#

initial_states, default_inputs = dt.initialize_rotor_states(
    rotor_parameters, freestream; wake_length=wake_length, wake_x_refine=wake_x_refine
)
nominal_strengths = dt.solve_rotor_only(initial_states, default_inputs)

nominal_Gamr, nominal_gamw, nominal_sigr = dt.extract_rotor_states(nominal_strengths, default_inputs)

#---------------------------------#
#    overwrite wake panels and    #
#      associated velocities      #
#---------------------------------#

## -- define grid points
# keep same x points
x_grid_points = default_inputs.x_grid_points

# update radial points
rscales = [0.5; 1.5; 1.0; 1.0]
xids = [
        findfirst(x -> x >= rotorzloc + wake_sections[i] * 2.0 * Rtip, x_grid_points[:,1]) for
    i in 1:length(wake_sections)
]
xidx = [1; xids; length(x_grid_points[:,1])]

r_wake = similar(x_grid_points[:,1]) .= Rtip
for i in 1:(length(xidx) - 1)
    r_wake[xidx[i]:xidx[i + 1]] = range(
        r_wake[xidx[i]], r_wake[xidx[i]] * rscales[i]; length=xidx[i + 1] - xidx[i] + 1
    )
end


r_grid_points = similar(x_grid_points) .= 0.0
for i in 1:length(r_wake)
    r_grid_points[i, :] = range(Rhub, r_wake[i]; length=length(r_grid_points[1,:]))
end

wake_vortex_panels = dt.generate_wake_panels(x_grid_points, r_grid_points)

##### ----- Wake to Rotor ----- #####
wake_to_rotor_mesh = [
    dt.generate_one_way_mesh(wake_vortex_panels[j], default_inputs.rotor_source_panels[i])
    for i in 1:length(default_inputs.rotor_source_panels), j in 1:length(wake_vortex_panels)
]

##### ----- Wake to Rotor ----- #####
A_wake_to_rotor = [
    dt.assemble_induced_velocity_matrices(
        wake_to_rotor_mesh[i, j],
        wake_vortex_panels[j],
        default_inputs.rotor_source_panels[i],
    ) for i in 1:length(default_inputs.rotor_source_panels),
    j in 1:length(wake_vortex_panels)
]

# - Axial - #
vx_rw = [
    A_wake_to_rotor[i, j][1] for i in 1:length(default_inputs.rotor_source_panels),
    j in 1:length(wake_vortex_panels)
]

# - Radial - #
vr_rw = [
    A_wake_to_rotor[i, j][2] for i in 1:length(default_inputs.rotor_source_panels),
    j in 1:length(wake_vortex_panels)
]

#actually overwrite the velocity and panel fields in the inputs.
inputs = (;
    default_inputs..., wake_vortex_panels=wake_vortex_panels, vx_rw=vx_rw, vr_rw=vr_rw
)

# SOLVE
strengths = dt.solve_rotor_only(initial_states, inputs)

Gamr, gamw, sigr = dt.extract_rotor_states(strengths, inputs)

## -- run ccblade -- ##
ccbouts = run_ccblade(Vinf; airfoil="test/data/rotorzloc_af_test.dat")

#---------------------------------#
#             PLOTS               #
#---------------------------------#

# - Plot Geometry - #
pgeom = plot(; aspectratio=1, xlabel="x", ylabel="r")
# Plot rotor position
plot!(
    pgeom,
    rotorzloc * ones(length(inputs.rotor_panel_edges)),
    inputs.rotor_panel_edges;
    label="",
)

# plot wake geometry
plot!(pgeom, x_grid_points, r_grid_points; color=mycolors[2], linewidth=0.5, label="")

savefig(project_dir * "/dev_debug_archive/rotor_only/wakeshapegeom.pdf")

# - Plot Circulation - #
pc = plot(; xlabel="Circulation", ylabel="r/Rtip")
# plot ccblade solution
plot!(pc, ccbouts.circ, r; label="ccblade")
# plot DuctAPE solution
plot!(pc, Gamr, inputs.rotor_panel_centers ./ Rtip; label="DuctAPE")
plot!(pc, nominal_Gamr, default_inputs.rotor_panel_centers ./ Rtip; label="DuctAPE straight wake")

savefig(project_dir * "/dev_debug_archive/rotor_only/wakeshapecirc.pdf")

plot(pgeom, pc; layout=(2, 1))

savefig(project_dir * "/dev_debug_archive/rotor_only/wakeshapecheck.pdf")
