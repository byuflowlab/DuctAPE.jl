
using DuctTAPE
const dt = DuctTAPE
using FLOWFoil
const ff = FLOWFoil
using ImplicitAD
const ia = ImplicitAD
using CCBlade
const ccb = CCBlade

# Load in Lilium-like Geometry
include("lilium_parameters.jl")

# Load convenience Parameterizations
include("parameterizations.jl")

# - My guestimates for other parameters needed - #
hub_nose_angle = pi / 4.0
hub_tail_angle = pi / 9.0
duct_wedge_angle = 2.0
duct_le_radius = 0.01
duct_te_camber = 15.0
duct_ctrlpt_x = 1.0 / 6.0

omega_rotor = 2500.0  * pi / 30.0 #2500 rpm to rad/s

x = [
    hub_le
    hub_te
    Rhub
    hub_nose_angle
    hub_tail_angle
    cyl_le
    cyl_te
    duct_le[1]
    duct_le[2]
    duct_te[1]
    duct_te[2]
    duct_wedge_angle
    duct_le_radius
    duct_te_camber
    duct_ctrlpt_x
    rotor_chord_guess
    rotor_root_twist_guess
    rotor_tip_twist_guess
    rotor_c4_pos
    stator_c4_pos
    stator_root_chord
    stator_tip_chord
    stator_twist_guess
    omega_rotor
]

hub_le = x[1]
hub_te = x[2]
Rhub = x[3]
hub_nose_angle = x[4]
hub_tail_angle = x[5]
cyl_le = x[6]
cyl_te = x[7]
duct_lex = x[8]
duct_ler = x[9]
duct_tex = x[10]
duct_ter = x[11]
duct_wedge_angle = x[12]
duct_le_radius = x[13]
duct_te_camber = x[14]
duct_ctrlpt_x = x[15]
rotor_chord_guess = x[16]
rotor_root_twist_guess = x[17]
rotor_tip_twist_guess = x[18]
rotor_c4_pos = x[19]
stator_c4_pos = x[20]
stator_root_chord = x[21]
stator_tip_chord = x[22]
stator_twist_guess = x[23]
omega_rotor = x[24]

method_body = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false, true])
method_vortex = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true])
method_source = ff.AxisymmetricProblem(Source(Constant()), Neumann(), [true])

problem_body = ff.define_problem(
    method_body, (duct_coordinates, hub_coordinates), 0.0, -1.0, -1.0
)

#---------------------------------#
#             GEOMETRY            #
#---------------------------------#

# - Get Hub Geometry - #
hx, hr = generate_hub_coordinates(;
    hub_le=hub_le,
    hub_te=hub_te,
    rmax=Rhub,
    nose_angle=hub_nose_angle, # pi/4 from trial and error
    tail_angle=hub_tail_angle, # pi/9 from trial and error
    cyl_x=[cyl_le; cyl_te],
    N=80,
    debug=true,
)

# - Get Duct Geometry - #
dxo, dro, dxi, dri = generate_duct_coordinates(;
    r_le=duct_le_radius, # try to match public images
    le_x=duct_lex,
    te_x=duct_tex,
    inlet_radius=duct_ler,
    outlet_radius=duct_ter,
    wedge_angle=duct_wedge_angle, # try to match public images
    te_camber_angle=duct_te_camber, # try to match public images
    ctrlpt_xpos=duct_ctrlpt_x, # try to match public images
    N=80,
)

# - Assemble Coordinates - #
duct_coordinates = [[reverse(dxi); dxo[2:end]] [reverse(dri); dro[2:end]]]
hub_coordinates = [hx hr]

body_geometry, body_panels = dt.generate_body_geometry(duct_coordinates, hub_coordinates)

# - Get Rotor Geometry - #

#choose number of blade elements
num_rotors = 1
nbe = 5
nbe_fine = 10

# - Chord - #
chords = range(rotor_chord_guess, rotor_chord_guess; length=nbe)

# - Twist - #
# note: need to make a guess, use degrees for input
twists = range(rotor_root_twist_guess, rotor_tip_twist_guess; length=nbe)

# - Non-dimensional Radial Locations - #
rotor_radial_positions = range(0.0, 1.0; length=nbe)

# - Airfoil Data - #
af = ccb.AlphaAF("test/data/naca_4412_extrapolated_rotated_APCshifted.dat")
# airfoils = [nothing for i in 1:nbe]
airfoils = fill(af, nbe)

# - Number of blades - #
B = 27

rotor_parameters = [(
    rotor_x_position=rotor_c4_pos,
    radial_positions=rotor_radial_positions,
    chords=chords,
    twists=twists,
    airfoils=airfoils,
    num_radial_stations=nbe,
    num_blades=B,
    Omega=omega_rotor,
)]

blade_elements = dt.initialize_blade_elements(
    rotor_parameters[1], body_geometry, num_rotors
)

######################################################################
#                                                                    #
#                           PLOT WAKE GRID                           #
#                                                                    #
######################################################################
using Plots

# - Initialize Grid Points - #
x_grid_points, r_grid_points, nx, nr, rotoridxs = dt.initialize_grid_points(
    body_geometry,
    rotor_parameters[1].rotor_x_position,
    blade_elements[1].radial_positions;
    wake_length=1.0,
    debug=false,
)

plot(; aspectratio=:equal)
plot!(dxi, dri; color=:black, linestyle=:dash, label="body")
plot!(hx, hr; color=:black, linestyle=:dash, label="")
plot!(
    x_grid_points,
    r_grid_points;
    seriestype=:scatter,
    markercolor=:blue,
    markerstrokecolor=:blue,
    markersize=2,
    label="",
)
savefig("initial_grid_points.pdf")

# - Relax Grid Points - #
dt.relax_grid!(x_grid_points, r_grid_points, nx, nr)
plot!(
    x_grid_points,
    r_grid_points;
    seriestype=:scatter,
    markercolor=:red,
    markerstrokecolor=:red,
    markersize=1,
    label="",
)
savefig("relaxed_grid_points.pdf")

# - Generate Panels - #
# TODO: just put the wake panels in the middle here rather than elsewhere
wake_panels = dt.generate_wake_panels(
    (x_grid_points[:, 1:(end - 1)] .+ x_grid_points[:, 2:end]) / 2.0,
    (r_grid_points[:, 1:(end - 1)] .+ r_grid_points[:, 2:end]) / 2.0,
    nr - 1;
    method=method_vortex,
)

for i in 1:length(wake_panels)
    plot!(
        wake_panels[i].panel_center[:, 1],
        wake_panels[i].panel_center[:, 2];
        seriestype=:scatter,
        markercolor=:green,
        markerstrokecolor=:green,
        markersize=1,
        label="",
    )
end
savefig("wake_panel_centers.pdf")
