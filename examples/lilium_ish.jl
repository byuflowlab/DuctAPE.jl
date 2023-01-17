#=

Example based on Lilium-like Geometry

=#

using DuctTAPE
const dt = DuctTAPE

# Load in Lilium-like Geometry
include("lilium_parameters.jl")

# Load convenience Parameterizations
include("parameterizations.jl")

#---------------------------------#
#             GEOMETRY            #
#---------------------------------#

# - Get Hub Geometry - #
hx, hr, ctrlpts, knot_points = generate_hub_coordinates(;
    hub_le=hub_le,
    hub_te=hub_te,
    rmax=Rhub,
    nose_angle=pi / 4.0, # pi/4 from trial and error
    tail_angle=pi / 9.0, # pi/9 from trial and error
    cyl_x=[cyl_le; cyl_te],
    degree=3,
    weights=nothing,
    N=80,
    debug=true,
)

# - Get Duct Geometry - #
dxo, dro, dxi, dri = generate_duct_coordinates(;
    r_le=0.01, # try to match public images
    le_x=duct_le[1],
    te_x=duct_te[1],
    inlet_radius=duct_le[2],
    outlet_radius=duct_te[2],
    wedge_angle=2.0, # try to match public images
    te_camber_angle=15.0, # try to match public images
    ctrlpt_xpos=1.0 / 6.0, # try to match public images
    N=80,
)

# - Assemble Coordinates - #
duct_coordinates = [[reverse(dxi); dxo[2:end]] [reverse(dri); dro[2:end]]]
hub_coordinates = [hx hr]

body_geometry, body_panels = dt.generate_body_geometry(duct_coordinates, hub_coordinates)

# - Get Rotor Geometry - #

#choose number of blade elements
nbe = 10

# - Chord - #
# note: just use linear chord distribution
# TODO: probaby need to do some reverse engineering with chord and twist to get the values in the parameters file, since the figures they are based off of likely show chord and twist together.
chords = range(croot_rotor, croot_stator; length=nbe)

# - Twist - #
# note: need to make a guess, use degrees for input
twists = range(40.0, 10.0; length=nbe)

# - Non-dimensional Radial Locations - #
radial_positions = range(0.0, 1.0; length=nbe)

# - TODO: add airfoil data - #
airfoils = nothing

blade_elements, rotor_panels = dt.generate_blade_elements(
    rotor_c4_pos, radial_positions, chords, twists, airfoils, body_geometry
)

#---------------------------------#
#             Meshes              #
#---------------------------------#

mesh_bodies_to_rotor = dt.generate_one_way_mesh(body_panels, rotor_panels)

#---------------------------------#
#       Coefficient Matrices      #
#---------------------------------#

A_bodies_to_rotor = dt.assemble_one_way_source_matrix(
    mesh_bodies_to_rotor, body_panels, rotor_panels
)
