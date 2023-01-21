#=

Example based on Lilium-like Geometry

=#

using DuctTAPE
const dt = DuctTAPE
using FLOWFoil
const ff = FLOWFoil
using ImplicitAD
const ia = ImplicitAD

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

x0 = [
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
]

"""
wrap everything in a function
"""
function wrapper(x; debug=false)

    #---------------------------------#
    #             Unwrap              #
    #---------------------------------#

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

    body_geometry, body_panels = dt.generate_body_geometry(
        duct_coordinates, hub_coordinates
    )

    # - Get Rotor Geometry - #

    #choose number of blade elements
    nbe = 10

    # - Chord - #
    # note: just use linear chord distribution
    # TODO: probaby need to do some reverse engineering with chord and twist to get the values in the parameters file, since the figures they are based off of likely show chord and twist together.
    chords = range(rotor_chord_guess, rotor_chord_guess; length=nbe)

    # - Twist - #
    # note: need to make a guess, use degrees for input
    twists = range(rotor_root_twist_guess, rotor_tip_twist_guess; length=nbe)

    # - Non-dimensional Radial Locations - #
    radial_positions = range(0.0, 1.0; length=nbe)

    # - TODO: add airfoil data - #
    airfoils = nothing

    # - Number of blades - #
    B = 5

    blade_elements, rotor_panels = dt.generate_blade_elements(
        rotor_c4_pos, radial_positions, chords, twists, airfoils, B, body_geometry
    )

    #---------------------------------#
    #             Poblem              #
    #---------------------------------#

    method = ff.AxisymmetricProblem(Vortex(Constant()), Neumann(), [false, true])
    problem = ff.define_problem(
        method, (duct_coordinates, hub_coordinates), 0.0, -1.0, -1.0
    )

    #---------------------------------#
    #             Meshes              #
    #---------------------------------#

    mesh_body_to_body = ff.generate_mesh(method, body_panels)

    mesh_bodies_to_rotor = dt.generate_one_way_mesh(body_panels, rotor_panels)

    #---------------------------------#
    #       Coefficient Matrices      #
    #---------------------------------#

    body_system = ff.generate_inviscid_system(method, body_panels, mesh_body_to_body)
    A_body_to_body = body_system.A
    bc_freestream_to_body = body_system.b

    A_bodies_to_rotor = dt.assemble_one_way_vortex_matrix(
        mesh_bodies_to_rotor, body_panels, rotor_panels
    )

    #---------------------------------#
    #      Singularity Strenghts      #
    #---------------------------------#

    vortex_strengths = ia.implicit_linear(A_body_to_body, bc_freestream_to_body)
    gamma_bodies = vortex_strengths[1:(end - 1)]

    #---------------------------------#
    #        Induced Velocities       #
    #---------------------------------#

    Vinf = 1.0
    body_induced_rotor_velocity = A_bodies_to_rotor * (gamma_bodies .* Vinf) .+ Vinf

    if debug
        return body_induced_rotor_velocity, ff.post_process(method,problem,
                                                            body_panels, mesh_body_to_body, ff.solve(body_system))
    else
        return body_induced_rotor_velocity
    end
end
