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
    rotor_c4_pos
    stator_c4_pos
    stator_root_chord
    stator_tip_chord
    stator_twist_guess
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
    rotor_c4_pos = x[19]
    stator_c4_pos = x[20]
    stator_root_chord = x[21]
    stator_tip_chord = x[22]
    stator_twist_guess = x[23]

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
    nbe = 5
    nbe_fine = 10

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
    airfoils = [nothing for i in 1:nbe]

    # - Number of blades - #
    B = 27

    rotor_blade_elements, rotor_panels = dt.generate_blade_elements(
        rotor_c4_pos, radial_positions, chords, twists, airfoils, nbe_fine, B, body_geometry
    )

    # - Get Stator Geometry - #
    # - Chord - #
    # note: just use linear chord distribution
    chords = range(stator_root_chord, stator_tip_chord; length=nbe)

    # - Twist - #
    # note: need to make a guess, use degrees for input
    twists = stator_twist_guess * ones(nbe)

    # - Non-dimensional Radial Locations - #
    radial_positions = range(0.0, 1.0; length=nbe)

    # - TODO: add airfoil data - #
    airfoils = [nothing for i in 1:nbe]

    # - Number of blades - #
    B = 8

    #TODO: actually put this after the wake generation
    #TODO: also need to update this function or add another that takes in updated radial positions rather than just a number of blade elements to refine to.
    stator_blade_elements, stator_panels = dt.generate_blade_elements(
        stator_c4_pos,
        radial_positions,
        chords,
        twists,
        airfoils,
        nbe_fine,
        B,
        body_geometry,
    )

    #---------------------------------#
    #             Poblem              #
    #---------------------------------#

    method = ff.AxisymmetricProblem(Vortex(Constant()), Neumann(), [false, true])
    problem = ff.define_problem(
        method, (duct_coordinates, hub_coordinates), 0.0, -1.0, -1.0
    )

    #---------------------------------#
    #       Initial Wake Points       #
    #---------------------------------#

    x_grid_points, r_grid_points, nx, nr, rotoridxs, wake_panels = dt.generate_wake_grid(
        body_geometry,
        [rotor_blade_elements; stator_blade_elements]; #TODO: actually just need x position for each rotor and radial locations for front most rotor
        wake_length=1.0,
        debug=false,
    )

    #TODO: need to re-interpolate the stator blade elements based on the grid radial positions at the rotoridx for the stator.

    #---------------------------------#
    #             Meshes              #
    #---------------------------------#

    # - Body -> Body - #
    mesh_body_to_body = ff.generate_mesh(method, body_panels)

    # - Body -> Rotor - #
    mesh_bodies_to_rotor = dt.generate_one_way_mesh(body_panels, rotor_panels)

    # - Body -> Stator - #
    mesh_bodies_to_stator = dt.generate_one_way_mesh(body_panels, stator_panels)

    # - Wake -> Body - #
    mesh_wake_to_body = dt.generate_one_way_mesh(wake_panels, body_panels)

    # - Wake -> Rotor - #
    mesh_wake_to_rotor = dt.generate_one_way_mesh(wake_panels, rotor_panels)

    # - Wake -> Stator - #
    mesh_wake_to_stator = dt.generate_one_way_mesh(wake_panels, stator_panels)

    # - Rotor -> Body - #
    mesh_rotor_to_body = dt.generate_one_way_mesh(
        rotor_panels, body_panels; singularity="source"
    )

    # - Rotor -> Rotor - #
    mesh_rotor_to_rotor = dt.generate_one_way_mesh(
        rotor_panels, rotor_panels; singularity="source"
    )

    # - Rotor -> Stator - #
    mesh_rotor_to_stator = dt.generate_one_way_mesh(
        rotor_panels, stator_panels; singularity="source"
    )

    # - Stator -> Body - #
    mesh_stator_to_body = dt.generate_one_way_mesh(
        stator_panels, body_panels; singularity="source"
    )

    # - Stator -> Rotor - #
    mesh_stator_to_rotor = dt.generate_one_way_mesh(
        stator_panels, rotor_panels; singularity="source"
    )

    # - Stator -> Stator - #
    mesh_stator_to_stator = dt.generate_one_way_mesh(
        stator_panels, stator_panels; singularity="source"
    )

    #---------------------------------#
    #       Coefficient Matrices      #
    #---------------------------------#

    # - Body -> Body - #
    body_system = ff.generate_inviscid_system(method, body_panels, mesh_body_to_body)
    A_body_to_body = body_system.A
    bc_freestream_to_body = body_system.b

    # - Body -> Rotor - #
    A_bodies_to_rotor = dt.assemble_one_way_coefficient_matrix(
        mesh_bodies_to_rotor, body_panels, rotor_panels
    )

    # - Body -> Stator - #
    A_bodies_to_stator = dt.assemble_one_way_coefficient_matrix(
        mesh_bodies_to_stator, body_panels, stator_panels
    )

    # - Wake -> Body - #
    A_wake_to_bodies = dt.assemble_one_way_coefficient_matrix(
        mesh_wake_to_body, wake_panels, body_panels
    )

    # - Wake -> Rotor - #
    A_wake_to_rotor = dt.assemble_one_way_coefficient_matrix(
        mesh_wake_to_rotor, wake_panels, stator_panels
    )

    # - Wake -> Stator - #
    A_wake_to_stator = dt.assemble_one_way_coefficient_matrix(
        mesh_wake_to_stator, wake_panels, stator_panels
    )

    # - Rotor -> Body - #
    A_rotor_to_body = dt.assemble_one_way_coefficient_matrix(
        mesh_rotor_to_body, rotor_panels, body_panels; singularity="source"
    )

    # - Rotor -> Rotor - #
    A_rotor_to_rotor = dt.assemble_one_way_coefficient_matrix(
        mesh_rotor_to_rotor, rotor_panels, rotor_panels; singularity="source"
    )

    # - Rotor -> Stator - #
    A_rotor_to_stator = dt.assemble_one_way_coefficient_matrix(
        mesh_rotor_to_stator, rotor_panels, stator_panels; singularity="source"
    )

    # - Stator -> Body - #
    A_stator_to_body = dt.assemble_one_way_coefficient_matrix(
        mesh_stator_to_body, stator_panels, body_panels; singularity="source"
    )

    # - Stator -> Rotor - #
    A_stator_to_rotor = dt.assemble_one_way_coefficient_matrix(
        mesh_stator_to_rotor, stator_panels, rotor_panels; singularity="source"
    )

    # - Stator -> Stator - #
    A_stator_to_stator = dt.assemble_one_way_coefficient_matrix(
        mesh_stator_to_stator, stator_panels, stator_panels; singularity="source"
    )

    #---------------------------------#
    #      Singularity Strengths      #
    #---------------------------------#

    body_vortex_strengths = ia.implicit_linear(A_body_to_body, bc_freestream_to_body)
    gamma_bodies = body_vortex_strengths[1:(end - 1)]

    #---------------------------------#
    #        Induced Velocities       #
    #---------------------------------#

    Vinf = 1.0
    body_induced_rotor_velocity = A_bodies_to_rotor * (gamma_bodies .* Vinf) .+ Vinf
    body_induced_stator_velocity = A_bodies_to_stator * (gamma_bodies .* Vinf) .+ Vinf

    if debug
        return x_grid_points, r_grid_points
    else
        return wake_panels[2].panel_center[:, 2]
    end
end
