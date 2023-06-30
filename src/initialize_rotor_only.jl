"""
rotor_parameters = rotor_paramters
freestream = freestream

currently only capable of a single rotor
"""
function initialize_rotor_states(
    rotor_parameters, freestream; wake_length=5.0, wake_x_refine=0.01
)
    #---------------------------------#
    #              Rotor              #
    #---------------------------------#

    nrotor = length(rotor_parameters)
    nbe = rotor_parameters[1].nwake_sheets - 1
    ##### ----- Panels ----- #####
    # - Rotor Panel Edges - #
    # rotor_panel_edges = range(rotor_parameters.Rhub, rotor_parameters.Rtip, rotor_parameters.nbe .+ 1)
    rotor_panel_edges = [
        range(rotor_parameters[i].Rhub, rotor_parameters[i].Rtip, nbe + 1) for i in 1:nrotor
    ]
    rotor_panel_edges = reduce(hcat, rotor_panel_edges)

    # - Generate Panels - #
    # note there will be nbe panels, but we require nbe+1 panel edges.
    rotor_source_panels = [
        generate_rotor_panels(rotor_parameters[i].xrotor, rotor_panel_edges[:, i]) for
        i in 1:nrotor
    ]

    #get rotor panel radial center points for convenience
    rotor_panel_centers = rotor_source_panels[1].panel_center[:, 2]
    nbe = length(rotor_panel_centers)

    ##### ----- Blade Elements ----- #####
    blade_elements = [
        generate_blade_elements(
            rotor_parameters[i].B,
            rotor_parameters[i].Omega,
            rotor_parameters[i].xrotor,
            rotor_parameters[i].r,
            rotor_parameters[i].chords,
            rotor_parameters[i].twists,
            rotor_parameters[i].airfoils,
            rotor_parameters[i].Rtip,
            rotor_parameters[i].Rhub,
            rotor_source_panels[i].panel_center[:, 2],
        ) for i in 1:nrotor
    ]

    #---------------------------------#
    #              Wake               #
    #---------------------------------#
    #=
    Trailing wake sheets extend from between the blade elements, in other words, from the rotor panel edge locations.
    =#

    ##### ----- General Geometry ----- #####

    # - Choose blade element spacing - #
    dr = wake_x_refine * rotor_parameters[1].Rtip

    # - Rotor Diameter - #
    D = 2.0 * rotor_parameters[1].Rtip

    # - Define x grid spacing - #
    #note that by using step rather than length, we aren't guarenteed to have the wake as long as we want, but it will be close.
    xrange = range(
        rotor_parameters[1].xrotor, rotor_parameters[1].xrotor + wake_length * D; step=dr
    )

    # - Put together the initial grid point matrices - #
    #=
    For the grid relaxation (which adjusts the grid points such that the grids lie along streamlines) the function expects the grid points to be formatted such that there are x rows, and r columns.
    In addition, we pass in 2 matrices, one for x and r such that all the x's are in one matrix, and all the r's in another.

    note, however, that we are not yet using the wake relaxation because we do not have any duct bodies to deform the streamlines, so we will assume that the wake extends straight back for now.
    NOTE THAT THIS MEANS WE AREN'T MODELING WAKE CONTRACTION, WHICH WILL LIKELY CAUSE DISCREPANCIES BETWEEN EXPERIMENT AND OUR MODEL, BUT THEY SHOULD BE REASONABLE.
    =#
    x_grid_points = repeat(xrange; outer=(1, nbe + 1))
    r_grid_points = repeat(rotor_panel_edges'; outer=(length(xrange), 1))

    ##### ----- Panels ----- #####

    #=
    In order to make sure there aren't erronious panels from the end of one wake sheet to the beginning of the next, we have a panel object for each wake sheet, rather than for the entire wake.
    Thus the wake_vortex_panels object is a vector of Panel structs, where each struct is comprised of a single row (sheet) of wake panels.
    We will need to remember this in indexing such that we call wake_vortex_panels[row].property[column]
    Note that there are nbe+1 trailing wake surfaces
    =#
    wake_vortex_panels = generate_wake_panels(x_grid_points, r_grid_points)

    ######################################################################
    #                                                                    #
    #                  PRECOMPUTE UNIT INDUCED VELOCITIES                #
    #                                                                    #
    ######################################################################

    #---------------------------------#
    #             MESHES              #
    #---------------------------------#
    #=
    Meshes contain the relative geometry used to calculate the unit induced velocities of the panels doing the influencing to the panels being affected.
    the first input is the influencing panel arrays
    the second input is the affected panel arrays
    the output is a matrix of [affect index, influence index] for the various panel arrays
    =#

    ##### ----- Wake to Rotor ----- #####
    wake_to_rotor_mesh = [
        generate_one_way_mesh(wake_vortex_panels[j], rotor_source_panels[i]) for
        i in 1:length(rotor_source_panels), j in 1:length(wake_vortex_panels)
    ]

    ##### ----- Rotor to Rotor ----- #####
    rotor_to_rotor_mesh = [
        generate_one_way_mesh(rotor_source_panels[j], rotor_source_panels[i]) for
        i in 1:length(rotor_source_panels), j in 1:length(rotor_source_panels)
    ]

    #---------------------------------#
    #           Velocities            #
    #---------------------------------#
    #=
    To get the induced velocities, we take the meshes just created, then input the influencing panel arrays again, then the affected panel arrays.
    The output is a matrix again with index [affected, influencing]

    We then separate out the inputs into the x and r components
    =#

    ##### ----- Wake to Rotor ----- #####
    A_wake_to_rotor = [
        assemble_induced_velocity_matrices(
            wake_to_rotor_mesh[i, j], wake_vortex_panels[j], rotor_source_panels[i]
        ) for i in 1:length(rotor_source_panels), j in 1:length(wake_vortex_panels)
    ]

    # - Axial - #
    vx_rw = [
        A_wake_to_rotor[i, j][1] for i in 1:length(rotor_source_panels),
        j in 1:length(wake_vortex_panels)
    ]

    # - Radial - #
    vr_rw = [
        A_wake_to_rotor[i, j][2] for i in 1:length(rotor_source_panels),
        j in 1:length(wake_vortex_panels)
    ]

    ##### ----- Rotor to Rotor ----- #####
    A_rotor_to_rotor = [
        assemble_induced_velocity_matrices(
            rotor_to_rotor_mesh[i, j],
            rotor_source_panels[j],
            rotor_source_panels[i];
            singularity="source",
        ) for i in 1:length(rotor_source_panels), j in 1:length(rotor_source_panels)
    ]

    # - Axial - #
    vx_rr = [
        A_rotor_to_rotor[i, j][1] for i in 1:length(rotor_source_panels),
        j in 1:length(rotor_source_panels)
    ]

    # - Radial - #
    vr_rr = [
        A_rotor_to_rotor[i, j][2] for i in 1:length(rotor_source_panels),
        j in 1:length(rotor_source_panels)
    ]

    #################################

    # - Initialize with freestream only - #
    Wθ = -rotor_panel_centers .* rotor_parameters.Omega'
    # use freestream magnitude as meridional velocity at each blade section
    Wm = similar(Wθ) .= freestream.Vinf
    # magnitude is simply freestream and rotation
    W = sqrt.(Wθ .^ 2 .+ Wm .^ 2)
    # initialize circulation and source panel strengths
    Gamma, sigr = calculate_gamma_sigma(blade_elements, Wm, Wθ, W)

    gamma_wake = initialize_wake_vortex_strengths(
        freestream.Vinf,
        Gamma,
        rotor_parameters.Omega,
        rotor_parameters.B,
        rotor_panel_edges,
    )

    # - Set Up States and Parameters - #
    # for ir in 1:length(rotor_panel_centers)
    #     for irotor in 1:length(rp)
    #         Wθ[ir,irotor] = -rotor_panel_centers[ir] * rp[irotor].Omega[1]
    #     end
    # end

    # use freestream magnitude as meridional velocity at each blade section
    Wm = similar(Wθ) .= fs.Vinf
    # magnitude is simply freestream and rotation
    W = sqrt.(Wθ .^ 2 .+ Wm .^ 2)
    # initialize circulation and source panel strengths
    Gamma, sigr = calculate_gamma_sigma(blade_elements, Wm, Wθ, W)

    gamma_wake = initialize_wake_vortex_strengths(
        fs.Vinf, Gamma, rp.Omega, rp.B, rotor_panel_edges
    )

    # - Set Up States and Parameters - #
    # initialize rotor source strengths to zero for now
    states = [Gamma; gamma_wake; sigr]

    params = (
        converged=[false],
        Vinf=freestream.Vinf,
        nxwake=length(xrange) - 1,
        vx_rw=vx_rw,
        vr_rw=vr_rw,
        vx_rr=vx_rr,
        vr_rr=vr_rr,
        blade_elements,
        num_rotors=1,
        rotor_panel_edges=rotor_panel_edges,
        rotor_panel_centers=rotor_panel_centers,
        wake_vortex_panels,
        rotor_source_panels,
        xrange,
        x_grid_points,
        r_grid_points,
    )

    return states, params
end
