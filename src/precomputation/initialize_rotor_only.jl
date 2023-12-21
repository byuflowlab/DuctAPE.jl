"""
rotor_parameters = rotor_paramters
freestream = freestream

currently only capable of a single rotor
"""
function precompute_inputs_rotor_only(
    rotor_parameters, freestream; wake_length=5.0, wake_x_refine=0.05
)

    # TODO: udpate if needed:
    TF = Float64
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
        generate_rotor_panels(rotor_parameters[i].rotorzloc, rotor_panel_edges[:, i]) for
        i in 1:nrotor
    ]

    #get rotor panel radial center points for convenience
    rotor_panel_centers = rotor_source_panels[1].controlpoint[:, 2]
    nbe = length(rotor_panel_centers)

    ##### ----- Blade Elements ----- #####
    blade_elements = [
        generate_blade_elements(
            rotor_parameters[i].B,
            rotor_parameters[i].Omega,
            rotor_parameters[i].rotorzloc,
            rotor_parameters[i].r,
            rotor_parameters[i].chords,
            rotor_parameters[i].twists,
            rotor_parameters[i].airfoils,
            rotor_parameters[i].Rtip,
            rotor_parameters[i].Rhub,
            rotor_source_panels[i].controlpoint[:, 2],
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
        rotor_parameters[1].rotorzloc, rotor_parameters[1].rotorzloc + wake_length * D; step=dr
    )

    num_wake_x_panels = length(xrange) - 1

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

    wakeK = get_wake_k(wake_vortex_panels)

    # Go through the wake panels and determine the index of the aftmost rotor infront and the blade node from which the wake strength is defined.
    rotorwakeid = ones(Int, wake_vortex_panels.totpanel, 2)
    for i in 1:(rotor_parameters[1].nwake_sheets)
        rotorwakeid[(1 + (i - 1) * num_wake_x_panels):(i * num_wake_x_panels), 1] .= i
    end
    for (i, cp) in enumerate(eachrow(wake_vortex_panels.controlpoint))
        rotorwakeid[i, 2] = findlast(x -> x < cp[1], rotor_parameters.rotorzloc)
    end

    ######################################################################
    #                                                                    #
    #                  PRECOMPUTE UNIT INDUCED VELOCITIES                #
    #                                                                    #
    ######################################################################

    #---------------------------------#
    #           Velocities            #
    #---------------------------------#

    # - wake to rotor - #
    v_rw = [
        influencefromvortexpanels(
            rotor_source_panels[i].controlpoint,
            wake_vortex_panels.controlpoint,
            wake_vortex_panels.len,
            ones(TF, wake_vortex_panels.totpanel),
        ) for i in 1:length(rotor_source_panels)
    ]

    # axial components
    vz_rw = [v_rw[i][:, :, 1] for i in 1:length(rotor_source_panels)]

    # radial components
    vr_rw = [v_rw[i][:, :, 2] for i in 1:length(rotor_source_panels)]

    ##### ----- Rotor to Rotor ----- #####
    # - rotor to rotor - #
    v_rr = [
        influencefromsourcepanels(
            rotor_source_panels[i].controlpoint,
            rotor_source_panels[j].controlpoint,
            rotor_source_panels[j].len,
            ones(TF, rotor_source_panels[j].totpanel),
        ) for i in 1:length(rotor_source_panels), j in 1:length(rotor_source_panels)
    ]

    # axial components
    vz_rr = [
        v_rr[i, j][:, :, 1] for i in 1:length(rotor_source_panels),
        j in 1:length(rotor_source_panels)
    ]

    # radial components
    vr_rr = [
        v_rr[i, j][:, :, 2] for i in 1:length(rotor_source_panels),
        j in 1:length(rotor_source_panels)
    ]

    # - rotor to wake - #
    v_wr = [
        influencefromsourcepanels(
            wake_vortex_panels.controlpoint,
            rotor_source_panels[j].controlpoint,
            rotor_source_panels[j].len,
            ones(TF, rotor_source_panels[j].totpanel),
        ) for j in 1:length(rotor_source_panels)
    ]

    # axial components
    vz_wr = [v_wr[j][:, :, 1] for j in 1:length(rotor_source_panels)]

    # radial components
    vr_wr = [v_wr[j][:, :, 2] for j in 1:length(rotor_source_panels)]

    # - wake to wake - #
    v_ww = influencefromvortexpanels(
        wake_vortex_panels.controlpoint,
        wake_vortex_panels.controlpoint,
        wake_vortex_panels.len,
        ones(TF, wake_vortex_panels.totpanel),
    )

    # axial components
    vz_ww = v_ww[:, :, 1]

    # radial components
    vr_ww = v_ww[:, :, 2]

    inputs = (;
        converged=[false],
        Vinf=freestream.Vinf,
        freestream,
        # nxwake=length(xrange) - 1,
        num_wake_x_panels,
        wakeK,
        rotorwakeid,
        vz_rw=vz_rw,
        vr_rw=vr_rw,
        vz_rr=vz_rr,
        vr_rr=vr_rr,
        vz_wr=vz_wr,
        vr_wr=vr_wr,
        vz_ww=vz_ww,
        vr_ww=vr_ww,
        blade_elements,
        num_rotors=1,
        rotor_panel_edges=rotor_panel_edges,
        rotor_panel_centers=rotor_panel_centers,
        wake_vortex_panels,
        rotor_source_panels,
        xrange,
        x_grid_points,
        r_grid_points,
        nrotor_panels=sum(rotor_source_panels.totpanel),
        nwake_panels=wake_vortex_panels.totpanel,
        ductwakeinterfaceid=nothing,
        hubwakeinterfaceid=nothing,
    )

    return inputs
end

function initilize_states_rotor_only(inputs)

    # - Initialize with freestream only - #
    Wtheta = -inputs.rotor_panel_centers .* inputs.blade_elements.Omega'
    # use freestream magnitude as meridional velocity at each blade section
    Wm = similar(Wtheta) .= inputs.Vinf
    # magnitude is simply freestream and rotation
    W = sqrt.(Wtheta .^ 2 .+ Wm .^ 2)

    # initialize circulation and source panel strengths
    Gamr, sigr = calculate_gamma_sigma(inputs.blade_elements, Wm, Wtheta, W, inputs.freestream)

    nwake = inputs.wake_vortex_panels.totpanel
    gamw = zeros(nwake)
    calculate_wake_vortex_strengths!(
        gamw, Gamr, inputs.Vinf * ones(length(gamw)), inputs; debug=false
    )

    # gamw = initialize_wake_vortex_strengths(
    #     inputs.Vinf,
    #     Gamr,
    #     inputs.blade_elements.Omega,
    #     inputs.blade_elements.B,
    #     inputs.rotor_panel_edges,
    #     inputs.num_wake_x_panels,
    # )

    # - Set Up States and Parameters - #
    # initialize rotor source strengths to zero for now
    return [Gamr; gamw; sigr]
end
