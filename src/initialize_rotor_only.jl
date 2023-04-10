"""
rp = rotor_paramters
fs = freestream

currently only capable of a single rotor
"""
function initialize_rotor_states(rp, fs)
    #---------------------------------#
    #              Rotor              #
    #---------------------------------#

    ##### ----- Panels ----- #####
    # - Rotor Panel Edges - #
    rotor_panel_edges = range(rp.Rhub, rp.Rtip, rp.nbe + 1)

    # - Generate Panels - #
    # note there will be nbe panels, but we require nbe+1 panel edges.
    rotor_source_panels = [generate_rotor_panels(rp.xrotor, rotor_panel_edges)]

    #get rotor panel radial center points for convenience
    rbe = rotor_source_panels[1].panel_center[:, 2]
    nbe = length(rbe)

    ##### ----- Blade Elements ----- #####
    blade_elements = generate_blade_elements(
        rp.B,
        rp.Omega,
        rp.xrotor,
        rp.r,
        rp.chord,
        rp.twist,
        rp.airfoil,
        rp.Rtip,
        rp.Rhub,
        rbe,
    )

    #---------------------------------#
    #              Wake               #
    #---------------------------------#
    #=
    Trailing wake sheets extend from between the blade elements, in other words, from the rotor panel edge locations.
    =#

    ##### ----- General Geometry ----- #####

    # - Choose blade element spacing - #
    dr = 0.001

    # - Rotor Diameter - #
    D = 2.0 * rp.Rtip

    # - Define x grid spacing - #
    #note that by using step rather than length, we aren't guarenteed to have the wake as long as we want, but it will be close.
    xrange = range(0.0, 5.0 * D; step=dr)

    # - Put together the initial grid point matrices - #
    #=
    For the grid relaxation (which adjusts the grid points such that the grids lie along streamlines) the function expects the grid points to be formatted such that there are x rows, and r columns.
    In addition, we pass in 2 matrices, one for x and r such that all the x's are in one matrix, and all the r's in another.

    note, however, that we are not yet using the wake relaxation because we do not have any duct bodies to deform the streamlines, so we will assume that the wake extends straight back for now.
    NOTE THAT THIS MEANS WE AREN'T MODELING WAKE CONTRACTION, WHICH WILL LIKELY CAUSE DISCREPANCIES BETWEEN EXPERIMENT AND OUR MODEL, BUT THEY SHOULD BE REASONABLE.
    =#
    x_grid_points = repeat(xrange; outer=(1, rp.nbe + 1))
    r_grid_points = repeat(rotor_panel_edges'; outer=(length(xrange), 1))

    ##### ----- Panels ----- #####

    # For the wake panels, we want to use Constant Vortex Panels.
    # The boundary condition and body of revolution boolean don't matter, should probably TODO: add a constructor that doesn't require those when not needed.
    wake_method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true])

    #=
    In order to make sure there aren't erronious panels from the end of one wake sheet to the beginning of the next, we have a panel object for each wake sheet, rather than for the entire wake.
    Thus the wake_vortex_panels object is a vector of Panel structs, where each struct is comprised of a single row (sheet) of wake panels.
    We will need to remember this in indexing such that we call wake_vortex_panels[row].property[column]
    Note that there are nbe+1 trailing wake surfaces
    =#
    wake_vortex_panels = generate_wake_panels(
        x_grid_points, r_grid_points; method=wake_method
    )

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
    vxd_wake_to_rotor = [
        A_wake_to_rotor[i, j][1] for i in 1:length(rotor_source_panels),
        j in 1:length(wake_vortex_panels)
    ]

    # - Radial - #
    vrd_wake_to_rotor = [
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
    vxd_rotor_to_rotor = [
        A_rotor_to_rotor[i, j][1] for i in 1:length(rotor_source_panels),
        j in 1:length(rotor_source_panels)
    ]

    # - Radial - #
    vrd_rotor_to_rotor = [
        A_rotor_to_rotor[i, j][2] for i in 1:length(rotor_source_panels),
        j in 1:length(rotor_source_panels)
    ]

    #################################

    # - Initialize with freestream only - #
    W = sqrt.(fs.Vinf^2 .+ (rp.Omega * rbe) .^ 2)
    #TODO: need to use the blade element airfoil stuff here (look at rotor aero file for usage)
    aoa = atan.(fs.Vinf, rp.Omega * rbe)

    cl = zeros(rp.nbe)
    cd = zeros(rp.nbe)
    for i in 1:(rp.nbe)
        clin, cdin = search_polars(blade_elements.inner_airfoil[i], aoa[i])
        clout, cdout = search_polars(blade_elements.outer_airfoil[i], aoa[i])
        # linearly interpolate between those two values at your blade element location
        cl[i] = fm.linear([0.0; 1.0], [clin, clout], blade_elements.inner_fraction[i])
        cd[i] = fm.linear([0.0; 1.0], [cdin, cdout], blade_elements.inner_fraction[i])
    end

    Gamma_init = 0.5 .* W .* cl .* blade_elements.chords
    Vmr = fs.Vinf * ones(length(W))

    Gamma_tilde_init = rp.B .* Gamma_init
    H_tilde_init = rp.Omega * rp.B * Gamma_init / (2.0 * pi)

    # - initialize from open rotor - #
    Vm_init = vm_from_vinf(fs.Vinf, Gamma_tilde_init, H_tilde_init, rbe)
    gamma_wake_init = [
        Vm_init[1] - fs.Vinf
        Vm_init[2:end] .- Vm_init[1:(end - 1)]
        fs.Vinf - Vm_init[end]
    ]

    # - Set Up States and Parameters - #
    # initialize rotor source strengths to zero for now
    states = [Gamma_init; gamma_wake_init; zeros(nbe)]
    # states = [Gamma_init; gamma_wake_init]

    params = (
        converged=[false],
        Vinf=fs.Vinf,
        nxwake=length(xrange) - 1,
        vx_rw=vxd_wake_to_rotor,
        vr_rw=vrd_wake_to_rotor,
        vx_rr=vxd_rotor_to_rotor,
        vr_rr=vrd_rotor_to_rotor,
        blade_elements=[blade_elements],
        num_rotors=1,
        rotor_panel_edges=rotor_panel_edges,
        rotor_panel_centers=rbe,
    )

    return states, params
end

"""
get meridional velocities using the marching method for when Vinf outside of wake is known.
"""
function vm_from_vinf(Vinf, Gamma_tilde, H_tilde, radial_stations)

    #rename for convenience
    nr = length(Gamma_tilde)

    #initialize
    Vm = zeros(nr)

    #march from just outside the tip to the hub
    for i in nr:-1:1
        if i == nr

            #this is the change at the tip, so we need to set Vm2 to Vinf, and set Gamma dn H 2's to zeros
            radical =
                Vinf^2 + (1.0 / (2.0 * pi * radial_stations[i]))^2 * (-Gamma_tilde[i]^2) -
                2.0 * (-H_tilde[i])
            # if radical > 0.0
            Vm[i] = sqrt(radical)
            # else
            #     Vm[i] = 0.0
            # end
        else

            #otherwise, we just take the differences inside as-is
            radical =
                Vm[i + 1]^2 +
                (1.0 / (2.0 * pi * radial_stations[i]))^2 *
                (Gamma_tilde[i + 1]^2 - Gamma_tilde[i]^2) -
                2.0 * (H_tilde[i + 1] - H_tilde[i])
            # if radical > 0.0
            Vm[i] = sqrt(radical)
            # else
            #     Vm[i] = 0.0
            # end
        end
    end

    #Vm should now be in the order of hub to tip
    return Vm
end
