println("\nPRE-PROCESSING TESTS")

@testset "Parameter Precomputation" begin
    include("data/basic_two_rotor_for_test_NEW.jl")

    pd = dt.get_problem_dimensions(paneling_constants)

    # discretize wake
    zwake, rotor_indices_in_wake = dt.discretize_wake(
        duct_coordinates,
        centerbody_coordinates,
        rotorstator_parameters.rotorzloc, # rotor axial locations
        paneling_constants.wake_length,
        paneling_constants.npanels,
        paneling_constants.dte_minus_cbte;
    )

    @test zwake == [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
    @test rotor_indices_in_wake == [1, 3]

    # reinterpolate_bodies

    rp_duct_coordinates = zeros(2, pd.ndn)
    rp_centerbody_coordinates = zeros(2, pd.ncbn)
    dt.reinterpolate_bodies!(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        duct_coordinates,
        centerbody_coordinates,
        zwake,
        paneling_constants.ncenterbody_inlet,
        paneling_constants.nduct_inlet;
        finterp=FLOWMath.linear,
    )

    @test size(rp_duct_coordinates, 1) < size(rp_duct_coordinates, 2)
    @test size(rp_centerbody_coordinates, 1) < size(rp_centerbody_coordinates, 2)
    @test size(rp_duct_coordinates, 2) == 2 * size(duct_coordinates, 1) - 1
    @test size(rp_centerbody_coordinates, 2) == 2 * size(centerbody_coordinates, 1) - 1

    @test rp_duct_coordinates == [
        1.0 0.75 0.5 0.25 0.0 0.25 0.5 0.75 1.0
        2.0 1.75 1.5 1.75 2.0 2.25 2.5 2.25 2.0
    ]
    @test rp_centerbody_coordinates == [
        0.0 0.25 0.5 0.75 1.0
        0.0 0.25 0.5 0.25 0.0
    ]

    rpb4 = copy(rp_duct_coordinates)

    dt.place_duct!(
        rp_duct_coordinates,
        rotorstator_parameters.Rtip[1],
        rotorstator_parameters.rotorzloc[1],
        rotorstator_parameters.tip_gap[1],
    )

    @test rp_duct_coordinates[1, :] == rpb4[1, :]
    @test rp_duct_coordinates[2, :] == rpb4[2, :] .- 0.75

    # get_blade_ends_from_body_geometry
    Rtips, Rhubs = dt.get_blade_ends_from_body_geometry(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        rotorstator_parameters.tip_gap,
        rotorstator_parameters.rotorzloc,
    )

    @test all(Rtips .== 1.0)
    @test all(Rhubs .== 0.25)

    # generate_wake_grid
    grid = dt.generate_wake_grid(
        pd,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        Rhubs[1],
        Rtips[1],
        rotorstator_parameters.tip_gap[1],
        zwake;
        wake_max_iter=100,
        wake_nlsolve_ftol=1e-9,
        verbose=false,
        silence_warnings=true,
    )

    @test isapprox(
        grid[:, :, 1],
        [
            0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0
            0.25 0.5 0.25 0.0 0.0 0.0 0.0 0.0
        ],
    )
    @test isapprox(
        grid[:, :, 2],
        [
            0.25 0.49999934235473875 0.7491867275686039 1.0003511147926796 1.2509282467245646 1.4999107373121212 1.7498120063838813 1.9999729168852773
            0.625 0.6026658944047609 0.6063118606432654 0.6760730328966827 0.6953337735090016 0.7118731817894995 0.7241605597510481 0.7283183433975138
        ],
        atol=1e-6,
    )
    @test grid[:, :, 3] == [
        0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0
        1.0 0.75 1.0 1.25 1.25 1.25 1.25 1.25
    ]

    # - All together - #
    wake_grid, rp_duct_coordinates, rp_centerbody_coordinates, rotor_indices_in_wake = dt.reinterpolate_geometry(
        pd,
        duct_coordinates,
        centerbody_coordinates,
        rotorstator_parameters,
        paneling_constants;
        wake_max_iter=100,
        wake_nlsolve_ftol=1e-9,
        autoshiftduct=true,
        finterp=FLOWMath.linear,
        verbose=false,
        silence_warnings=true,
    )

    @test all(grid .== wake_grid)
    @test rotor_indices_in_wake == [1, 3]

    # all panels together
    body_vortex_panels, rotor_source_panels, wake_vortex_panels = dt.generate_all_panels(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        paneling_constants.nwake_sheets,
        rotor_indices_in_wake,
        rotorstator_parameters.rotorzloc,
        wake_grid;
        itcpshift=0.05,
        axistol=1e-15,
        tegaptol=1e1 * eps(),
        silence_warnings=true,
    )

    #TODO: add some sort of test here, eventually add the in place version and compare probably

    # - Unit Induced Velocities - #
    ivr, ivw, ivb = dt.calculate_unit_induced_velocities(
        pd, (; body_vortex_panels, rotor_source_panels, wake_vortex_panels)
    )

    # - Initialize Linear System - #
    linsys = dt.initialize_linear_system(
        ivb,
        body_vortex_panels,
        rotor_source_panels,
        wake_vortex_panels,
        operating_point.Vinf,
    )

    @test isapprox(
        linsys.b_bf,
        operating_point.Vinf[1] * [body_vortex_panels.normal[1, :]; 1.0; zeros(3)],
    )

    # - Blade element interpolation - #

    blade_elements = dt.interpolate_blade_elements(
        rotorstator_parameters, Rtips, Rhubs, rotor_source_panels.controlpoint[2, :], pd.nbe
    )

    @test isapprox(
        blade_elements.inner_fraction,
        [0.75 0.7126237212865307; 0.25 0.21262372128653073],
        atol=1e-6,
    )
    @test blade_elements.stagger == 70 * ones(2, 2) * pi / 180
    @test blade_elements.rotor_panel_centers ==
        reshape(rotor_source_panels.controlpoint[2, :], (2, 2))
    @test blade_elements.B == [2, 4]
    @test blade_elements.chords == 0.1 * ones(2, 2)
    @test blade_elements.twists == 20 * ones(2, 2) * pi / 180
    @test isapprox(
        blade_elements.solidity,
        [0.07275654541343787 0.14868876670453912; 0.039176601376466544 0.07926477889699948],
        atol=1e-6,
    )
    @test blade_elements.fliplift == [0.0, 0.0]
    @test blade_elements.Rhub == rotorstator_parameters.Rhub
    @test blade_elements.Rtip == rotorstator_parameters.Rtip

    # - Index Mapping - #
    idmaps = dt.set_index_maps(
        paneling_constants.npanels,
        paneling_constants.ncenterbody_inlet,
        paneling_constants.nwake_sheets,
        paneling_constants.dte_minus_cbte,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.endnodeidxs,
        pd.nwp,
        pd.nwsp,
        pd.nbn,
        rotor_indices_in_wake,
        pd.nrotor,
    )

    @test idmaps.body_totnodes == body_vortex_panels.totnode
    @test idmaps.wake_nodemap == wake_vortex_panels.nodemap
    @test idmaps.wake_endnodeidxs == wake_vortex_panels.endnodeidxs
    @test idmaps.rotor_indices_in_wake == rotor_indices_in_wake
    @test idmaps.wake_node_sheet_be_map[:, 1] ==
        [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3]
    @test idmaps.wake_node_sheet_be_map[:, 2] ==
        [1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2]
    @test idmaps.wake_node_ids_along_casing_wake_interface == [15, 16, 17, 18, 19]
    @test idmaps.wake_node_ids_along_centerbody_wake_interface == [1, 2, 3, 4]

    # - EVERYTHING all together - #
    tivr, tivw, tivb, tlinsys, tblade_elements, tidmaps, tpanels = dt.precompute_parameters_iad(
        propulsor;
        wake_max_iter=100,
        wake_nlsolve_ftol=1e-9,
        autoshiftduct=true,
        itcpshift=0.05,
        axistol=1e-15,
        tegaptol=1e1 * eps(),
        finterp=FLOWMath.linear,
        silence_warnings=true,
        verbose=false,
    )

    @test all(compare_namedtuples(tivr, ivr))
    @test all(compare_namedtuples(tivw, ivw))
    @test all(compare_namedtuples(tivb, ivb))
    @test all(compare_namedtuples(tlinsys, linsys))
    @test all(compare_namedtuples(tblade_elements, blade_elements))
    @test all(compare_namedtuples(tidmaps, idmaps))
    @test all(compare_namedtuples(tpanels.body_vortex_panels, body_vortex_panels))
    @test all(compare_namedtuples(tpanels.rotor_source_panels, rotor_source_panels))
    @test all(compare_namedtuples(tpanels.wake_vortex_panels, wake_vortex_panels))
end

@testset "other idmap cases" begin

    # hub is longer than duct
    npanels = [2, 1, 1] #npanels
    ncenterbody_inlet = 2 # paneling_constants.ncenterbody_inlet
    nwake_sheets = 3 # paneling_constants.nwake_sheets
    dte_minus_cbte = -1 # paneling_constants.dte_minus_cbte
    wake_nodemape = [
        1 2 3 4 6 7 8 9 11 12 13 14
        2 3 4 5 7 8 9 10 12 13 14 15
    ]# wake_vortex_panels.nodemap
    wake_endnodeidxs = [1 6 11; 5 10 15] # wake_vortex_panels.endnodeidxs
    nwp = 12 # pd.nwp
    nwsp = 4 # pd.nwsp
    nbn = 12 # pd.nbn
    rotor_indices_in_wake = [1] # rotor_indices_in_wake
    nrotor = 1 # pd.nrotor

    idmaps = dt.set_index_maps(
        npanels,
        ncenterbody_inlet,
        nwake_sheets,
        dte_minus_cbte,
        wake_nodemap,
        wake_endnodeidxs,
        nwp,
        nwsp,
        nbn,
        rotor_indices_in_wake,
        nrotor,
    )

    @test idmaps.wake_node_sheet_be_map[:, 1] ==
        [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3]
    @test all(idmaps.wake_node_sheet_be_map[:, 2] .== 1)
    @test idmaps.wake_panel_sheet_be_map[:, 1] == [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
    @test all(idmaps.wake_panel_sheet_be_map[:, 2] .== 1)

    @test idmaps.wake_node_ids_along_casing_wake_interface == [11, 12, 13]
    @test idmaps.wake_node_ids_along_centerbody_wake_interface == [1, 2, 3, 4]

    @test idmaps.wake_panel_ids_along_casing_wake_interface == [9, 10]
    @test idmaps.wake_panel_ids_along_centerbody_wake_interface == [1, 2, 3]

    @test idmaps.centerbody_panel_ids_along_centerbody_wake_interface == [3, 4, 5]
    @test idmaps.duct_panel_ids_along_centerbody_wake_interface == [2, 1]

    @test idmaps.id_of_first_casing_panel_aft_of_each_rotor == [2]
    @test idmaps.id_of_first_centerbody_panel_aft_of_each_rotor == [3]

    # hub is shorter than duct
    npanels = [2, 1, 1] #npanels
    ncenterbody_inlet = 2 # paneling_constants.ncenterbody_inlet
    nwake_sheets = 3 # paneling_constants.nwake_sheets
    dte_minus_cbte = 1 # paneling_constants.dte_minus_cbte
    wake_nodemape = [
        1 2 3 4 6 7 8 9 11 12 13 14
        2 3 4 5 7 8 9 10 12 13 14 15
    ]# wake_vortex_panels.nodemap
    wake_endnodeidxs = [1 6 11; 5 10 15] # wake_vortex_panels.endnodeidxs
    nwp = 12 # pd.nwp
    nwsp = 4 # pd.nwsp
    nbn = 12 # pd.nbn
    rotor_indices_in_wake = [1] # rotor_indices_in_wake
    nrotor = 1 # pd.nrotor

    idmaps = dt.set_index_maps(
        npanels,
        ncenterbody_inlet,
        nwake_sheets,
        dte_minus_cbte,
        wake_nodemap,
        wake_endnodeidxs,
        nwp,
        nwsp,
        nbn,
        rotor_indices_in_wake,
        nrotor,
    )

    @test idmaps.wake_node_sheet_be_map[:, 1] ==
        [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3]
    @test all(idmaps.wake_node_sheet_be_map[:, 2] .== 1)
    @test idmaps.wake_panel_sheet_be_map[:, 1] == [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
    @test all(idmaps.wake_panel_sheet_be_map[:, 2] .== 1)

    @test idmaps.wake_node_ids_along_casing_wake_interface == [11, 12, 13, 14]
    @test idmaps.wake_node_ids_along_centerbody_wake_interface == [1, 2, 3]

    @test idmaps.wake_panel_ids_along_casing_wake_interface == [9, 10, 11]
    @test idmaps.wake_panel_ids_along_centerbody_wake_interface == [1, 2]

    @test idmaps.centerbody_panel_ids_along_centerbody_wake_interface == [3, 4]
    @test idmaps.duct_panel_ids_along_centerbody_wake_interface == [3, 2, 1]

    @test idmaps.id_of_first_casing_panel_aft_of_each_rotor == [3]
    @test idmaps.id_of_first_centerbody_panel_aft_of_each_rotor == [3]

    # hub and duct have same TE axial coordinate
    npanels = [3, 1] #npanels
    ncenterbody_inlet = 2 # paneling_constants.ncenterbody_inlet
    nwake_sheets = 3 # paneling_constants.nwake_sheets
    dte_minus_cbte = 0 # paneling_constants.dte_minus_cbte
    wake_nodemape = [
        1 2 3 4 6 7 8 9 11 12 13 14
        2 3 4 5 7 8 9 10 12 13 14 15
    ]# wake_vortex_panels.nodemap
    wake_endnodeidxs = [1 6 11; 5 10 15] # wake_vortex_panels.endnodeidxs
    nwp = 12 # pd.nwp
    nwsp = 4 # pd.nwsp
    nbn = 12 # pd.nbn
    rotor_indices_in_wake = [1] # rotor_indices_in_wake
    nrotor = 1 # pd.nrotor

    idmaps = dt.set_index_maps(
        npanels,
        ncenterbody_inlet,
        nwake_sheets,
        dte_minus_cbte,
        wake_nodemap,
        wake_endnodeidxs,
        nwp,
        nwsp,
        nbn,
        rotor_indices_in_wake,
        nrotor,
    )

    @test idmaps.wake_node_sheet_be_map[:, 1] ==
        [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3]
    @test all(idmaps.wake_node_sheet_be_map[:, 2] .== 1)
    @test idmaps.wake_panel_sheet_be_map[:, 1] == [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
    @test all(idmaps.wake_panel_sheet_be_map[:, 2] .== 1)

    @test idmaps.wake_node_ids_along_casing_wake_interface == [11, 12, 13, 14]
    @test idmaps.wake_node_ids_along_centerbody_wake_interface == [1, 2, 3, 4]

    @test idmaps.wake_panel_ids_along_casing_wake_interface == [9, 10, 11]
    @test idmaps.wake_panel_ids_along_centerbody_wake_interface == [1, 2, 3]

    @test idmaps.centerbody_panel_ids_along_centerbody_wake_interface == [3, 4, 5]
    @test idmaps.duct_panel_ids_along_centerbody_wake_interface == [3, 2, 1]

    @test idmaps.id_of_first_casing_panel_aft_of_each_rotor == [3]
    @test idmaps.id_of_first_centerbody_panel_aft_of_each_rotor == [3]
end

# @testset "Velocity Initialization" begin
#     dt.initialize_velocities(
#         operating_point, blade_elements, linsys, ivr, ivw, body_totnodes, wake_node_sheet_be_map
#     )
# end

