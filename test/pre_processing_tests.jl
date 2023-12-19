#=
Tests for all the components of precomputed_inputs function
=#
println("\nPRECOMPUTED INPUTS TESTS")

#simple geometry to work with:
# include("data/basic_two_rotor_for_test.jl")

@testset "Wake Discretization" begin
    # get input data
    include("data/basic_two_rotor_for_test.jl")

    zwake, rotor_indices_in_wake, ductTE_index, hubTE_index = dt.discretize_wake(
        duct_coordinates,
        hub_coordinates,
        rotorstator_parameters.rotorzloc,
        paneling_constants.wake_length,
        paneling_constants.npanels,
    )

    @test zwake == [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
    @test rotor_indices_in_wake == [1, 3]
    @test ductTE_index == 4
    @test hubTE_index == 4
end

@testset "Body Repanleing" begin
    # get input data
    include("data/basic_two_rotor_for_test.jl")

    rp_duct_coordinates, rp_hub_coordinates = dt.repanel_bodies(
        duct_coordinates,
        hub_coordinates,
        zwake,
        paneling_constants.nhub_inlet,
        paneling_constants.nduct_inlet;
        finterp=fm.linear,
    )

    @test size(rp_duct_coordinates, 1) < size(rp_duct_coordinates, 2)
    @test size(rp_hub_coordinates, 1) < size(rp_hub_coordinates, 2)
    @test size(rp_duct_coordinates, 2) == 2 * size(duct_coordinates, 1) - 1
    @test size(rp_hub_coordinates, 2) == 2 * size(hub_coordinates, 1) - 1

    @test rp_duct_coordinates == [
        1.0 0.75 0.5 0.25 0.0 0.25 0.5 0.75 1.0
        2.0 1.75 1.5 1.75 2.0 2.25 2.5 2.25 2.0
    ]
    @test rp_hub_coordinates == [
        0.0 0.25 0.5 0.75 1.0
        0.0 0.25 0.5 0.25 0.0
    ]
end

@testset "Duct Auto Placement" begin
    # get input data
    include("data/basic_two_rotor_for_test.jl")

    zwake, rotor_indices_in_wake, ductTE_index, hubTE_index = dt.discretize_wake(
        duct_coordinates,
        hub_coordinates,
        rotorstator_parameters.rotorzloc,
        paneling_constants.wake_length,
        paneling_constants.npanels,
    )

    rp_duct_coordinates, rp_hub_coordinates = dt.repanel_bodies(
        duct_coordinates,
        hub_coordinates,
        zwake,
        paneling_constants.nhub_inlet,
        paneling_constants.nduct_inlet;
        finterp=fm.linear,
    )

    rpb4 = copy(rp_duct_coordinates)

    dt.place_duct!(
        rp_duct_coordinates,
        rotorstator_parameters[1].Rtip,
        rotorstator_parameters[1].rotorzloc,
        rotorstator_parameters[1].tip_gap,
    )

    @test rp_duct_coordinates[1, :] == rpb4[1, :]
    @test rp_duct_coordinates[2, :] == rpb4[2, :] .- 0.75
end

@testset "Blade Extremity Calculation" begin
    # get input data
    include("data/basic_two_rotor_for_test.jl")

    zwake, rotor_indices_in_wake, ductTE_index, hubTE_index = dt.discretize_wake(
        duct_coordinates,
        hub_coordinates,
        rotorstator_parameters.rotorzloc,
        paneling_constants.wake_length,
        paneling_constants.npanels,
    )

    rp_duct_coordinates, rp_hub_coordinates = dt.repanel_bodies(
        duct_coordinates,
        hub_coordinates,
        zwake,
        paneling_constants.nhub_inlet,
        paneling_constants.nduct_inlet;
        finterp=fm.linear,
    )

    dt.place_duct!(
        rp_duct_coordinates,
        rotorstator_parameters[1].Rtip,
        rotorstator_parameters[1].rotorzloc,
        rotorstator_parameters[1].tip_gap,
    )

    Rtips, Rhubs = dt.get_blade_ends_from_body_geometry(
        rp_duct_coordinates,
        rp_hub_coordinates,
        rotorstator_parameters.tip_gap,
        rotorstator_parameters.rotorzloc,
    )

    @test all(Rtips .== 1.0)
    @test all(Rhubs .== 0.25)
end

@testset "Wake Initialization" begin
    # get input data
    include("data/basic_two_rotor_for_test.jl")

    zwake, rotor_indices_in_wake, ductTE_index, hubTE_index = dt.discretize_wake(
        duct_coordinates,
        hub_coordinates,
        rotorstator_parameters.rotorzloc,
        paneling_constants.wake_length,
        paneling_constants.npanels,
    )

    rp_duct_coordinates, rp_hub_coordinates = dt.repanel_bodies(
        duct_coordinates,
        hub_coordinates,
        zwake,
        paneling_constants.nhub_inlet,
        paneling_constants.nduct_inlet;
        finterp=fm.linear,
    )

    dt.place_duct!(
        rp_duct_coordinates,
        rotorstator_parameters[1].Rtip,
        rotorstator_parameters[1].rotorzloc,
        rotorstator_parameters[1].tip_gap,
    )

    Rtips, Rhubs = dt.get_blade_ends_from_body_geometry(
        rp_duct_coordinates,
        rp_hub_coordinates,
        rotorstator_parameters.tip_gap,
        rotorstator_parameters.rotorzloc,
    )

    rwake = range(Rhubs[1], Rtips[1]; length=paneling_constants.nwake_sheets)

    # Check grid initialization
    grid = dt.initialize_wake_grid(rp_duct_coordinates, rp_hub_coordinates, zwake, rwake)

    @test grid[:, :, 1] == [
        0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0
        0.25 0.5 0.25 0.0 0.0 0.0 0.0 0.0
    ]
    @test isapprox(
        grid[:, :, 2],
        [
            0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0
            0.625 0.599479 0.625 0.73951 0.73951 0.73951 0.73951 0.73951
        ],
        atol=1e-6,
    )
    @test grid[:, :, 3] == [
        0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0
        1.0 0.75 1.0 1.25 1.25 1.25 1.25 1.25
    ]

    # check wake paneling
    wake_vortex_panels = dt.generate_wake_panels(
        grid[1, :, 1:length(rwake)], grid[2, :, 1:length(rwake)]
    )

    # # TODO: add any checks for things that get used.
    # @test wake_vortex_panels.node[1,:] == reduce(vcat,grid[1,:,:])
    # @test wake_vortex_panels.node[2,:] == reduce(vcat,grid[2,:,:])
    # @test wake_vortex_panels.endpanelidxs == [1 8 15; 7 14 21]
    # @test wake_vortex_panels.endnodeidxs == [1 9 17; 8 16 24]

    # check wakeK calcualtion
    wakeK = dt.get_wake_k(wake_vortex_panels)

    @test isapprox(
        wakeK,
        [
            -0.20264236728467555
            -0.05066059182116889
            -0.20264236728467555
            0.0
            0.0
            0.0
            0.0
            0.0
            -0.03242277876554809
            -0.03524215083211749
            -0.03242277876554809
            -0.0231591276896772
            -0.0231591276896772
            -0.0231591276896772
            -0.0231591276896772
            -0.0231591276896772
            -0.012665147955292222
            -0.02251581858718617
            -0.012665147955292222
            -0.008105694691387022
            -0.008105694691387022
            -0.008105694691387022
            -0.008105694691387022
            -0.008105694691387022
        ],
    )

    # just make sure it converges
    dt.relax_grid!(grid; max_iterations=100, tol=1e-9, verbose=false)
end

@testset "Rotor Initialization" begin
    # get input data
    include("data/basic_two_rotor_for_test.jl")

    zwake, rotor_indices_in_wake, ductTE_index, hubTE_index = dt.discretize_wake(
        duct_coordinates,
        hub_coordinates,
        rotorstator_parameters.rotorzloc,
        paneling_constants.wake_length,
        paneling_constants.npanels,
    )

    rp_duct_coordinates, rp_hub_coordinates = dt.repanel_bodies(
        duct_coordinates,
        hub_coordinates,
        zwake,
        paneling_constants.nhub_inlet,
        paneling_constants.nduct_inlet;
        finterp=fm.linear,
    )

    dt.place_duct!(
        rp_duct_coordinates,
        rotorstator_parameters[1].Rtip,
        rotorstator_parameters[1].rotorzloc,
        rotorstator_parameters[1].tip_gap,
    )

    Rtips, Rhubs = dt.get_blade_ends_from_body_geometry(
        rp_duct_coordinates,
        rp_hub_coordinates,
        rotorstator_parameters.tip_gap,
        rotorstator_parameters.rotorzloc,
    )

    num_rotors = length(rotorstator_parameters)

    # rotor source panel objects
    rotor_source_panels = [
        dt.generate_rotor_panels(
            rotorstator_parameters[i].rotorzloc,
            grid[2, rotor_indices_in_wake[i], 1:length(rwake)],
        ) for i in 1:num_rotors
    ]

    # rotor blade element objects
    blade_elements = [
        dt.generate_blade_elements(
            rotorstator_parameters[i].B,
            rotorstator_parameters[i].Omega,
            rotorstator_parameters[i].rotorzloc,
            rotorstator_parameters[i].r,
            rotorstator_parameters[i].chords,
            rotorstator_parameters[i].twists,
            rotorstator_parameters[i].airfoils,
            Rtips[i],
            Rhubs[i],
            rotor_source_panels[i].controlpoint[2, :];
            fliplift=rotorstator_parameters[i].fliplift,
        ) for i in 1:num_rotors
    ]

    @test blade_elements[1].inner_fraction == [0.75, 0.25]
    @test blade_elements[1].stagger == [70, 70] * pi / 180
    @test blade_elements[1].rbe == [0.4375, 0.8125]
    @test blade_elements[1].B == 2
    @test blade_elements[1].chords == [0.1, 0.1]
    @test blade_elements[1].twists == [20, 20] * pi / 180
    @test blade_elements[1].solidity == [
        0.07275654541343787
        0.039176601376466544
    ]
    @test blade_elements[1].rotorzloc == 0.25
    @test blade_elements[1].fliplift == false
    @test blade_elements[1].Omega == rotorstator_parameters[1].Omega
    @test blade_elements[1].Rhub == rotorstator_parameters[1].Rhub
    @test blade_elements[1].Rtip == rotorstator_parameters[1].Rtip

    @test blade_elements[2].inner_fraction == [
        0.7334602398026477
        0.23346023980264796
    ]
    @test blade_elements[2].stagger == [70, 70] * pi / 180
    @test blade_elements[2].rbe == [
        0.43336505995066193
        0.808365059950662
    ]
    @test blade_elements[2].B == 4
    @test blade_elements[2].chords == [0.1, 0.1]
    @test blade_elements[2].twists == [20, 20] * pi / 180
    @test blade_elements[2].solidity == [
        0.1469014997286721
        0.07875399419247996
    ]
    @test blade_elements[2].rotorzloc == 0.75
    @test blade_elements[2].fliplift == false
    @test blade_elements[2].Omega == rotorstator_parameters[2].Omega
    @test blade_elements[2].Rhub == rotorstator_parameters[2].Rhub
    @test blade_elements[2].Rtip == rotorstator_parameters[2].Rtip
end

@testset "Dimension, Index, and Double Checks" begin
    # get input data
    include("data/basic_two_rotor_for_test.jl")

    inputs = dt.precomputed_inputs(
        duct_coordinates,
        hub_coordinates,
        paneling_constants,
        rotorstator_parameters, #vector of named tuples
        freestream,
        reference_parameters;
        finterp=fm.akima,
        autoshiftduct=true,
    )

    # - Dimension Checks - #
    # check all the induced velocity matrix dimensions to make sure they're oriented correctly for later access.

    # - Index (Bookkeeping) Checks - #
    # check rotor indices on hub and duct
    # check rotorwakeid
    # check hubwakeinterfaceid
    @test inputs.hubwakeinterfaceid == 1:3
    # check ductwakeinterfaceid
    @test inputs.ductwakeinterfaceid == 15:17

    @test inputs.num_wake_x_panels == 7

    # - Double Check Everything - #
    # want to make sure it works in the full function as well as just the pieces
    @test inputs.zwake == [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
    @test inputs.rotor_indices_in_wake == [1, 3]
    @test inputs.ductTE_index == 4
    @test inputs.hubTE_index == 4
    @test inputs.isapprox(
        wakeK,
        [
            -0.20264236728467555
            -0.05066059182116889
            -0.20264236728467555
            0.0
            0.0
            0.0
            0.0
            0.0
            -0.03242277876554809
            -0.03524215083211749
            -0.03242277876554809
            -0.0231591276896772
            -0.0231591276896772
            -0.0231591276896772
            -0.0231591276896772
            -0.0231591276896772
            -0.012665147955292222
            -0.02251581858718617
            -0.012665147955292222
            -0.008105694691387022
            -0.008105694691387022
            -0.008105694691387022
            -0.008105694691387022
            -0.008105694691387022
        ],
    )
    @test inputs.blade_elements[1].inner_fraction == [0.75, 0.25]
    @test inputs.blade_elements[1].stagger == [70, 70] * pi / 180
    @test inputs.blade_elements[1].rbe == [0.4375, 0.8125]
    @test inputs.blade_elements[1].B == 2
    @test inputs.blade_elements[1].chords == [0.1, 0.1]
    @test inputs.blade_elements[1].twists == [20, 20] * pi / 180
    @test inputs.blade_elements[1].solidity == [
        0.07275654541343787
        0.039176601376466544
    ]
    @test inputs.blade_elements[1].rotorzloc == 0.25
    @test inputs.blade_elements[1].fliplift == false
    @test inputs.blade_elements[1].Omega == rotorstator_parameters[1].Omega
    @test inputs.blade_elements[1].Rhub == rotorstator_parameters[1].Rhub
    @test inputs.blade_elements[1].Rtip == rotorstator_parameters[1].Rtip

    @test inputs.blade_elements[2].inner_fraction == [
        0.7334602398026477
        0.23346023980264796
    ]
    @test inputs.blade_elements[2].stagger == [70, 70] * pi / 180
    @test inputs.blade_elements[2].rbe == [
        0.43336505995066193
        0.808365059950662
    ]
    @test inputs.blade_elements[2].B == 4
    @test inputs.blade_elements[2].chords == [0.1, 0.1]
end

