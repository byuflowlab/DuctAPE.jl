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

@testset "Body Repaneling" begin
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

    rwake = range(Rhubs[1], Rtips[1]; length=paneling_constants.nwake_sheets)

    # Check grid initialization
    grid = dt.initialize_wake_grid(rp_duct_coordinates, rp_hub_coordinates, zwake, rwake)

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
    @test isapprox(
        blade_elements[1].solidity,
        [
            0.07275654541343787
            0.039176601376466544
        ],
        atol=1e-6,
    )
    @test blade_elements[1].rotorzloc == 0.25
    @test blade_elements[1].fliplift == false
    @test blade_elements[1].Omega == rotorstator_parameters[1].Omega
    @test blade_elements[1].Rhub == rotorstator_parameters[1].Rhub
    @test blade_elements[1].Rtip == rotorstator_parameters[1].Rtip

    @test blade_elements[2].inner_fraction == [0.75, 0.25]
    @test blade_elements[2].stagger == [70, 70] * pi / 180
    @test blade_elements[2].rbe == [0.4375, 0.8125]
    @test blade_elements[2].B == 4
    @test blade_elements[2].chords == [0.1, 0.1]
    @test blade_elements[2].twists == [20, 20] * pi / 180
    @test isapprox(
        blade_elements[2].solidity, [0.14551309082687575, 0.07835320275293309], atol=1e-6
    )
    @test blade_elements[2].rotorzloc == 0.75
    @test blade_elements[2].fliplift == false
    @test blade_elements[2].Omega == rotorstator_parameters[2].Omega
    @test blade_elements[2].Rhub == rotorstator_parameters[2].Rhub
    @test blade_elements[2].Rtip == rotorstator_parameters[2].Rtip
end

#TODO: geometry generation function test
@testset "Geometry Generation" begin
    include("data/basic_two_rotor_for_test.jl")

    sg = dt.generate_geometry(
        duct_coordinates,
        hub_coordinates,
        paneling_constants,
        rotorstator_parameters; #vector of named tuples
        finterp=fm.linear,
        autoshiftduct=false,
    )

    @test sg.duct_coordinates == [
        1.0 0.75 0.5 0.25 0.0 0.25 0.5 0.75 1.0
        2.0 1.75 1.5 1.75 2.0 2.25 2.5 2.25 2.0
    ]

    @test sg.hub_coordinates == [0.0 0.25 0.5 0.75 1.0; 0.0 0.25 0.5 0.25 0.0]

    @test sg.Rhubs == [0.25, 0.25]
    @test sg.Rtips == [1.75, 1.75]
    @test sg.noduct == false
    @test sg.nohub == false
    @test sg.rotoronly == false
    @test sg.rwake == sg.rpe
    @test sg.rpe == 0.25:0.75:1.75
    @test sg.zwake == [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
    @test sg.rotor_indices_in_wake == [1, 3]

    testgrid = similar(sg.grid) .= 0.0
    testgrid[:, :, 1] = [
        0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0
        0.25 0.5 0.25 0.0 0.0 0.0 0.0 0.0
    ]
    testgrid[:, :, 2] = [
        0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0
        1.0 0.9861636413819727 1.0008941893344476 1.0314740183248352 1.056241323544226 1.0718388705767528 1.080459189206699 1.0833330706503477
    ]
    testgrid[:, :, 3] = [
        0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0
        1.75 1.5 1.75 2.0 2.0 2.0 2.0 2.0
    ]

    @test isapprox(sg.grid, testgrid, atol=1e-12)
end

@testset "Dimension and Indexing Checks" begin
    # get input data
    include("data/basic_two_rotor_for_test.jl")

    inputs = dt.precomputed_inputs(
        duct_coordinates,
        hub_coordinates,
        paneling_constants,
        rotorstator_parameters, #vector of named tuples
        freestream,
        reference_parameters;
        finterp=fm.linear,
        autoshiftduct=true,
    )

    # - Dimension Checks - #
    # check all the induced velocity matrix dimensions to make sure they're oriented correctly for later access.
    sysN = inputs.body_vortex_panels.totnode + 2
    npan = inputs.body_vortex_panels.totpanel
    nnode = inputs.body_vortex_panels.totnode
    nrnode = inputs.rotor_source_panels.totnode
    nrpan = inputs.rotor_source_panels.totpanel
    nrotor = length(nrnode)
    nwpan = inputs.wake_vortex_panels.totpanel
    nwnode = inputs.wake_vortex_panels.totnode
    @test size(inputs.converged) == (1,)
    @test size(inputs.b_bf) == (sysN,)
    @test size(inputs.blade_elements) == (2,)
    @test size(inputs.rotor_source_panels) == (2,)
    @test size(inputs.A_bb) == (sysN, sysN)
    @test size(inputs.vz_rb) == (2,)
    @test size(inputs.vz_rb[1]) == (2, nnode)
    @test size(inputs.vr_rb) == (2,)
    @test size(inputs.vr_rb[1]) == (2, nnode)
    @test size(inputs.A_br) == (npan, nrnode[1], nrotor)
    @test size(inputs.vz_rr) == (2, 2)
    @test size(inputs.vz_rr[1]) == (2, nrnode[1])
    @test size(inputs.vr_rr) == (2, 2)
    @test size(inputs.vr_rr[1]) == (2, nrnode[1])
    @test size(inputs.A_bw) == (npan, nwnode)
    @test size(inputs.vz_rw) == (2,)
    @test size(inputs.vz_rw[1]) == (2, nwnode)
    @test size(inputs.vr_rw) == (2,)
    @test size(inputs.vr_rw[1]) == (2, nwnode)

    # - Index (Bookkeeping) Checks - #
    # check rotor indices on hub and duct
    # check rotorwakeid
    @test size(inputs.rotorwakeid, 2) == 2
    @test size(inputs.rotorwakeid, 1) == inputs.wake_vortex_panels.totnode
    # check hubwakeinterfaceid
    @test inputs.hubwakeinterfaceid == 1:3
    # check ductwakeinterfaceid
    @test inputs.ductwakeinterfaceid == 15:17

    @test inputs.num_wake_z_panels == 7
    @test all(
        inputs.rotorwakeid .== [
            1 1
            1 1
            1 2
            1 2
            1 2
            1 2
            1 2
            1 2
            2 1
            2 1
            2 2
            2 2
            2 2
            2 2
            2 2
            2 2
            3 1
            3 1
            3 2
            3 2
            3 2
            3 2
            3 2
            3 2
        ],
    )
end

@testset "Rotor/Wake Aero Initialization" begin

    # - Set up Rotor raw geometry - #
    # set up ccblade example
    Rtip = 10 / 2.0 * 0.0254  # inches to meters
    Rhub = 0.10 * Rtip
    B = [2; 4]  # number of blades

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
        1.00 0.041 8.99
    ]

    r = propgeom[:, 1] * Rtip
    chords = propgeom[:, 2] * Rtip
    twists = propgeom[:, 3] * pi / 180
    af = dt.c4b.AlphaAF("data/naca4412.dat")

    Vinf = 5.0
    Omega = [5400 * pi / 30; 0.0]  # convert to rad/s
    rhoinf = 1.225
    muinf = 1.81e-5
    asound = 343.0

    rotorzloc = [0.0, 0.5]
    nwake_sheets = length(r)
    num_rotors = length(Omega)
    fliplift = [false, false]

    # - put things in terms of blade elements and such - #
    rotorstator_parameters = [
        (;
            rotorzloc=rotorzloc[i],
            nwake_sheets,
            r=r ./ Rtip, #non-dimensionalize
            chords=chords,
            twists=twists,
            airfoils=fill(af, length(r)),
            Rtip,
            Rhub,
            tip_gap=0.0,
            B=B[i],
            Omega=Omega[i],
            fliplift=fliplift[i],
        ) for i in 1:length(Omega)
    ]
    paneling_constants = (; npanels, nhub_inlet, nduct_inlet, wake_length, nwake_sheets)
    freestream = (; rhoinf, muinf, asound, Vinf)

    # define rotor panels and blade elements
    rotor_source_panels = [
        dt.generate_rotor_panels(
            rotorstator_parameters[i].rotorzloc, rotorstator_parameters[i].r * Rtip
        ) for i in 1:num_rotors
    ]

    rotor_panel_centers = reduce(hcat, [r.controlpoint[2, :] for r in rotor_source_panels])

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
            Rtip,
            Rhub,
            rotor_source_panels[i].controlpoint[2, :];
            fliplift=rotorstator_parameters[i].fliplift,
        ) for i in 1:num_rotors
    ]

    # wake grid
    zwake = [
        range(rotorzloc[1], rotorzloc[2]; step=0.125)[1:(end - 1)]
        range(rotorzloc[2], 1.5; step=0.125)
    ]
    rwake = rotor_source_panels[1].node[2, :]
    grid = zeros(2, length(zwake), length(rwake))
    grid[1, :, :] .= zwake
    grid[2, :, :] .= rwake'

    # generate wake sheet paneling
    wake_vortex_panels = dt.generate_wake_panels(grid[1, :, :], grid[2, :, :])

    # calculate radius dependent "constant" for wake strength calcualtion
    wakeK = dt.get_wake_k(wake_vortex_panels)

    # Go through the wake panels and determine the index of the aftmost rotor infront and the blade node from which the wake strength is defined.
    rotorwakepanelid = ones(Int, wake_vortex_panels.totpanel, 2)
    num_wake_x_panels = length(zwake) - 1
    for i in 1:(length(rwake))
        rotorwakepanelid[(1 + (i - 1) * num_wake_x_panels):(i * num_wake_x_panels), 1] .= i
    end
    for (i, wn) in enumerate(eachcol(wake_vortex_panels.controlpoint))
        rotorwakepanelid[i, 2] = findlast(x -> x <= wn[1], rotorzloc)
    end

    # Go through the wake panels and determine the index of the aftmost rotor infront and the blade node from which the wake strength is defined.
    rotorwakeid = ones(Int, wake_vortex_panels.totnode, 2)
    num_wake_x_nodes = length(zwake)
    for i in 1:(length(rwake))
        rotorwakeid[(1 + (i - 1) * num_wake_x_nodes):(i * num_wake_x_nodes), 1] .= i
    end
    for (i, wn) in enumerate(eachcol(wake_vortex_panels.node))
        rotorwakeid[i, 2] = findlast(x -> x <= wn[1], rotorzloc)
    end

    inputs = (;
        wakeK,
        rotorwakeid,
        rotorwakepanelid,
        freestream,
        blade_elements,
        rotor_panel_centers,
        wake_vortex_panels,
        ductwakeinterfaceid=nothing,
        hubwakeinterfaceid=nothing,
        num_wake_x_nodes,
        num_wake_x_panels,
        Vconv=[0.0],
    )

    # - initialize outputs - #
    sigr = zeros(size(reduce(hcat, rotorstator_parameters.r)))
    Gamr = similar(sigr, size(sigr, 1) - 1, size(sigr, 2)) .= 0.0
    gamw = zeros(wake_vortex_panels.totnode)

    # - run init function - #
    dt.initialize_rotorwake_aero!(Gamr, sigr, gamw, inputs)

    # - test - #
    @test all(Gamr[:, 1] .> 0.0)
    @test all(Gamr[:, 2] .< 0.0)
    @test all(sigr .!= 0.0)
    gamwcheck = reshape(gamw, num_wake_x_nodes, length(rwake))
    for i in 1:3
        @test all(gamwcheck[i, :] .== gamwcheck[i + 1, :])
    end
    for i in 6:(length(zwake) - 1)
        @test all(gamwcheck[i, :] .== gamwcheck[i + 1, :])
    end
    @test all(gamw .!= 0.0)
end
