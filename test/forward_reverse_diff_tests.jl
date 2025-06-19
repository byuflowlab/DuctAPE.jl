println("\nAUTOMATIC DERIVATION TESTS")

@testset "Compare Forward and Revese Diff" begin

    #TODO: things to test
    #=
    # possible design variable inputs:
    # - duct geometry
    #   - inlet area
    #   - outlet area
    #   - duct length
    #   - duct thickness
    #   - hub length
    #   - hub radius
    # - rotor geometry
    #   - twist
    #   - chord/solidity
    #   - tip radius
    #   - rotor position
    # - operating conditions
    #   - Vinf
    #   - Omega
    #
    # Final outputs
    # - Rotor CT, CP, CQ
    # - Body thrust
    # - efficiencies
    #
    #
    # TODO: start with post-processing and work backward from there so you don't have to look at a ton of outputs.
    # TODO: do test with DFDC example geometry for now.
    =#

end

@testset "temp" begin
    datapath = "data/dfdc_init_iter1/"
    include(datapath * "element_indices.jl")

    # get ranges of element idices
    nidx = dfdc_eid[:, 2:3]
    pidx = dfdc_eid[:, 4:5]
    nidranges = nidx[:, 2] .- nidx[:, 1] .+ 1
    pidranges = pidx[:, 2] .- pidx[:, 1] .+ 1

    # - read in geometry file to help find indices - #
    include(datapath * "dfdc_geometry.jl")

    # get wake sheet length
    num_wake_z_panels = num_wake_z_nodes - 1

    # get number of nodes and panels on bodies interface with wake
    nhub_interface_nodes = num_wake_z_nodes - nidranges[4]
    nduct_interface_nodes = num_wake_z_nodes - nidranges[end]
    nhub_interface_panels = num_wake_z_panels - pidranges[4]
    nduct_interface_panels = num_wake_z_panels - pidranges[end]

    # - read in functions to reformat - #
    include(datapath*"reformatting_functions.jl")

    include(datapath * "iter1_relaxed_bgam.jl")
    Gamr1 = reformat_circulation(bgamr1, B)

    include(datapath*"iter2_sigr.jl")

    include(datapath*"iter1_relaxed_gamw.jl")
    gamw1 = reformat_gamw(relaxed_GTH, nidx, nhub_interface_nodes, nduct_interface_nodes)

    input_vector = [Gamr1; sigr2; -gamw1]
    @time "initial run" begin
        output = dt_post_wrapper(input_vector)
    end
    @time "ForwardDiff jacobian" begin
        fordiff_j = ForwardDiff.jacobian(dt_post_wrapper, input_vector)
    end
end

"""
"""
function dt_post_wrapper(input_vector)
    TF = eltype(input_vector)

    datapath = "test/data/dfdc_init_iter1/"
    include(datapath * "ductape_parameters.jl")
    include(datapath * "ductape_formatted_dfdc_geometry.jl")

    # generate inputs
    inputs = dt.precomputed_inputs(
        system_geometry,
        rotor, #vector of named tuples
        freestream,
        reference_parameters;
        debug=false,
    )

    nbe = rotor[1].num_wake_sheets - 1

    Gamr = zeros(TF, nbe, 1) .= input_vector[1:nbe]
    sigr = zeros(TF, nbe + 1, 1) .= input_vector[(nbe + 1):(nbe * 2 + 1)]
    gamw = input_vector[(nbe * 2 + 2):end]

    # calculate outputs
    outs = dt.post_process(Gamr, sigr, gamw, inputs; write_outputs=false)

    # put relevant outputs in a vector
    outputs = [
        outs.totals.CT
        outs.totals.total_efficiency
        outs.rotors.blade_normal_force_per_unit_span
        outs.rotors.blade_tangential_force_per_unit_span
    ]

    return outputs
end
