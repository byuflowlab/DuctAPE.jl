println("\nSOLVER CHECKS")

@testset "Wake Vorticity Calculations" begin
    # - Single Rotor - #
    include("data/dfdc_init_iter1/reformat_wake_data.jl")
    gamw = similar(gamw_init) .= 0

    # generate inputs
    include("data/dfdc_init_iter1/ductape_formatted_dfdc_geometry.jl")
    include("data/dfdc_init_iter1/ductape_parameters.jl")
    inputs = dt.precomputed_inputs(
        system_geometry,
        rotor, #vector of named tuples
        freestream,
        reference_parameters;
        debug=false,
    )

    # calculate wake strengths
    # TODO: input here is not Wm_avg, but Wm_panel
    gamw_out, deltaGamma2, deltaH = dt.calculate_wake_vortex_strengths!(
        gamw, Gamr_init, Wm_avg_init, inputs; post=true
    )
end

@testset "" begin
end

