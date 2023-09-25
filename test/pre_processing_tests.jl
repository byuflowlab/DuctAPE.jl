@testset "Book-keeping Checks" begin

# load basic 2-rotor geometry that can be hand counted
include("data/basic_two_rotor_for_test.jl")

# get inputs
inputs = dt.precomputed_inputs(
    duct_coordinates,
    hub_coordinates,
    paneling_constants,
    rotor_parameters,
    freestream,
    reference_parameters,
)

# check interface indexing and such
@test inputs.hubwakeinterfaceid == 1:3
@test inputs.ductwakeinterfaceid == 15:17
@test inputs.num_wake_x_panels == 7

# put in wake discritization function test instead of here. they are used inside the precomputed inputs function, but should not be passed out.
# @test rotor_indices_in_wake == [1; 3]
# @test hubTE_index = 4
# @test ductTE_index = 4
end
