#---------------------------------#
#              SETUP              #
#---------------------------------#

project_dir = dirname(dirname(@__FILE__))
if project_dir == ""
    project_dir = "."
end

savepath = joinpath(project_dir, "validation_scripts/", "outputs/")

using DelimitedFiles
using Statistics
using FLOWMath
const fm = FLOWMath
using DuctTAPE
const dt = DuctTAPE

# plotting defaults
include("../plots_default_new.jl")
pyplot()

#---------------------------------#
#            Initialize           #
#---------------------------------#
# - Read in Parameters - #
include("ducttape_parameters.jl")

# initialize various inputs used in analysis
inputs = dt.precomputed_inputs(
    duct_coordinates,
    hub_coordinates,
    paneling_constants,
    rotorstator_parameters,
    freestream,
    reference_parameters;
)

# TODO: run state initialization manually and get the flow angles at the very beginning.

# calculate initial guess for state variables
states = dt.initialize_states(inputs)

# extract states
mub1, gamw1, Gamr1, sigr1 = dt.extract_state_variables(states, inputs)
mub0 = similar(mub1) .= 0.0
gamw0 = similar(gamw1) .= 0.0
Gamr0 = similar(Gamr1) .= 0.0
sigr0 = similar(sigr1) .= 0.0

#---------------------------------#
#            Get Angles           #
#---------------------------------#

# Get velocities at rotor
vx_rotor, vr_rotor, vtheta_rotor, Wx_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor = dt.calculate_rotor_velocities(
    Gamr0, gamw0, sigr0, mub0, inputs
)

# get the inflow and attack angles
rotor_inflow, rotor_aoa = dt.calculate_inflow_angles(
    Wm_rotor[:, 1], Wtheta_rotor[:, 1], inputs.blade_elements[1].twists
)

stator_inflow, stator_aoa = dt.calculate_inflow_angles(
    Wm_rotor[:, 2], Wtheta_rotor[:, 2], inputs.blade_elements[2].twists
)

#---------------------------------#
#              PLOT               #
#---------------------------------#
rotor_r = inputs.rotor_panel_centers[:, 1] ./ inputs.reference_parameters.Rref
stator_r = inputs.rotor_panel_centers[:, 2] ./ inputs.reference_parameters.Rref

pr = plot(;
    title="Rotor Angles, rotor RPM=$RPM_rotor",
    xlabel="Angles (degrees)",
    ylabel="Normalized Radial Position",
)
plot!(pr, inputs.blade_elements[1].twists * 180 / pi, rotor_r; label="Twist")
plot!(pr, rotor_inflow * 180 / pi, rotor_r; label="Inflow")
plot!(pr, rotor_aoa * 180 / pi, rotor_r; label="Attack")
savefig(pr, savepath * "rotor_initial_angles$RPM_rotor.pdf")
savefig(pr, savepath * "rotor_initial_angles$RPM_rotor.png")

ps = plot(;
    title="Stator Angles, rotor RPM=$RPM_rotor",
    xlabel="Angles (degrees)",
    ylabel="Normalized Radial Position",
)
plot!(ps, inputs.blade_elements[2].twists * 180 / pi, stator_r; label="Twist")
plot!(ps, stator_inflow * 180 / pi, stator_r; label="Inflow")
plot!(ps, stator_aoa * 180 / pi, stator_r; label="Attack")
savefig(ps, savepath * "stator_initial_angles$RPM_rotor.pdf")
savefig(ps, savepath * "stator_initial_angles$RPM_rotor.png")
