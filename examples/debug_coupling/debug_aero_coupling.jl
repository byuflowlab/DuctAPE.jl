#=

Manually run through everything and plot everything to check for bugs

=#

#---------------------------------#
#             Includes            #
#---------------------------------#
using DuctTAPE
const dt = DuctTAPE

# CCBlade used for it's airfoils function objects here.
using CCBlade
const ccb = CCBlade

using FLOWMath
const fm = FLOWMath

include("../../plots_default.jl")
include("../../test/compare_objects.jl")

#plot flags
plot_inputs = true
plot_initialization = true
plot_analysis = true

######################################################################
#                                                                    #
#                       CHECK SEPARATE PIECES                        #
#                                                                    #
######################################################################
#---------------------------------#
#          user inputs            #
#---------------------------------#
include("debug_user_inputs.jl")

#---------------------------------#
#        constant inputs          #
#---------------------------------#
if !plot_inputs

    # initialize various inputs used in analysis
    inputs = dt.precomputed_inputs(
        duct_coordinates,
        # hub_coordinates,
        nothing,
        paneling_constants,
        rotor_parameters,
        freestream;
        # finterp=FLOWMath.linear,
    )

else
    include("debug_input_precomputation.jl")

    # check that inputs were calculated correctly
    @assert all(
        compare_namedtuples(
            inputs,
            dt.precomputed_inputs(
                duct_coordinates,
                # hub_coordinates,
                nothing,
                paneling_constants,
                rotor_parameters,
                freestream;
                # finterp=FLOWMath.linear,
            ),
        ),
    )
end

#---------------------------------#
#        Initialize States        #
#---------------------------------#
if !plot_initialization
    states = dt.initialize_states(inputs)
else
    include("debug_state_initialization.jl")

    @assert all(states .== dt.initialize_states(inputs))
end

#---------------------------------#
#           Run Analysis          #
#---------------------------------#

if !plot_analysis
    strengths = dt.analyze_propulsor(
        duct_coordinates,
        # hub_coordinates,
        nothing,
        paneling_constants,
        rotor_parameters,
        freestream;
        tol=1e-8,
        maxiter=100,
    )
else
    if !plot_initialization
        gamb, gamw, Gamr, sigr = dt.extract_state_variables(states, inputs)
        dp = inputs.body_panels[1].panel_center[:, 1]
        # hp = inputs.body_panels[2].panel_center[:, 1]
        gamd = 1.0 .- (gamb[1:length(dp)] ./ Vinf) .^ 2
        # gamh = 1.0 .- (gamb[(length(dp) + 1):end]./Vinf).^2

        pb = plot(
            dp,
            gamd;
            yflip=true,
            xlabel="x",
            ylabel="cp",
            label="initial duct surface pressure",
        )
        # plot!(pb, hp, gamh ; label="initial hub surface pressure")
        pG = plot(
            Gamr, inputs.rotor_panel_centers; xlabel=L"\Gamma", ylabel="r", label="initial"
        )

        ps = plot(
            sigr, inputs.rotor_panel_centers; xlabel=L"\sigma", ylabel="r", label="initial"
        )

        pw = plot(
            gamw,
            inputs.rotor_panel_edges;
            xlabel=L"\gamma_\theta",
            ylabel="r",
            label="initial",
        )
    end
    include("debug_analysis.jl")
end

strengths = dt.analyze_propulsor(
    duct_coordinates,
    # hub_coordinates,
    nothing,
    paneling_constants,
    rotor_parameters,
    freestream;
)
gamb, gamw, Gamr, sigr = dt.extract_state_variables(states, inputs)
gamd = 1.0 .- (gamb[1:length(dp)] ./ Vinf) .^ 2
# gamh = 1.0 .- (gamb[(length(dp) + 1):end]./Vinf).^2

plot!(pb, dp, gamd; xlabel="x", ylabel="cp", label="converged duct surface pressure")
# plot!(pb, hp, gamh ; label="iter hub surface pressure")

## -- check rotor circulation and source initial strengths -- ##
plot!(pG, Gamr, inputs.rotor_panel_centers; xlabel=L"\Gamma", ylabel="r", label="converged")

plot!(ps, sigr, inputs.rotor_panel_centers; xlabel=L"\sigma", ylabel="r", label="converged")

plot!(
    pw,
    gamw,
    inputs.rotor_panel_edges;
    xlabel=L"\gamma_\theta",
    ylabel="r",
    label="converged",
)

savefig(pb, "examples/debug_coupling/body-strengths-converged.pdf")
savefig(pG, "examples/debug_coupling/circulation-converged.pdf")
savefig(ps, "examples/debug_coupling/sources-converged.pdf")
savefig(pw, "examples/debug_coupling/wake-strengths-converged.pdf")

