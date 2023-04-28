#---------------------------------#
#             Includes            #
#---------------------------------#
project_dir = dirname(dirname(@__FILE__))

using DuctTAPE
const dt = DuctTAPE
include("duct_only/run_duct.jl")

# CCBlade used for it's airfoils function objects here.
using CCBlade
const ccb = CCBlade
include("rotor_only/run_ccblade.jl")

using FLOWMath
const fm = FLOWMath

include("../plots_default.jl")

function full_clean(;
    tip_gap=[0.1],
    xrotor=[0.25],
    npanels_inlet=10,
    wake_length=1.0,
    discscale=1,
    duct_file=project_dir * "/test/data/naca_662-015.jl",
    airfoil_file=project_dir * "/test/data/naca4412.dat",
)
    #---------------------------------#
    #         ROTOR Geometry          #
    #---------------------------------#

    # Blade Tip Radius, in meters
    Rtip = 10 / 2.0 * 0.0254  # inches to meters

    # Blade Hub radius, in meters
    Rhub = 0.10 * Rtip

    # number of blades
    B = 2

    # x position of rotor
    # xr = 0.25

    # Blade section non-dimensional radial positions, chords lengths, and local twists angles in degrees
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
        # 1.00 0.041 8.99
    ]

    # extract non-dimensional radial positions
    rnondim = propgeom[:, 1]
    # Dimensionalize chords
    chords = propgeom[:, 2] * Rtip
    # convert twists to radians
    twists = propgeom[:, 3] * pi / 180

    # use a NACA 4412 airfoils
    airfoils = fill(ccb.AlphaAF(airfoil_file), length(rnondim))

    #---------------------------------#
    #       Operating Conditions      #
    #---------------------------------#

    #Vinf
    Vinf = 5.0

    # rotor rotation rate in rad/s
    Omega = 5400 * pi / 30  # convert from RPM to rad/s

    # freestream conditions
    rho = 1.225 #kg/m^3
    mu = 1.81e-5 # kg/(m⋅s)
    asound = 341.0 #m/s

    #---------------------------------#
    #      Define BODY Coordinates    #
    #---------------------------------#
    scale = 0.5

    # - Duct Coordinates - #
    include(duct_file)
    duct_coordinates = [x_duct r_duct] .* scale

    # - Hub Coordinates - #
    # include("../test/data/bodyofrevolutioncoords.jl")
    # hub_coordinates = [x_hub[1:(end - 1)] ./ 2.0 r_hub[1:(end - 1)] * Rhub / maximum(r_hub)]
    hub_coordinates = nothing

    #---------------------------------#
    #         Paneling Options        #
    #---------------------------------#

    nwake_sheets = 15

    # non-dimensional wake length
    # wake_length = 1.0

    #= paneling for which zero tip gap case converges
    npanels = [10; 25]
    nhub_inlet = 20
    nduct_inlet = 20
    wakeoutratio = 25 / 10
    outletinletratio = 0.5
    =#

    for xr in xrotor
        ductle = minimum(duct_coordinates[:, 1])
        ductte = maximum(duct_coordinates[:, 1])
        ductchord = maximum(duct_coordinates[:, 1]) - minimum(duct_coordinates[:, 1])
        outletinletratio = (ductte - xr) / (xr - ductle)

        nhub_inlet = round(Int, npanels_inlet * discscale)

        nduct_inlet = round(Int, npanels_inlet * discscale)

        nduct_outlet = round(Int, nduct_inlet * outletinletratio)

        nwake = round(Int, (nduct_inlet + nduct_outlet) * wake_length)

        npanels = [nduct_outlet, nwake]

        nducttot = (nduct_inlet + nduct_outlet) * 2
        npersheet = nduct_outlet + nwake

        ###### LOOP ######

        for i in 1:length(tip_gap)
            #---------------------------------#
            #          Define Inputs          #
            #---------------------------------#

            # Rotor Parameters
            rotor_parameters = [(;
                xrotor=xr,
                nwake_sheets,
                r=rnondim,
                chords,
                twists,
                airfoils,
                Rtip,
                Rhub,
                tip_gap=tip_gap[i],
                B,
                Omega,
            )]

            # Paneling Parameters
            paneling_constants = (;
                npanels, nhub_inlet, nduct_inlet, wake_length, nwake_sheets
            )

            # Freestream Parameters
            freestream = (; rho, mu, asound, Vinf)

            ##### ----- Run Analysis ----- #####
            strengths, inputs, initials, convergeflag = dt.analyze_propulsor(
                duct_coordinates,
                # hub_coordinates,
                nothing,
                paneling_constants,
                rotor_parameters,
                freestream;
                debug=true,
            )

            println("solution converged? ", convergeflag)
            if convergeflag
                convlabel = "Converged"
            else
                convlabel = "NOT converged"
            end

            #extract solution
            gamb, gamw, Gamr, sigr = dt.extract_state_variables(strengths, inputs)
            gamb_init, gamw_init, Gamr_init, sigr_init = dt.extract_state_variables(
                initials, inputs
            )

            #---------------------------------#
            #          Run Isolated           #
            #---------------------------------#

            ccbouts = run_ccblade(Vinf)

            ductouts = run_duct(inputs.duct_coordinates)

            #---------------------------------#
            #              PLOTS              #
            #---------------------------------#
            println("PLOTTING...\n")

            ##### ----- Plot rotor circulation distribution ----- #####
            # initialize plot
            pG = plot(;
                title="nbodypanel$(nducttot)_npersheet$(npersheet)_tipgap$(tip_gap[i])"
            )

            # plot solution
            plot!(
                pG,
                Gamr,
                inputs.rotor_panel_centers;
                xlabel=L"\Gamma",
                ylabel="r",
                label=convlabel,
            )

            # plot initials
            plot!(
                pG,
                Gamr_init,
                inputs.rotor_panel_centers;
                xlabel=L"\Gamma",
                ylabel="r",
                label="Initial",
                color=mycolors[1],
                linestyle=:dash,
            )

            # plot ccblade isolated rotor solution
            plot!(
                pG,
                ccbouts.circ,
                ccbouts.r;
                color=mycolors[2],
                label="CCBlade, isolated rotor",
            )

            #save
            savefig(
                pG, "examples/rotorcirculation-check_xrotor$(xr)_tipgap$(tip_gap[i]).pdf"
            )

            ##### ----- Plot duct surface velocity ----- #####
            #prepare outputs
            dp = inputs.body_panels[1].panel_center[:, 1]
            _, leidx = findmin(dp)
            gamd = 1.0 .- (gamb[1:length(dp)] ./ Vinf) .^ 2
            gamd = gamb[1:length(dp)] ./ Vinf
            #split into inner and outer surfaces
            dpinner = dp[1:leidx]
            dpouter = dp[(leidx + 1):end]
            gamdinner = gamd[1:leidx]
            gamdouter = gamd[(leidx + 1):end]

            # initialize plot
            pb = plot(;
                title="nbodypanel$(nducttot)_npersheet$(npersheet)_tipgap$(tip_gap[i])",
                xlabel="x",
                ylabel=L"\frac{V_s}{V_\infty}",
            )

            # plot solution
            plot!(pb, dpinner, gamdinner; label=convlabel * " inner surface, with rotor")

            plot!(pb, dpouter, gamdouter; label=convlabel * " outer surface, with rotor")

            ##prepare initials
            #gamd_init = 1.0 .- (gamb_init[1:length(dp)] ./ Vinf) .^ 2
            #gamd_init = gamb_init[1:length(dp)] ./ Vinf
            ##split into inner and outer surfaces
            #gamdinner_init = gamd_init[1:leidx]
            #gamdouter_init = gamd_init[(leidx + 1):end]

            ## plot solution
            #plot!(pb, dpinner, gamdinner_init; label="initial inner surface, with rotor")
            #plot!(pb, dpouter, gamdouter_init; label="initial outer surface, with rotor")

            # plot isolated duct
            # spit inner/outer
            _, iductle = findmin(ductouts.x)
            idxin = ductouts.x[1:iductle]
            idxout = ductouts.x[(iductle + 1):end]
            idvsin = ductouts.vs_duct[1:iductle]
            idvsout = ductouts.vs_duct[(iductle + 1):end]

            plot!(
                pb,
                idxin,
                idvsin;
                label="isolated duct inner surface (initial state)",
                color=mycolors[1],
                linestyle=:dash,
            )

            plot!(
                pb,
                idxout,
                idvsout;
                label="isolated duct outer surface (initial state)",
                color=mycolors[2],
                linestyle=:dash,
            )

            #plot rotor location
            plot!(
                pb,
                xr * ones(2),
                [minimum(gamd); maximum(gamd)];
                # linewidth=0.25,
                linestyle=:dash,
                color=mycolors[3],
                label="rotor location",
            )

            #save
            savefig(pb, "examples/body-velocity_xrotor$(xr)_tipgap$(tip_gap[i]).pdf")

            ##### ----- Plot Wake Strengths ----- #####
            pg = plot(;
                title="nbodypanel$(nducttot)_npersheet$(npersheet)_tipgap$(tip_gap[i])",
                xlabel=L"\gamma_\theta",
                ylabel="r",
            )

            # plot solution
            plot!(pg, gamw, inputs.rotor_panel_edges; label=convlabel)
            # plot initial
            plot!(
                pg,
                gamw_init,
                inputs.rotor_panel_edges;
                label="initial",
                linestyle=:dash,
                color=mycolors[1],
            )

            #save
            savefig(pg, "examples/wake-strength_xrotor$(xr)_tipgap$(tip_gap[i]).pdf")

            ##### ----- Plot Source Strengths ----- #####
            ps = plot(;
                title="nbodypanel$(nducttot)_npersheet$(npersheet)_tipgap$(tip_gap[i])",
                xlabel=L"\sigma",
                ylabel="r",
            )

            # plot solution
            plot!(ps, sigr, inputs.rotor_panel_centers; label=convlabel)
            # plot initial
            plot!(
                ps,
                sigr_init,
                inputs.rotor_panel_centers;
                label="initial",
                linestyle=:dash,
                color=mycolors[1],
            )

            #save
            savefig(
                ps, "examples/source-strength-check_xrotor$(xr)_tipgap$(tip_gap[i]).pdf"
            )

            ##### ----- Plot GEOMETRY ----- #####
            #initialize plot
            pgeom = plot(;
                aspectratio=1,
                title="nbodypanel$(nducttot)_npersheet$(npersheet)_tipgap$(tip_gap[i])",
            )
            plot!(
                pgeom,
                xr * ones(length(inputs.rotor_panel_edges)),
                inputs.rotor_panel_edges;
                color=mycolors[2],
                linewidth=0.25,
                markersize=0.5,
                markershape=:rect,
                label="",
            )

            plot!(
                pgeom,
                inputs.rotor_source_panels[1].panel_center[:, 1],
                inputs.rotor_source_panels[1].panel_center[:, 2];
                color=mycolors[3],
                seriestype=:scatter,
                markersize=0.75,
                markershape=:circle,
                label="",
            )

            for iw in 1:nwake_sheets
                plot!(
                    pgeom,
                    inputs.wakexgrid[:, iw],
                    inputs.wakergrid[:, iw];
                    linewidth=0.25,
                    markersize=0.5,
                    markershape=:rect,
                    color=:black,
                    label="",
                )

                plot!(
                    pgeom,
                    inputs.wake_vortex_panels[iw].panel_center[:, 1],
                    inputs.wake_vortex_panels[iw].panel_center[:, 2];
                    seriestype=:scatter,
                    markersize=0.75,
                    markershape=:circle,
                    color=mycolors[2],
                    label="",
                )
            end

            savefig(
                "examples/precomputed-rotor-wake-geometry_xrotor$(xr)_tipgap$(tip_gap[i]).pdf",
            )

            # plot body panels
            plot!(
                pgeom,
                inputs.duct_coordinates[:, 1],
                inputs.duct_coordinates[:, 2];
                linewidth=0.25,
                markersize=0.5,
                markershape=:rect,
                color=mycolors[3],
                label="",
            )

            # Plot body panel centers
            # for ib in 1:2
            for ib in 1:1
                plot!(
                    pgeom,
                    inputs.body_panels[ib].panel_center[:, 1],
                    inputs.body_panels[ib].panel_center[:, 2];
                    color=mycolors[1],
                    seriestype=:scatter,
                    markersize=0.75,
                    label="",
                )
            end

            savefig(
                pgeom,
                "examples/precomputed-full-geometry_xrotor$(xr)_tipgap$(tip_gap[i]).pdf",
            )
        end #for tip gap
    end #for xrotor

    println("DONE.")

    return nothing
end

full_clean(; xrotor=[0.25], tip_gap=[0.0], discscale=1, npanels_inlet=10, wake_length=1.0)

# full_clean(; tip_gap=[0.1; 0.01; 0.001; 0.0], discscale=1)

# ds = [1; 2; 4]
# # ds = [1]
# for d in 1:length(ds)
#     full_clean(; tip_gap=[0.01], discscale=ds[d])
# end
