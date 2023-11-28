#---------------------------------#
#             Includes            #
#---------------------------------#

project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

include(project_dir * "/plots_default.jl")

using DuctAPE
const dt = DuctAPE
include(project_dir * "/dev_debug_archive/duct_only/run_duct.jl")

# CCBlade used for it's airfoils function objects here.
using CCBlade
const ccb = CCBlade
include(project_dir * "/dev_debug_archive/rotor_only/run_ccblade.jl")

using FLOWMath
const fm = FLOWMath

include(project_dir * "/dev_debug_archive/debug_together/debug_outputs.jl")

function step_through(;
    iterations=1,
    tip_gap=0.0,
    xrotor=0.25,
    npanels_inlet=10,
    wake_length=1.0,
    discscale=1,
    duct_file=project_dir * "/test/data/naca_662-015.jl",
    airfoil_file=project_dir * "/test/data/xrotor_af_test.dat",
    debug=true,
)

    #---------------------------------#
    #         ROTOR Geometry          #
    #---------------------------------#

    # Blade Tip Radius, in meters
    Rtip = 10 / 2.0 * 0.0254  # inches to meters
    tg = tip_gap .* Rtip / 100.0

    # Blade Hub radius, in meters
    Rhub = 0.10 * Rtip

    # number of blades
    B = 2

    # x position of rotor
    # xrotor = 0.25

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

    nwake_sheets = 10

    # non-dimensional wake length
    # wake_length = 1.0

    #= paneling for which zero tip gap case converges
    npanels = [10; 25]
    nhub_inlet = 20
    nduct_inlet = 20
    wakeoutratio = 25 / 10
    outletinletratio = 0.5
    =#

    ductle = minimum(duct_coordinates[:, 1])
    ductte = maximum(duct_coordinates[:, 1])
    ductchord = maximum(duct_coordinates[:, 1]) - minimum(duct_coordinates[:, 1])
    outletinletratio = (ductte - xrotor) / (xrotor - ductle)

    nhub_inlet = round(Int, npanels_inlet * discscale)

    nduct_inlet = round(Int, npanels_inlet * discscale)

    nduct_outlet = round(Int, nduct_inlet * outletinletratio)

    nwake = round(Int, (nduct_inlet + nduct_outlet) * wake_length)

    npanels = [nduct_outlet, nwake]

    nducttot = (nduct_inlet + nduct_outlet) * 2
    npersheet = nduct_outlet + nwake

    #---------------------------------#
    #          Define Inputs          #
    #---------------------------------#

    # Rotor Parameters
    rotor_parameters = [(;
        xrotor=xrotor,
        nwake_sheets,
        r=rnondim,
        chords,
        twists,
        airfoils,
        Rtip,
        Rhub,
        tip_gap,
        B,
        Omega,
    )]

    # Paneling Parameters
    paneling_constants = (; npanels, nhub_inlet, nduct_inlet, wake_length, nwake_sheets)

    # Freestream Parameters
    freestream = (; rho, mu, asound, Vinf)

    ##### ----- Setup and make copy of states that won't get changed ----- #####
    # initialize various inputs used in analysis
    inputs = dt.precomputed_inputs(
        duct_coordinates,
        hub_coordinates,
        paneling_constants,
        rotor_parameters,
        freestream;
        debug=debug,
    )

    println("ductTE_index = ", inputs.ductTE_index)
    println("hubTE_index = ", inputs.hubTE_index)

    initial_states = dt.initialize_states(inputs)
    initials = copy(initial_states)
    # - Extract initial states - #
    gamb_init, gamw_init, Gamr_init, sigr_init = dt.extract_state_variables(
        initials, inputs
    )

    #---------------------------------------------
    #run isolated cases
    ccbouts = run_ccblade(Vinf)
    ductouts = run_duct(inputs.duct_coordinates)
    #---------------------------------------------

    ## -- Plot Initial States (and set up plot objects) -- ##
    # - Body Strengths - #
    #split inner/outer
    dp = inputs.body_panels[1].panel_center[:, 1]
    _, leidx = findmin(dp)
    dpinner = dp[1:leidx]
    dpouter = dp[(leidx + 1):end]
    gamb_init_inner = gamb_init[1:leidx]
    gamb_init_outer = gamb_init[(leidx + 1):end]

    #plot initial body strengths
    pgb = plot(; xlabel="x", ylabel=L"\gamma_\theta^{body}")

    plot!(
        pgb,
        dpinner,
        gamb_init_inner;
        color=mycolors[1],
        linestyle=:dot,
        label="initial, inner surface",
    )
    plot!(
        pgb,
        dpouter,
        gamb_init_outer;
        color=mycolors[2],
        linestyle=:dot,
        label="initial, outer surface",
    )

    #plot rotor location for reference
    plot!(
        pgb,
        xrotor * ones(2),
        [minimum(gamb_init); maximum(gamb_init)];
        linestyle=:dot,
        color=mycolors[3],
        label="rotor location",
    )

    # - Wake Strengths - #
    pgw = plot(; xlabel=L"\gamma_\theta^{wake}", ylabel="r")

    # plot initial
    plot!(
        pgw,
        gamw_init,
        inputs.rotor_panel_edges;
        label="initial",
        linestyle=:dash,
        color=mycolors[1],
    )

    # - Circulation Strengths - #
    # initialize plot
    pG = plot(; xlabel=L"\Gamma", ylabel="r")

    # plot initial
    plot!(
        pG,
        Gamr_init,
        inputs.rotor_panel_centers;
        label="initial",
        linestyle=:dash,
        color=mycolors[1],
    )

    # - Source Strengths - #
    ps = plot(; xlabel=L"\sigma", ylabel="r")

    # plot initial
    plot!(
        ps,
        sigr_init,
        inputs.rotor_panel_centers;
        label="initial",
        linestyle=:dash,
        color=mycolors[1],
    )

    ## -- Initialize other plots -- ##
    # - Induced velocity plots -- #
    pvx = plot(; xlabel="Induced Axial Velocity", ylabel="r")
    pvr = plot(; xlabel="Induced Radial Velocity", ylabel="r")
    pvt = plot(; xlabel="Induced Tangential Velocity", ylabel="r")
    # - Total velocity plots -- #
    pwm = plot(; xlabel="Total Meridional Velocity", ylabel="r")
    plot!(
        pwm,
        inputs.Vinf .* ones(length(inputs.rotor_panel_centers)),
        inputs.rotor_panel_centers;
        label="freestream for reference",
        linestyle=:dot,
        color=:black,
    )
    pwmag = plot(; xlabel="Velocity Magnitude", ylabel="r")
    pwt = plot(; xlabel="Total Tangential Velocity", ylabel="r")
    plot!(
        pwt,
        -inputs.blade_elements[1].Omega .* inputs.rotor_panel_centers,
        inputs.rotor_panel_centers;
        label="rotational for reference",
        linestyle=:dot,
        color=:black,
    )
    # - Wake constituent plots -- #
    pGT = plot(; xlabel=L"\widetilde{\Gamma}", ylabel="r")
    pHT = plot(; xlabel=L"\widetilde{H}", ylabel="r")

    for i in 1:iterations
        println("\n\nITERATION $i\n")

        ##### ----- Enter Residual Function ----- #####
        updated_states = copy(initial_states)

        ### --- Enter update function --- ###
        blade_elements = inputs.blade_elements
        rpc = inputs.rotor_panel_centers
        Vinf = inputs.Vinf

        # - Extract states - #
        gamb, gamw, Gamr, sigr = dt.extract_state_variables(updated_states, inputs)

        # - Fill out wake strengths - #
        wake_vortex_strengths = dt.fill_out_wake_strengths(
            gamw,
            inputs.rotor_indices,
            inputs.num_wake_x_panels;
            ductTE_index=10,
            interface="hard",
        )

        println("\tCalculating and Plotting Velocities")
        # - Get the induced velocities at the rotor plane - #
        vx_rotor, vr_rotor, vtheta_rotor, vxb, vrb, vxw, vrw, vxr, vrr = dt.calculate_induced_velocities_on_rotors(
            blade_elements,
            Gamr,
            inputs.vx_rw,
            inputs.vr_rw,
            wake_vortex_strengths,
            inputs.vx_rr,
            inputs.vr_rr,
            sigr,
            inputs.vx_rb,
            inputs.vr_rb,
            gamb;
            debug=true,
        )

        # the axial component also includes the freestream velocity ( see eqn 1.87 in dissertation)
        Wx_rotor = vx_rotor .+ inputs.Vinf

        # the tangential also includes the negative of the rotation rate (see eqn 1.87 in dissertation)
        Wtheta_rotor = vtheta_rotor .- inputs.blade_elements[1].Omega .* rpc

        # meridional component
        Wm_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2)

        # Get the inflow magnitude at the rotor as the combination of all the components
        Wmag_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2 .+ Wtheta_rotor .^ 2)

        #---------------------------------#
        #               PLOT              #
        #---------------------------------#

        #plot axial induced velocity
        plot!(
            pvx,
            vx_rotor,
            inputs.rotor_panel_centers;
            color=mycolors[i],
            label="total iter $i",
        )
        plot!(
            pvx,
            vxb,
            inputs.rotor_panel_centers;
            color=mycolors[i],
            linestyle=:dash,
            label="from body",
        )
        plot!(
            pvx,
            vxw,
            inputs.rotor_panel_centers;
            color=mycolors[i],
            linestyle=:dot,
            label="from wake",
        )

        #plot radial induced velocity
        plot!(
            pvr,
            vr_rotor,
            inputs.rotor_panel_centers;
            color=mycolors[i],
            label="total iter $i",
        )
        plot!(
            pvr,
            vrb,
            inputs.rotor_panel_centers;
            color=mycolors[i],
            linestyle=:dash,
            label="from body",
        )
        plot!(
            pvr,
            vrw,
            inputs.rotor_panel_centers;
            color=mycolors[i],
            linestyle=:dot,
            label="from wake",
        )

        #plot tangential induced velocity
        plot!(
            pvt,
            vtheta_rotor,
            inputs.rotor_panel_centers;
            color=mycolors[i],
            label="total iter $i",
        )

        #plot meridional velocity
        plot!(
            pwm,
            Wm_rotor,
            inputs.rotor_panel_centers;
            color=mycolors[i],
            label="total iter $i",
        )

        #plot tangential velocity
        plot!(
            pwt,
            Wtheta_rotor,
            inputs.rotor_panel_centers;
            color=mycolors[i],
            label="total iter $i",
        )

        #plot magnitude of velocity
        plot!(
            pwmag,
            Wmag_rotor,
            inputs.rotor_panel_centers;
            color=mycolors[i],
            label="total iter $i",
        )

        #---------------------------------#
        #            CONTINUE             #
        #---------------------------------#

        println("\tCalculating and Plotting Body Strengths")
        # - Calculate body vortex strengths (before updating state dependencies) - #
        dt.calculate_body_vortex_strengths!(
            gamb,
            inputs.A_bb,
            inputs.b_bf,
            inputs.kutta_idxs,
            inputs.A_bw,
            wake_vortex_strengths,
            inputs.A_br,
            sigr,
        )

        #---------------------------------#
        #               PLOT              #
        #---------------------------------#

        plot!(pgb, dpinner, gamb[1:leidx]; color=mycolors[i], label="inner surface iter $i")
        plot!(
            pgb,
            dpouter,
            gamb[(leidx + 1):end];
            color=mycolors[i],
            linestyle=:dash,
            label="outer surface iter $i",
        )

        #---------------------------------#
        #            CONTINUE             #
        #---------------------------------#

        println("\tCalculating and Plotting Wake Strengths")

        # - Calculate net circulation and enthalpy jumps - #
        Gamma_tilde = dt.calculate_net_circulation(Gamr, blade_elements.B)
        H_tilde = dt.calculate_enthalpy_jumps(Gamr, blade_elements.Omega, blade_elements.B)

        # - update wake strengths - #
        dt.calculate_wake_vortex_strengths!(
            gamw, inputs.rotor_panel_edges, Wm_rotor, Gamma_tilde, H_tilde
        )

        #---------------------------------#
        #               PLOT              #
        #---------------------------------#

        #gamma tilde
        plot!(
            pGT, Gamma_tilde, inputs.rotor_panel_centers; label="iter $i", color=mycolors[i]
        )

        #h tildet
        plot!(pHT, H_tilde, inputs.rotor_panel_centers; label="iter $i", color=mycolors[i])

        # wake strengths
        plot!(pgw, gamw, inputs.rotor_panel_edges; label="iter $i", color=mycolors[i])

        #---------------------------------#
        #            CONTINUE             #
        #---------------------------------#

        println("\tCalculating and Plotting Circulation and Source Strengths")
        # - Update rotor circulation and source panel strengths - #
        dt.calculate_gamma_sigma!(
            Gamr, sigr, blade_elements, Wm_rotor, Wtheta_rotor, Wmag_rotor
        )

        #---------------------------------#
        #               PLOT              #
        #---------------------------------#

        plot!(pG, Gamr, inputs.rotor_panel_centers; label="iter $i", color=mycolors[i])
        plot!(ps, sigr, inputs.rotor_panel_centers; label="iter $i", color=mycolors[i])

        ### --- END OF UPDATE FUNCTION --- ###

        res = updated_states - initial_states
        println("\tIninity Norm of Residual = ", maximum(res))

        ##### ----- END OF RESIDUAL FUNCTION ----- #####
        # - update states for next iteration
        initial_states .= updated_states
    end
    ## -- Save Figures -- ##
    println("\nSaving Figures...")
    # - states - #
    savefig(pgb, project_dir * "/dev_debug_archive/step_through/body-strengths.pdf")
    savefig(pgw, project_dir * "/dev_debug_archive/step_through/wake-strengths.pdf")
    savefig(ps, project_dir * "/dev_debug_archive/step_through/source-strengths.pdf")
    savefig(pG, project_dir * "/dev_debug_archive/step_through/circulation-strengths.pdf")

    # - velocities - #
    savefig(pvx, project_dir * "/dev_debug_archive/step_through/vx.pdf")
    savefig(pvr, project_dir * "/dev_debug_archive/step_through/vr.pdf")
    savefig(pvt, project_dir * "/dev_debug_archive/step_through/vtheta.pdf")
    savefig(pwm, project_dir * "/dev_debug_archive/step_through/wm.pdf")
    savefig(pwmag, project_dir * "/dev_debug_archive/step_through/wmag.pdf")
    savefig(pwt, project_dir * "/dev_debug_archive/step_through/wtheta.pdf")

    return nothing
end

step_through(;
    xrotor=0.25, tip_gap=0.0, discscale=1, npanels_inlet=10, wake_length=1.0, iterations=7
)
