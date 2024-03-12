#---------------------------------#
#             Includes            #
#---------------------------------#
# project_dir = dirname(dirname(dirname(@__FILE__)))
# if project_dir == ""
#     project_dir = "."
# end

# include(project_dir * "/plots_default.jl")

# using DuctAPE
# const dt = DuctAPE

# using FLOWMath
# const fm = FLOWMath

function debug_outputs(states, inputs; suffix="")

    #extract outputs
    gamb, gamw, Gamr, sigr = dt.extract_state_variables(states, inputs)

    #---------------------------------#
    #         Plot Velocities         #
    #---------------------------------#
    println("Plotting Velocities")

    # - Fill out wake strengths - #
    wake_vortex_strengths = dt.fill_out_wake_strengths(
        gamw, inputs.rotor_indices, inputs.num_wake_x_panels
    )

    # get induced velocities
    vz_total, vr_total, vtheta_total, vz_frombody, vr_frombody, vz_fromwake, vr_fromwake, vz_fromrotor, vr_fromrotor = dt.calculate_induced_velocities_on_rotors(
        inputs.blade_elements,
        Gamr,
        inputs.vz_rw,
        inputs.vr_rw,
        wake_vortex_strengths,
        inputs.vz_rr,
        inputs.vr_rr,
        sigr,
        inputs.vz_rb,
        inputs.vr_rb,
        gamb;
        debug=true,
    )

    Wx_rotor = vz_total .+ inputs.Vinf

    Wtheta_rotor =
        vtheta_total .- inputs.blade_elements[1].Omega .* inputs.rotor_panel_centers

    Wm_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_total .^ 2)

    Wmag_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_total .^ 2 .+ Wtheta_rotor .^ 2)

    ##### ----- Plot vx distribution ----- #####
    # initialize plot
    pvx = plot(; xlabel=L"v_x", ylabel="r")

    # plot solution
    plot!(pvx, vz_total, inputs.rotor_panel_centers; label="Total")
    plot!(pvx, vz_frombody, inputs.rotor_panel_centers; label="Due to Body")
    plot!(pvx, vz_fromwake, inputs.rotor_panel_centers; label="Due to Wake")
    plot!(pvx, vz_fromrotor, inputs.rotor_panel_centers; label="Due to Rotor (sources)")

    #save
    savefig(pvx, "dev_debug_archive/debug_together/vxdist" * suffix * ".pdf")

    ##### ----- Plot vtheta distribution ----- #####
    # initialize plot
    pvtheta = plot(; xlabel=L"v_\theta", ylabel="r")

    # plot solution
    plot!(pvtheta, vtheta_total, inputs.rotor_panel_centers; label="")

    #save
    savefig(pvtheta, "dev_debug_archive/debug_together/vthetadist" * suffix * ".pdf")

    ##### ----- Plot vr distribution ----- #####
    # initialize plot
    pvr = plot(; xlabel=L"v_r", ylabel="r")

    # plot solution
    plot!(pvr, vr_total, inputs.rotor_panel_centers; label="Total")
    plot!(pvx, vr_frombody, inputs.rotor_panel_centers; label="Due to Body")
    plot!(pvx, vr_fromwake, inputs.rotor_panel_centers; label="Due to Wake")
    plot!(pvx, vr_fromrotor, inputs.rotor_panel_centers; label="Due to Rotor (sources)")

    #save
    savefig(pvr, "dev_debug_archive/debug_together/vrdist" * suffix * ".pdf")

    ##### ----- Plot Wx distribution ----- #####
    # initialize plot
    pwx = plot(; xlabel=L"W_x", ylabel="r")

    # plot solution
    plot!(pwx, Wx_rotor, inputs.rotor_panel_centers; label="")

    #save
    savefig(pwx, "dev_debug_archive/debug_together/Wxdist" * suffix * ".pdf")

    ##### ----- Plot Wtheta distribution ----- #####
    # initialize plot
    pwt = plot(; xlabel=L"W_\theta", ylabel="r")

    # plot solution
    plot!(pwt, Wtheta_rotor, inputs.rotor_panel_centers; label="")

    #save
    savefig(pwt, "dev_debug_archive/debug_together/Wthetadist" * suffix * ".pdf")

    ##### ----- Plot Wm distribution ----- #####
    # initialize plot
    pwm = plot(; xlabel=L"W_m", ylabel="r")

    # plot solution
    plot!(pwm, Wm_rotor, inputs.rotor_panel_centers; label="")
    #
    #save
    savefig(pwm, "dev_debug_archive/debug_together/Wmdist" * suffix * ".pdf")

    #---------------------------------#
    #  Plot Circulation Constituents  #
    #---------------------------------#
    println("Plotting Circulation Constituents")

    phi = atan.(Wm_rotor, -Wtheta_rotor)
    alpha = inputs.blade_elements[1].twists .- phi

    cl = zeros(length(alpha))
    cd = zeros(length(alpha))

    for ir in 1:length(alpha)
        # look up lift and drag data for the nearest two input sections
        clin, cdin = dt.search_polars(inputs.blade_elements[1].inner_airfoil[ir], alpha[ir])
        clout, cdout = dt.search_polars(
            inputs.blade_elements[1].outer_airfoil[ir], alpha[ir]
        )
        # linearly interpolate between those two values at your blade element location
        cl[ir] = FLOWMath.linear(
            [0.0; 1.0], [clin, clout], inputs.blade_elements[1].inner_fraction[ir]
        )
        cd[ir] = FLOWMath.linear(
            [0.0; 1.0], [cdin, cdout], inputs.blade_elements[1].inner_fraction[ir]
        )
    end

    ##### ----- Plot angle of attack distribution ----- #####
    # initialize plot
    pa = plot(; xlabel="Angles", ylabel="r")

    # plot twist
    plot!(
        pa,
        inputs.blade_elements.twists * 180.0 / pi,
        inputs.rotor_panel_centers;
        label="Blade Twist",
    )
    # plot inflow
    plot!(pa, phi * 180.0 / pi, inputs.rotor_panel_centers; label="Inflow Angle")
    # plot aoa
    plot!(pa, alpha * 180.0 / pi, inputs.rotor_panel_centers; label="Angle of Attack")

    #save
    savefig(pa, "dev_debug_archive/debug_together/angledist" * suffix * ".pdf")

    # plot chord distribution
    pc = plot(; xlabel="Chords", ylabel="r")
    plot!(pc, inputs.blade_elements.chords, inputs.rotor_panel_centers)
    savefig(pc, "dev_debug_archive/debug_together/chorddist" * suffix * ".pdf")

    ### --- generate cl data plot --- ###
    # plot airfoil data!
    aoas = range(minimum(alpha), maximum(alpha), length(alpha) * 5)
    clrange, cdrange = dt.search_polars(inputs.blade_elements[1].outer_airfoil[1], aoas)
    pafcl = plot(; xlabel="Angle of Attack", ylabel=L"c_\ell")
    plot!(pafcl, aoas * 180.0 / pi, clrange)
    savefig(pafcl, "dev_debug_archive/debug_together/cldata" * suffix * ".pdf")
    pafcd = plot(; xlabel="Angle of Attack", ylabel=L"c_d")
    plot!(pafcd, aoas * 180.0 / pi, cdrange)
    savefig(pafcd, "dev_debug_archive/debug_together/cddata" * suffix * ".pdf")
    #plot cl and cd

    ##### ----- Plot cl distribution ----- #####
    # initialize plot
    pcl = plot(; xlabel=L"c_\ell", ylabel="r")

    # plot solution
    plot!(pcl, cl, inputs.rotor_panel_centers; label="")

    #save
    savefig(pcl, "dev_debug_archive/debug_together/cldist" * suffix * ".pdf")

    ##### ----- Plot cd distribution ----- #####
    # initialize plot
    pcd = plot(; xlabel=L"c_d", ylabel="r")

    # plot solution
    plot!(pcd, cd, inputs.rotor_panel_centers; label="")

    #save
    savefig(pcd, "dev_debug_archive/debug_together/cddist" * suffix * ".pdf")

    #---------------------------------#
    #      Plot Influences on Body    #
    #---------------------------------#
    println("Plotting Influences on Body (RHS)")

    bfree, bwake, brotor = dt.calculate_body_vortex_strengths!(
        gamb,
        inputs.A_bb,
        inputs.b_bf,
        inputs.kutta_idxs,
        inputs.A_bw,
        wake_vortex_strengths,
        inputs.A_br,
        sigr;
        debug=true,
    )

    ##### ----- Plot duct surface velocity ----- #####

    #prepare outputs
    dp = inputs.body_panels[1].panel_center[:, 1]
    _, leidx = findmin(dp)
    gamd = 1.0 .- (gamb[1:length(dp)] ./ inputs.Vinf) .^ 2
    gamd = gamb[1:length(dp)] ./ inputs.Vinf
    #split into inner and outer surfaces
    dpinner = dp[1:leidx]
    dpouter = dp[(leidx + 1):end]
    bfreeinner = bfree[1:leidx]
    bfreeouter = bfree[(leidx + 1):end]
    bwakeinner = bwake[1:leidx]
    bwakeouter = bwake[(leidx + 1):end]
    brotorinner = brotor[1:leidx]
    brotorouter = brotor[(leidx + 1):end]

    # initialize plot
    pb = plot(; xlabel="x", ylabel="Linear System RHS components")

    # plot solution
    plot!(pb, dpinner, bfreeinner; label="inner surface, freestream")
    plot!(pb, dpouter, bfreeouter; label="outer surface, freestream")
    plot!(pb, dpinner, bwakeinner; label="inner surface, wake")
    plot!(pb, dpouter, bwakeouter; label="outer surface, wake")
    plot!(pb, dpinner, brotorinner; label="inner surface, rotor")
    plot!(pb, dpouter, brotorouter; label="outer surface, rotor")

    savefig(pb, "dev_debug_archive/debug_together/linearsystemRHS" * suffix * ".pdf")

    return nothing
end

