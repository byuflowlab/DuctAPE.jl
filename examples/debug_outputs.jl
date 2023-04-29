#---------------------------------#
#             Includes            #
#---------------------------------#
project_dir = dirname(dirname(@__FILE__))

using DuctTAPE
const dt = DuctTAPE

using FLOWMath
const fm = FLOWMath

include("../plots_default.jl")

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
    vx_total, vr_total, vtheta_total, vx_frombody, vr_frombody, vx_fromwake, vr_fromwake, vx_fromrotor, vr_fromrotor = dt.calculate_induced_velocities_on_rotors(
        inputs.blade_elements,
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

    Wx_rotor = vx_total .+ inputs.Vinf

    Wtheta_rotor =
        vtheta_total .- inputs.blade_elements[1].Omega .* inputs.rotor_panel_centers

    Wm_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_total .^ 2)

    Wmag_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_total .^ 2 .+ Wtheta_rotor .^ 2)

    ##### ----- Plot vx distribution ----- #####
    # initialize plot
    pvx = plot(; xlabel=L"v_x", ylabel="r")

    # plot solution
    plot!(pvx, vx_total, inputs.rotor_panel_centers; label="Total")
    plot!(pvx, vx_frombody, inputs.rotor_panel_centers; label="Due to Body")
    plot!(pvx, vx_fromwake, inputs.rotor_panel_centers; label="Due to Wake")
    plot!(pvx, vx_fromrotor, inputs.rotor_panel_centers; label="Due to Rotor (sources)")

    #save
    savefig(pvx, "examples/vxdist" * suffix * ".pdf")

    ##### ----- Plot vtheta distribution ----- #####
    # initialize plot
    pvtheta = plot(; xlabel=L"v_\theta", ylabel="r")

    # plot solution
    plot!(pvtheta, vtheta_total, inputs.rotor_panel_centers; label="")

    #save
    savefig(pvtheta, "examples/vthetadist" * suffix * ".pdf")

    ##### ----- Plot vr distribution ----- #####
    # initialize plot
    pvr = plot(; xlabel=L"v_r", ylabel="r")

    # plot solution
    plot!(pvr, vr_total, inputs.rotor_panel_centers; label="Total")
    plot!(pvx, vr_frombody, inputs.rotor_panel_centers; label="Due to Body")
    plot!(pvx, vr_fromwake, inputs.rotor_panel_centers; label="Due to Wake")
    plot!(pvx, vr_fromrotor, inputs.rotor_panel_centers; label="Due to Rotor (sources)")

    #save
    savefig(pvr, "examples/vrdist" * suffix * ".pdf")

    ##### ----- Plot Wx distribution ----- #####
    # initialize plot
    pwx = plot(; xlabel=L"W_x", ylabel="r")

    # plot solution
    plot!(pwx, Wx_rotor, inputs.rotor_panel_centers; label="")

    #save
    savefig(pwx, "examples/Wxdist" * suffix * ".pdf")

    ##### ----- Plot Wtheta distribution ----- #####
    # initialize plot
    pwt = plot(; xlabel=L"W_\theta", ylabel="r")

    # plot solution
    plot!(pwt, Wtheta_rotor, inputs.rotor_panel_centers; label="")

    #save
    savefig(pwt, "examples/Wthetadist" * suffix * ".pdf")

    ##### ----- Plot Wm distribution ----- #####
    # initialize plot
    pwm = plot(; xlabel=L"W_m", ylabel="r")

    # plot solution
    plot!(pwm, Wm_rotor, inputs.rotor_panel_centers; label="")
    #
    #save
    savefig(pwm, "examples/Wmdist" * suffix * ".pdf")

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
        cl[ir] = fm.linear(
            [0.0; 1.0], [clin, clout], inputs.blade_elements[1].inner_fraction[ir]
        )
        cd[ir] = fm.linear(
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
    savefig(pa, "examples/angledist" * suffix * ".pdf")

    # plot chord distribution
    pc = plot(; xlabel="Chords", ylabel="r")
    plot!(pc, inputs.blade_elements.chords, inputs.rotor_panel_centers)
    savefig(pc, "examples/chorddist" * suffix * ".pdf")

    ### --- generate cl data plot --- ###
    # plot airfoil data!
    aoas = range(minimum(alpha), maximum(alpha), length(alpha) * 5)
    clrange, cdrange = dt.search_polars(airfoils[1], aoas)
    pafcl = plot(; xlabel="Angle of Attack", ylabel=L"c_\ell")
    plot!(pafcl, aoas * 180.0 / pi, clrange)
    savefig(pafcl, "examples/cldata" * suffix * ".pdf")
    pafcd = plot(; xlabel="Angle of Attack", ylabel=L"c_d")
    plot!(pafcd, aoas * 180.0 / pi, cdrange)
    savefig(pafcd, "examples/cddata" * suffix * ".pdf")
    #plot cl and cd

    ##### ----- Plot cl distribution ----- #####
    # initialize plot
    pcl = plot(; xlabel=L"c_\ell", ylabel="r")

    # plot solution
    plot!(pcl, cl, inputs.rotor_panel_centers; label="")

    #save
    savefig(pcl, "examples/cldist" * suffix * ".pdf")

    ##### ----- Plot cd distribution ----- #####
    # initialize plot
    pcd = plot(; xlabel=L"c_d", ylabel="r")

    # plot solution
    plot!(pcd, cd, inputs.rotor_panel_centers; label="")

    #save
    savefig(pcd, "examples/cddist" * suffix * ".pdf")

    return nothing
end

