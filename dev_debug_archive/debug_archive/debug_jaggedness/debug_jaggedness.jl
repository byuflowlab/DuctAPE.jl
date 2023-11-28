#---------------------------------#
#             Includes            #
#---------------------------------#
project_dir = dirname(dirname(dirname(@__FILE__)))

using DuctAPE
const dt = DuctAPE

# CCBlade used for it's airfoils function objects here.
using CCBlade
const ccb = CCBlade

using FLOWMath
const fm = FLOWMath

include(project_dir*"/plots_default.jl")

"""
"""
function initwake(inputs, blade_elements; vx_rotor=0.0, vtheta_rotor=0.0, vr_rotor=0.0, prefix="")

    # get dimensions
    nr = length(inputs.blade_elements[1].rbe)
    nrotor = length(inputs.blade_elements)

    # initialize states
    Gamr = zeros(nr, nrotor)
    sigr = zeros(nr, nrotor)
    gamw = zeros(nr + 1, nrotor)

    # Velocities
    Wx_rotor = vx_rotor .+ inputs.Vinf .* ones(nr)

    Wtheta_rotor =
        vtheta_rotor .- inputs.blade_elements[1].Omega .* inputs.rotor_panel_centers

    Wm_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2)

    Wmag_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2 .+ Wtheta_rotor .^ 2)

    #---------------------------------#
    #         Plot Velocities         #
    #---------------------------------#
    println("Plotting Velocities")

    ##### ----- Plot Wx distribution ----- #####
    # initialize plot
    pwx = plot(; xlabel=L"W_x", ylabel="r")

    # plot solution
    plot!(pwx, Wx_rotor, inputs.rotor_panel_centers; label="")

    #save
    savefig(pwx, "dev_debug_archive/debug_jaggedness/"*prefix*"Wxdist.pdf")

    ##### ----- Plot Wtheta distribution ----- #####
    # initialize plot
    pwt = plot(; xlabel=L"W_\theta", ylabel="r")

    # plot solution
    plot!(pwt, Wtheta_rotor, inputs.rotor_panel_centers; label="")

    #save
    savefig(pwt, "dev_debug_archive/debug_jaggedness/"*prefix*"Wthetadist.pdf")

    ##### ----- Plot Wm distribution ----- #####
    # initialize plot
    pwm = plot(; xlabel=L"W_m", ylabel="r")

    # plot solution
    plot!(pwm, Wm_rotor, inputs.rotor_panel_centers; label="")

    #save
    savefig(pwm, "dev_debug_archive/debug_jaggedness/"*prefix*"Wmdist.pdf")

    # initialize circulation and source panel strengths
    Gamr, sigr, phi, alpha, cl, cd = dt.calculate_gamma_sigma!(
        Gamr, sigr, inputs.blade_elements, Wm_rotor, Wtheta_rotor, Wmag_rotor; debug=true
    )

    #---------------------------------#
    #  Plot Circulation Constituents  #
    #---------------------------------#
    println("Plotting Circulation Constituents")

    ##### ----- Plot angle of attack distribution ----- #####
    # initialize plot
    pa = plot(; xlabel="Angles", ylabel="r")

    # plot twist
    plot!(
        pa,
        blade_elements.twists * 180.0 / pi,
        inputs.rotor_panel_centers;
        label="Blade Twist",
    )
    # plot inflow
    plot!(pa, phi * 180.0 / pi, inputs.rotor_panel_centers; label="Inflow Angle")
    # plot aoa
    plot!(pa, alpha * 180.0 / pi, inputs.rotor_panel_centers; label="Angle of Attack")

    #save
    savefig(pa, "dev_debug_archive/debug_jaggedness/"*prefix*"angledist.pdf")

    # plot chord distribution
    pc = plot(; xlabel="Chords", ylabel="r")
    plot!(pc, blade_elements.chords, inputs.rotor_panel_centers)
    savefig(pc, "dev_debug_archive/debug_jaggedness/"*prefix*"chorddist.pdf")

    ### --- generate cl data plot --- ###
    # plot airfoil data!
    aoas = range(minimum(alpha), maximum(alpha), length(alpha) * 5)
    clrange, cdrange = dt.search_polars(airfoils[1], aoas)
    pafcl = plot(; xlabel="Angle of Attack", ylabel=L"c_\ell")
    plot!(pafcl, aoas * 180.0 / pi, clrange)
    savefig(pafcl, "dev_debug_archive/debug_jaggedness/"*prefix*"cldata.pdf")
    pafcd = plot(; xlabel="Angle of Attack", ylabel=L"c_d")
    plot!(pafcd, aoas * 180.0 / pi, cdrange)
    savefig(pafcd, "dev_debug_archive/debug_jaggedness/"*prefix*"cddata.pdf")

    ##### ----- Plot cl distribution ----- #####
    # initialize plot
    pcl = plot(; xlabel=L"c_\ell", ylabel="r")

    # plot solution
    plot!(pcl, cl, inputs.rotor_panel_centers; label="")

    #save
    savefig(pcl, "dev_debug_archive/debug_jaggedness/"*prefix*"cldist.pdf")

    ##### ----- Plot cd distribution ----- #####
    # initialize plot
    pcd = plot(; xlabel=L"c_d", ylabel="r")

    # plot solution
    plot!(pcd, cd, inputs.rotor_panel_centers; label="")

    #save
    savefig(pcd, "dev_debug_archive/debug_jaggedness/"*prefix*"cddist.pdf")

    # - Calculate net circulation and enthalpy jumps - #
    Gamma_tilde = dt.calculate_net_circulation(Gamr, inputs.blade_elements.B)
    H_tilde = dt.calculate_enthalpy_jumps(
        Gamr, inputs.blade_elements.Omega, inputs.blade_elements.B
    )

    #---------------------------------#
    #      Plot Wake Constituents     #
    #---------------------------------#

    println("Plotting Wake Strength Constituents")

    ##### ----- Plot Enthalpy Jumps ----- #####
    ph = plot(; xlabel=L"\widetilde{H}", ylabel="r")

    plot!(ph, H_tilde, inputs.rotor_panel_centers; label="")

    #save
    savefig(ph, "dev_debug_archive/debug_jaggedness/"*prefix*"htildedist.pdf")

    ##### ----- Plot Net Circulation ----- #####
    pGt = plot(; xlabel=L"\widetilde{\Gamma}", ylabel="r")

    plot!(pGt, Gamma_tilde, inputs.rotor_panel_centers; label="")

    #save
    savefig(pGt, "dev_debug_archive/debug_jaggedness/"*prefix*"Gammatildedist.pdf")

    # - update wake strengths - #
    gamw, dhtilde, dgammatilde2, kdb = dt.calculate_wake_vortex_strengths!(
        gamw, inputs.rotor_panel_edges, Wm_rotor, Gamma_tilde, H_tilde; debug=true
    )

    #---------------------------------#
    #       Plot more wake stuff      #
    #---------------------------------#
    println("Plotting More Wake Stuff")

    pdh = plot(; xlabel=L"\Delta \widetilde{H}", ylabel="r")
    plot!(pdh, dhtilde, inputs.rotor_panel_edges; label="")
    savefig(pdh, "dev_debug_archive/debug_jaggedness/"*prefix*"DeltaHtildedist.pdf")

    pdg2 = plot(; xlabel=L"\Delta \widetilde{\Gamma}^2", ylabel="r")
    plot!(pdg2, dgammatilde2, inputs.rotor_panel_edges; label="")
    savefig(pdg2, "dev_debug_archive/debug_jaggedness/"*prefix*"DeltaGamtilde2dist.pdf")

    pk = plot(; xlabel="constant dgamma is multiplied by", ylabel="r")
    plot!(pk, kdb, inputs.rotor_panel_edges; label="")
    savefig(pk, "dev_debug_archive/debug_jaggedness/"*prefix*"K.pdf")

    pKGt = plot(; xlabel=L"K\Delta \widetilde{\Gamma}^2", ylabel="r")
    plot!(pKGt, kdb .* dgammatilde2, inputs.rotor_panel_edges; label="")
    savefig(pKGt, "dev_debug_archive/debug_jaggedness/"*prefix*"KdG2.pdf")

    #---------------------------------#
    #          Plot States            #
    #---------------------------------#
    println("Plotting Initial States")

    ##### ----- Plot rotor circulation distribution ----- #####
    # initialize plot
    pG = plot(; xlabel=L"\Gamma", ylabel="r")

    # plot solution
    plot!(pG, Gamr, inputs.rotor_panel_centers; label="")

    #save
    savefig(pG, "dev_debug_archive/debug_jaggedness/"*prefix*"circulationdist.pdf")

    ##### ----- Plot Wake Strengths ----- #####
    pg = plot(; xlabel=L"\gamma_\theta", ylabel="r")

    # plot solution
    plot!(pg, gamw, inputs.rotor_panel_edges; label="")

    #save
    savefig(pg, "dev_debug_archive/debug_jaggedness/"*prefix*"gammadist.pdf")

    ##### ----- Plot Source Strengths ----- #####
    ps = plot(; xlabel=L"\sigma", ylabel="r")

    plot!(ps, sigr, inputs.rotor_panel_centers; label="")

    #save
    savefig(ps, "dev_debug_archive/debug_jaggedness/"*prefix*"sigmadist.pdf")

    return nothing
end

#---------------------------------#
#            Constants            #
#---------------------------------#

# Blade Tip Radius, in meters
Rtip = 10 / 2.0 * 0.0254  # inches to meters

# Blade Hub radius, in meters
Rhub = 0.10 * Rtip

# number of blades
B = 2

# Blade section non-dimensional radial positions, chords lengths, and local twists angles in degrees
propgeom = [
    # 0.15 0.130 32.76
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
#
twists = propgeom[:, 3] * pi / 180

plot(rnondim, chords ./ Rtip; xlabel="r/R", ylabel="c/R")
savefig("dev_debug_archive/debug_jaggedness/"*prefix*"chord_dist_raw.pdf")
plot(rnondim, twists * 180 / pi; xlabel="r/R", ylabel="Twist (deg)")
savefig("dev_debug_archive/debug_jaggedness/"*prefix*"twist_dist_raw.pdf")

# use a NACA 4412 airfoils
airfoil_file = project_dir * "/test/data/xrotor_af_test.dat"
airfoils = fill(ccb.AlphaAF(airfoil_file), length(rnondim))

#Vinf
Vinf = 5.0

# rotor rotation rate in rad/s
Omega = 5400 * pi / 30  # convert from RPM to rad/s

nwake_sheets = 10

#---------------------------------#
#          Define Inputs          #
#---------------------------------#
# Rotor Parameters
rotor_parameters = [(;
    xrotor=0.0,
    nwake_sheets,
    r=rnondim,
    chords,
    twists,
    airfoils,
    Rtip,
    Rhub,
    tip_gap=0.0,
    B,
    Omega,
)]

# Freestream Parameters
freestream = (; Vinf)

_, inputs = dt.initialize_rotor_states(rotor_parameters, freestream; wake_length=1.0)

# look at initial case without body
initwake(inputs, inputs.blade_elements; vx_rotor=0.0, vtheta_rotor=0.0, vr_rotor=0.0, prefix="")

# look at converged case without body

# look at initial case with body nearby

