#---------------------------------#
#             Includes            #
#---------------------------------#
project_dir = dirname(dirname(@__FILE__))

using DuctTAPE
const dt = DuctTAPE

# CCBlade used for it's airfoils function objects here.
using CCBlade
const ccb = CCBlade

using FLOWMath
const fm = FLOWMath

include("../plots_default.jl")

"""
"""
function initwake(inputs, blade_elements; vx_rotor=0.0, vtheta_rotor=0.0, vr_rotor=0.0)

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
    savefig(pwx, "examples/Wxdist.pdf")

    ##### ----- Plot Wtheta distribution ----- #####
    # initialize plot
    pwt = plot(; xlabel=L"W_\theta", ylabel="r")

    # plot solution
    plot!(pwt, Wtheta_rotor, inputs.rotor_panel_centers; label="")

    #save
    savefig(pwt, "examples/Wthetadist.pdf")

    ##### ----- Plot Wm distribution ----- #####
    # initialize plot
    pwm = plot(; xlabel=L"W_m", ylabel="r")

    # plot solution
    plot!(pwm, Wm_rotor, inputs.rotor_panel_centers; label="")

    #save
    savefig(pwm, "examples/Wmdist.pdf")

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
    savefig(pa, "examples/angledist.pdf")

    ### --- generate cl data plot --- ###
    # plot airfoil data!
    aoas = range(minimum(alpha), maximum(alpha), length(alpha) * 5)
    clrange, cdrange = dt.search_polars(airfoils[1], aoas)
    pafcl = plot(; xlabel="Angle of Attack", ylabel=L"c_\ell")
    plot!(pafcl, aoas * 180.0 / pi, clrange)
    savefig(pafcl, "examples/cldata.pdf")
    pafcd = plot(; xlabel="Angle of Attack", ylabel=L"c_d")
    plot!(pafcd, aoas * 180.0 / pi, cdrange)
    savefig(pafcd, "examples/cddata.pdf")

    ##### ----- Plot cl distribution ----- #####
    # initialize plot
    pcl = plot(; xlabel=L"c_\ell", ylabel="r")

    # plot solution
    plot!(pcl, cl, inputs.rotor_panel_centers; label="")

    #save
    savefig(pcl, "examples/cldist.pdf")

    ##### ----- Plot cd distribution ----- #####
    # initialize plot
    pcd = plot(; xlabel=L"c_d", ylabel="r")

    # plot solution
    plot!(pcd, cd, inputs.rotor_panel_centers; label="")

    #save
    savefig(pcd, "examples/cddist.pdf")

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
    savefig(ph, "examples/htildedist.pdf")

    ##### ----- Plot Net Circulation ----- #####
    pGt = plot(; xlabel=L"\widetilde{\Gamma}", ylabel="r")

    plot!(pGt, Gamma_tilde, inputs.rotor_panel_centers; label="")

    #save
    savefig(pGt, "examples/Gammatildedist.pdf")

    # - update wake strengths - #
    dt.calculate_wake_vortex_strengths!(
        gamw, inputs.rotor_panel_edges, Wm_rotor, Gamma_tilde, H_tilde
    )

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
    savefig(pG, "examples/circulationdist.pdf")

    ##### ----- Plot Wake Strengths ----- #####
    pg = plot(; xlabel=L"\gamma_\theta", ylabel="r")

    # plot solution
    plot!(pg, gamw, inputs.rotor_panel_edges; label="")

    #save
    savefig(pg, "examples/gammadist.pdf")

    ##### ----- Plot Source Strengths ----- #####
    ps = plot(; xlabel=L"\sigma", ylabel="r")

    plot!(ps, sigr, inputs.rotor_panel_centers; label="")

    #save
    savefig(ps, "examples/sigmadist.pdf")

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
#
twists = propgeom[:, 3] * pi / 180

# use a NACA 4412 airfoils
airfoil_file = project_dir * "/test/data/naca4412.dat"
airfoils = fill(ccb.AlphaAF(airfoil_file), length(rnondim))

#Vinf
Vinf = 5.0

# rotor rotation rate in rad/s
Omega = 5400 * pi / 30  # convert from RPM to rad/s

nwake_sheets = 15

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

initwake(inputs, inputs.blade_elements; vx_rotor=0.0, vtheta_rotor=0.0, vr_rotor=0.0)
