#=

Example script for running the rotor-only case.

Rotor data comes from the APC 10x5 data available on the UIUC database.

Authors: Judd Mehr, Andrew Ning

=#

#---------------------------------#
#             Includes            #
#---------------------------------#
using DuctTAPE
const dt = DuctTAPE

# CCBlade used for it's airfoils function objects here.
using CCBlade
const ccb = CCBlade
include("run_ccblade.jl")

# using Plots
# pyplot()
# using LaTeXStrings
include("../../plots_default.jl")

#---------------------------------#
#             Geometry            #
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
r = propgeom[:, 1]
# Dimensionalize chords
chords = propgeom[:, 2] * Rtip
# convert twists to radians
twists = propgeom[:, 3] * pi / 180

# plot(r, chords, xlabel=L"r/R_\mathrm{tip}", ylabel=L"c/R_\mathrm{tip}", label="",ylim=(0.0,0.22),xlim=(0.0,1.0))
# savefig("examples/rotor_only/apc_chord.pdf")
# plot(r, twists, xlabel=L"r/R_\mathrm{tip}", ylabel="twist (deg)",label="",ylim=(0.0,39),xlim=(0.0,1.0))
# savefig("examples/rotor_only/apc_twist.pdf")

# use a NACA 4412 airfoils
#=
Note here we are using the CCBlade functionality to define the airfoils data function.
In addition, we are using the airfoils data file available from the CCBlade repository that has been extrapolated using the Viterna method as well as corrected for rotational effects as described in the CCBlade documentation.
=#
airfoils = fill(ccb.AlphaAF("test/data/naca4412.dat"), length(r))

##### ----- User Options ----- #####
# number of blade elements to use in analysis
#=
Note: the solver with interpolate the rotor data using the given number of blade element inputs
=#
nwake_sheets = 15

# x position of rotor
xrotor = 0.0

#---------------------------------#
#       Operating Conditions      #
#---------------------------------#

# rotor rotation rate in rad/s
Omega = 5400 * pi / 30  # convert from RPM to rad/s

# freestream conditions
rho = 1.225 #kg/m^3
mu = 1.81e-5 # kg/(m⋅s)
asound = 341.0 #m/s

#---------------------------------#
#          Define Inputs          #
#---------------------------------#

# Rotor Parameters
rotor_parameters = [(;
    xrotor, nwake_sheets, r, chords, twists, airfoils, Rtip, Rhub, B, Omega
)]

# Freestream Parameters
Vinf = 5.0
freestream = (; rho, mu, asound, Vinf)

#---------------------------------#
#        Experimental Data        #
#---------------------------------#
exp = [
    0.113 0.0912 0.0381 0.271
    0.145 0.0890 0.0386 0.335
    0.174 0.0864 0.0389 0.387
    0.200 0.0834 0.0389 0.429
    0.233 0.0786 0.0387 0.474
    0.260 0.0734 0.0378 0.505
    0.291 0.0662 0.0360 0.536
    0.316 0.0612 0.0347 0.557
    0.346 0.0543 0.0323 0.580
    0.375 0.0489 0.0305 0.603
    0.401 0.0451 0.0291 0.620
    0.432 0.0401 0.0272 0.635
    0.466 0.0345 0.0250 0.644
    0.493 0.0297 0.0229 0.640
    0.519 0.0254 0.0210 0.630
    0.548 0.0204 0.0188 0.595
    0.581 0.0145 0.0162 0.520
]
# extract advance ratios
Jexp = exp[:, 1]
# extract thrust coefficients
CTexp = exp[:, 2]
# extract power coefficients
CPexp = exp[:, 3]
# extract efficiencies
etaexp = exp[:, 4]

#---------------------------------#
#          Set Up Solves          #
#---------------------------------#
# get values needed for backing out freestream velocity from advance ratio
n = Omega / (2 * pi) #get revolutions per second
D = 2 * Rtip #rotor tip diameter

J = collect(range(0.1, 0.6; step=0.025))  # advance ratio
# J = Vinf / (D * n)
nJ = length(J)

# Initialize Inputs
wake_length=3.0
inputs, params = dt.initialize_rotor_states(rotor_parameters, freestream, wake_length)

Gamr, gamw, sigr = dt.extract_rotor_states(inputs, params)

rbe = params.blade_elements[1].rbe
# nwake_sheets = length(rbe)

# initialize outputs
eff = zeros(nJ)
CT = zeros(nJ)
CQ = zeros(nJ)
effccb = zeros(nJ)
CTccb = zeros(nJ)
CQccb = zeros(nJ)

# Loop through advance ratios
for i in 1:nJ
    println()
    println("Running for J = $(J[i])")
    println()

    # calculate freestream velocity for given advance ratio
    Vinf_sweep = J[i] * D * n

    @time begin
        # run solver
        # remember to update freestream velocity in parameters
        states = dt.solve_rotor_only(inputs, (; params..., Vinf=Vinf_sweep))
        Gamrconv, gamwconv, sigrconv = dt.extract_rotor_states(states, params)

        # pG = plot(
        #     Gamr, params.rotor_panel_centers; xlabel=L"\Gamma", ylabel="r", label="initial"
        # )
        # plot!(pG, Gamrconv, params.rotor_panel_centers; label="converged")
        # savefig("examples/rotor_only/Circulation_J$(J[i]).pdf")

        # pw = plot(
        #     gamw,
        #     params.rotor_panel_edges;
        #     xlabel=L"\gamma_\theta",
        #     ylabel="r",
        #     label="initial",
        # )
        # plot!(pw, gamwconv, params.rotor_panel_edges; label="converged")
        # savefig("examples/rotor_only/wake-strengths_J$(J[i]).pdf")

        # ps = plot(
        #     sigr, params.rotor_panel_centers; xlabel=L"\sigma", ylabel="r", label="initial"
        # )
        # plot!(ps, sigrconv, params.rotor_panel_centers; label="converged")
        # savefig("examples/rotor_only/source-strenghts_J$(J[i]).pdf")

        # - Post Process - #
        dtout = dt.states_to_outputs_rotor_only(states, (; params..., Vinf=Vinf_sweep))

        aero = dt.get_rotor_loads(
            dtout.W,
            dtout.phi,
            dtout.cl,
            dtout.cd,
            params.blade_elements[1],
            (; freestream..., Vinf=Vinf_sweep),
        )

        # states, dtout, aero = dt.analyze_propulsor(rotor_parameters, (; freestream..., Vinf=Vinf_sweep))
    end

    CT[i] = aero.CT
    CQ[i] = aero.CQ
    eff[i] = aero.eff

    # Run CCBlade:
    ccbouts = run_ccblade(Vinf_sweep)
    effccb[i] = ccbouts.eff
    CTccb[i] = ccbouts.CT
    CQccb[i] = ccbouts.CQ
    out = ccbouts.out

    #### --- PLOTS --- ###
    ##Uncomment to see all the details
    #println()
    #println("Plotting...")
    #println()

    ## double check geometry
    #plot(
    #    dtout.r,
    #    dtout.chord;
    #    xlabel="r",
    #    ylabel="chords",
    #    label="DuctTAPE",
    #    title="J = $(J[i])",
    #)
    #plot!(r * Rtip, chords; label="CCBlade")
    #savefig("examples/rotor_only/chord_J$(J[i]).pdf")

    #plot(
    #    dtout.r,
    #    dtout.twist * 180 / pi;
    #    ylabel="twists (deg)",
    #    xlabel="r",
    #    label="DuctTAPE",
    #    title="J = $(J[i])",
    #)
    #plot!(r * Rtip, twists * 180 / pi; label="CCBlade")
    #savefig("examples/rotor_only/twist_J$(J[i]).pdf")

    ## axial induced velocity
    #plot(
    #    dtout.vx_rotor,
    #    rbe / Rtip;
    #    xlabel=L"v_x",
    #    ylabel="r",
    #    label="DuctTAPE",
    #    title="J = $(J[i])",
    #)
    #plot!(out.u, r; label="CCBlade")
    #savefig("examples/rotor_only/vx_J$(J[i]).pdf")

    ## tangential induced velocity
    #plot(
    #    dtout.vtheta_rotor,
    #    rbe / Rtip;
    #    xlabel=L"v_\theta",
    #    ylabel="r",
    #    label="DuctTAPE",
    #    title="J = $(J[i])",
    #)
    #plot!(out.v, r; label="CCBlade")
    #savefig("examples/rotor_only/vtheta_J$(J[i]).pdf")

    ## tangential total velocity
    #plot(
    #    dtout.Wθ,
    #    rbe / Rtip;
    #    xlabel=L"W_\theta",
    #    ylabel="r",
    #    label="DuctTAPE",
    #    title="J = $(J[i])",
    #)
    #plot!(out.v .- Omega * r * Rtip, r; label="CCBlade")
    #savefig("examples/rotor_only/Wtheta_J$(J[i]).pdf")

    ## meridional total velocity
    #plot(
    #    dtout.Wm,
    #    rbe / Rtip;
    #    xlabel=L"W_m",
    #    ylabel="r",
    #    label="DuctTAPE",
    #    title="J = $(J[i])",
    #)
    #plot!(out.u .+ Vinf_sweep, r; label="CCBlade")
    #savefig("examples/rotor_only/Wm_J$(J[i]).pdf")

    ## inflow angle
    #plot(
    #    dtout.phi * 180 / pi,
    #    rbe / Rtip;
    #    xlabel=L"\phi~(deg)",
    #    ylabel="r",
    #    label="DuctTAPE",
    #    title="J = $(J[i])",
    #)
    #plot!(out.phi * 180 / pi, r; label="CCBlade")
    #savefig("examples/rotor_only/phi_J$(J[i]).pdf")

    ## angle of attack
    #plot(
    #    dtout.alpha * 180 / pi,
    #    rbe / Rtip;
    #    xlabel=L"\alpha~(deg)",
    #    ylabel="r",
    #    label="DuctTAPE",
    #    title="J = $(J[i])",
    #)
    #plot!(out.alpha * 180 / pi, r; label="CCBlade")
    #savefig("examples/rotor_only/alpha_J$(J[i]).pdf")

    ## inflow magnitude
    #plot(
    #    dtout.W, rbe / Rtip; xlabel=L"W", ylabel="r", label="DuctTAPE", title="J = $(J[i])"
    #)
    #plot!(out.W, r; label="CCBlade")
    #savefig("examples/rotor_only/W_J$(J[i]).pdf")

    ## Lift
    #plot(
    #    dtout.cl,
    #    rbe / Rtip;
    #    xlabel=L"c_\ell",
    #    ylabel="r",
    #    label="DuctTAPE",
    #    title="J = $(J[i])",
    #)
    #plot!(dtout.clin, rbe / Rtip; label="inner", linestyle=:dash, color=mycolors[1])
    #plot!(dtout.clout, rbe / Rtip; label="outer", linestyle=:dot, color=mycolors[1])
    #plot!(out.cl, r; label="CCBlade")
    #savefig("examples/rotor_only/cl_J$(J[i]).pdf")

    ## Drag
    #plot(
    #    dtout.cd,
    #    rbe / Rtip;
    #    xlabel=L"c_d",
    #    ylabel="r",
    #    label="DuctTAPE",
    #    title="J = $(J[i])",
    #)
    #plot!(dtout.cdin, rbe / Rtip; label="inner", linestyle=:dash, color=mycolors[1])
    #plot!(dtout.cdout, rbe / Rtip; label="outer", linestyle=:dot, color=mycolors[1])
    #plot!(out.cd, r; label="CCBlade")
    #savefig("examples/rotor_only/cd_J$(J[i]).pdf")

    ## normal coeff
    #plot(
    #    aero.cn,
    #    rbe / Rtip;
    #    xlabel=L"c_n",
    #    ylabel="r",
    #    label="DuctTAPE",
    #    title="J = $(J[i])",
    #)
    #plot!(out.cn, r; label="CCBlade")
    #savefig("examples/rotor_only/cn_J$(J[i]).pdf")

    ## tangential coeff
    #plot(
    #    aero.ct,
    #    rbe / Rtip;
    #    xlabel=L"c_t",
    #    ylabel="r",
    #    label="DuctTAPE",
    #    title="J = $(J[i])",
    #)
    #plot!(out.ct, r; label="CCBlade")
    #savefig("examples/rotor_only/ct_J$(J[i]).pdf")

    ## Distributed Loads
    #plot(
    #    rbe / Rtip,
    #    aero.Np;
    #    xlabel="r/Rtip",
    #    label="Normal load per unit length",
    #    title="J = $(J[i])",
    #)
    #plot!(r, out.Np; label="CCBlade Np")
    #plot!(rbe / Rtip, aero.Tp; label="Tangential load per unit length")
    #plot!(r, out.Tp; label="CCBlade Tp")
    #savefig("examples/rotor_only/distributed_loads_J$(J[i]).pdf")
end

#---------------------------------#
#              Plots              #
#---------------------------------#

plot(J, CT; xlabel=L"J", label=L"C_T~DuctTAPE", color=mycolors[1])
plot!(J, CQ * 2 * pi; label=L"C_P~DuctTAPE", color=mycolors[2])
plot!(J, CTccb; label=L"C_T~CCBlade", color=mycolors[1], linestyle=:dash)
plot!(J, CQccb * 2 * pi; label=L"C_P~CCBlade", color=mycolors[2], linestyle=:dash)
plot!(Jexp, CTexp; seriestype=:scatter, label="C_T~experimental", color=mycolors[1])
plot!(Jexp, CPexp; seriestype=:scatter, label="C_P~experimental", color=mycolors[2])
savefig("examples/rotor_only/rotor-only-thrust-and-power-validation.pdf")
# savefig("examples/rotor_only/rotor-only-thrust-and-power-validation.png")

plot(J, eff; xlabel=L"J", ylabel=L"\eta", label="DuctTAPE")
plot!(J, effccb; label="CCBlade")
plot!(Jexp, etaexp; seriestype=:scatter, label="experimental")
savefig("examples/rotor_only/rotor-only-efficiency-validation.pdf")
# savefig("examples/rotor_only/rotor-only-efficiency-validation.png")
