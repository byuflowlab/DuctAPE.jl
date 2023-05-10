
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

include(project_dir * "/plots_default.jl")

using DuctTAPE
const dt = DuctTAPE

# CCBlade used for it's airfoils function objects here.
using CCBlade
const ccb = CCBlade

using FLOWMath
const fm = FLOWMath

#---------------------------------#
#         Rotor Geometry          #
#---------------------------------#
#=
ROTOR
!       Xdisk        Nblds       NRPdef
  0.12000                5           11
!  #stations
    10
!           r        Chord         Beta
  0.50491E-01  0.89142E-01   69.012
  0.61567E-01  0.79785E-01   59.142
  0.72644E-01  0.71300E-01   51.825
  0.83721E-01  0.63979E-01   46.272
  0.94798E-01  0.57777E-01   41.952
  0.10587      0.52541E-01   38.509
  0.11695      0.48103E-01   35.699
  0.12803      0.44316E-01   33.354
  0.13911      0.41061E-01   31.349
  0.15018      0.38243E-01   29.596
ENDROTOR
=#

xrotor = 0.12
B = 5

#dimensional radius
r = [
    0.50491E-01
    0.61567E-01
    0.72644E-01
    0.83721E-01
    0.94798E-01
    0.10587
    0.11695
    0.12803
    0.13911
    0.15018
]

#dimensional chord
chords = [
    0.89142E-01
    0.79785E-01
    0.71300E-01
    0.63979E-01
    0.57777E-01
    0.52541E-01
    0.48103E-01
    0.44316E-01
    0.41061E-01
    0.38243E-01
]

#twist in degrees converted to radians
twists =
    [
        69.012
        59.142
        51.825
        46.272
        41.952
        38.509
        35.699
        33.354
        31.349
        29.596
    ] * pi / 180.0

airfoils = fill(ccb.AlphaAF("examples/dfdc_comp/dfdc_af_test1.dat"), length(r))

#---------------------------------#
#       Duct and Hub Geometry     #
#---------------------------------#

include(project_dir * "/examples/dfdc_comp/dfdc_ductgeom.jl")
_, duct_leidx = findmin(duct_coords[:, 1])
ductxin = reverse(duct_coords[1:duct_leidx, 1])
ductrin = reverse(duct_coords[1:duct_leidx, 2])

# load in duct and hub geometry, spline, and find out what the duct and hub radii are at the rotor positions to figure out what Rtip and Rhub are.
Rhub = FLOWMath.akima(hub_coords[:, 1], hub_coords[:, 2], xrotor)
Rtip = FLOWMath.akima(ductxin, ductrin, xrotor)

#---------------------------------#
#      Operating Conditions       #
#---------------------------------#
#=
OPER
!        Vinf         Vref          RPM
   0.0000       50.000       8000.0
!         Rho          Vso          Rmu           Alt
   1.2260       340.00      0.17800E-04   0.0000
!       XDwake        Nwake
  0.80000               20
!       Lwkrlx
            F
ENDOPER
=#

Vinf = 0.0 #TODO: probably want to change that...
Vref = 50.0 #TODO: look up how this is used.  is this why the Cp values are low in the dfdc example???
RPM = 8000
Omega = RPM * pi / 30  # convert from RPM to rad/s
rho = 1.226
asound = 340.0
mu = 1.78e-5

#---------------------------------#
#        Paneling Options         #
#---------------------------------#

wake_length = 0.8 #times duct chord?

nwake_sheets = 11 #note that nwake in the dfdc file is how many panels to have in the wake

npanels_inlet = 10
discscale = 1.0

ductle = minimum(duct_coords[:, 1])
ductte = maximum(duct_coords[:, 1])
ductchord = maximum(duct_coords[:, 1]) - minimum(duct_coords[:, 1])
outletinletratio = (ductte - xrotor) / (xrotor - ductle)

nhub_inlet = round(Int, npanels_inlet * discscale)

nduct_inlet = round(Int, npanels_inlet * discscale)

nduct_outlet = round(Int, nduct_inlet * outletinletratio)

nwake = round(Int, (nduct_inlet + nduct_outlet) * wake_length)

npanels = [nduct_outlet, nwake]

#---------------------------------#
#      Assemble Named Tuples      #
#---------------------------------#

# Rotor Parameters
rotor_parameters = [(;
    xrotor=xrotor,
    nwake_sheets,
    r=r ./ Rtip, #non-dimensionalize
    chords,
    twists,
    airfoils,
    Rtip,
    Rhub,
    tip_gap=0.0,
    B,
    Omega,
)]

# Paneling Parameters
paneling_constants = (; npanels, nhub_inlet, nduct_inlet, wake_length, nwake_sheets)

# Freestream Parameters
freestream = (; rho, mu, asound, Vinf)

#---------------------------------#
#           Plot Inputs           #
#---------------------------------#

pgeom = plot(; aspectratio=1, xlabel="x", rlabel="r")
plot!(
    pgeom,
    duct_coords[:, 1],
    duct_coords[:, 2];
    linewidth=1.0,
    marker=true,
    markersize=2.0,
    label="Duct Input Coordinates",
)
plot!(
    pgeom,
    hub_coords[:, 1],
    hub_coords[:, 2];
    linewidth=1.0,
    marker=true,
    markersize=2.0,
    label="Hub Input Coordinates",
)
plot!(
    pgeom,
    xrotor * ones(length(r)),
    r;
    linewidth=1.0,
    marker=true,
    markersize=2.0,
    label="Rotor Input Coordinates",
)

#---------------------------------#
#           Run Solver            #
#---------------------------------#

strengths, inputs, initials, convergeflag = dt.analyze_propulsor(
    duct_coords,
    hub_coords,
    paneling_constants,
    rotor_parameters,
    freestream;
    debug=false,
    maximum_linesearch_step_size=1e6,
    iteration_limit=100,
)

println("solution converged? ", convergeflag)
if convergeflag
    convlabel = "Converged"
else
    convlabel = "NOT converged"
end

#---------------------------------#
#    Extract Outputs and Plot     #
#---------------------------------#
gamb, gamw, Gamr, sigr = dt.extract_state_variables(strengths, inputs)

##### ----- Plot GEOMETRY ----- #####
#initialize plot
pgeom = plot(; aspectratio=1, xlabel="x", ylabel="r")
plot!(
    pgeom,
    xrotor * ones(length(inputs.rotor_panel_edges)),
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

savefig(project_dir * "/examples/dfdc_comp/precomputed-rotor-wake-geometry.pdf")

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

plot!(
    pgeom,
    inputs.hub_coordinates[:, 1],
    inputs.hub_coordinates[:, 2];
    linewidth=0.25,
    markersize=0.5,
    markershape=:rect,
    color=mycolors[3],
    label="",
)

# Plot body panel centers
for ib in 1:2
    # for ib in 1:1
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

savefig(pgeom, project_dir * "/examples/dfdc_comp/precomputed-full-geometry.pdf")

##### ----- Plot rotor circulation distribution ----- #####
# initialize plot
pG = plot(; xlabel=L"\Gamma", ylabel="r")

# plot solution
# TODO: figure out where to add a print statement in DFDC to compare this output
plot!(pG, Gamr, inputs.rotor_panel_centers; label=convlabel)

# # plot initials
# plot!(
#     pG,
#     Gamr_init,
#     inputs.rotor_panel_centers;
#     xlabel=L"\Gamma",
#     ylabel="r",
#     label="Initial",
#     color=mycolors[1],
#     linestyle=:dash,
# )

#save
savefig(pG, project_dir * "/examples/dfdc_comp/rotorcirculation-check.pdf")

##### ----- Plot duct surface velocity ----- #####

#prepare outputs
dp = inputs.body_panels[1].panel_center[:, 1]
_, leidx = findmin(dp)
#split into inner and outer surfaces
dpinner = dp[1:leidx]
dpouter = dp[(leidx + 1):end]

#TODO: NOTE THAT YOU'RE USING VREF HERE SEE NOTES ABOUT IT ABOVE.
gamdinner = gamb[1:leidx] / Vref
gamdouter = gamb[(leidx + 1):length(dp)] / Vref

# initialize plot
pb = plot(; xlabel="x", ylabel="Q/Qref")

# plot solution
plot!(pb, dpinner, gamdinner; label=convlabel * " inner surface, with rotor")

plot!(pb, dpouter, gamdouter; label=convlabel * " outer surface, with rotor")

###prepare initials
##gamd_init = 1.0 .- (gamb_init[1:length(dp)] ./ Vinf) .^ 2
##gamd_init = gamb_init[1:length(dp)] ./ Vinf
###split into inner and outer surfaces
##gamdinner_init = gamd_init[1:leidx]
##gamdouter_init = gamd_init[(leidx + 1):end]

# plot initials
##plot!(pb, dpinner, gamdinner_init; label="initial inner surface, with rotor")
##plot!(pb, dpouter, gamdouter_init; label="initial outer surface, with rotor")

## plot isolated duct
## spit inner/outer
#_, iductle = findmin(ductouts.x)
#idxin = ductouts.x[1:iductle]
#idxout = ductouts.x[(iductle + 1):end]
#idvsin = ductouts.vs_duct[1:iductle]
#idvsout = ductouts.vs_duct[(iductle + 1):end]

#plot!(
#    pb,
#    idxin,
#    idvsin;
#    label="isolated duct inner surface (initial state)",
#    color=mycolors[1],
#    linestyle=:dash,
#)

#plot!(
#    pb,
#    idxout,
#    idvsout;
#    label="isolated duct outer surface (initial state)",
#    color=mycolors[2],
#    linestyle=:dash,
#)

#plot rotor location
plot!(
    pb,
    xrotor * ones(2),
    [minimum([gamdinner; gamdouter]); maximum([gamdinner; gamdouter])];
    # linewidth=0.25,
    linestyle=:dash,
    color=mycolors[3],
    label="rotor location",
)

#save
savefig(pb, project_dir * "/examples/dfdc_comp/body-velocity.pdf")

##### ----- Plot Surface Pressure ----- #####
cpd = 1.0 .- (gamb[1:length(dp)] ./ Vref) .^ 2
cpdinner = cpd[1:leidx]
cpdouter = cpd[(leidx + 1):length(dp)]

pcp = plot(; xlabel="x", ylabel=L"C_p")

# plot solution
plot!(pcp, dpinner, cpdinner; label=convlabel * " inner surface, with rotor")

plot!(pcp, dpouter, cpdouter; label=convlabel * " outer surface, with rotor")

#plot rotor location
plot!(
    pcp,
    xrotor * ones(2),
    [minimum([cpdinner; cpdouter]); maximum([cpdinner; cpdouter])];
    # linewidth=0.25,
    linestyle=:dash,
    color=mycolors[3],
    label="rotor location",
)

#save
savefig(pcp, project_dir * "/examples/dfdc_comp/body-pressure.pdf")

##### ----- Plot Wake Strengths ----- #####
pg = plot(; xlabel=L"\gamma_\theta^{wake}", ylabel="r")

# plot solution
plot!(pg, gamw, inputs.rotor_panel_edges; label=convlabel)

## plot initial
#plot!(
#    pg,
#    gamw_init,
#    inputs.rotor_panel_edges;
#    label="initial",
#    linestyle=:dash,
#    color=mycolors[1],
#)

#save
savefig(pg, project_dir * "/examples/dfdc_comp/wake-strength.pdf")

##### ----- Plot Source Strengths ----- #####
ps = plot(; xlabel=L"\sigma", ylabel="r")

# plot solution
plot!(ps, sigr, inputs.rotor_panel_centers; label=convlabel)

# # plot initial
# plot!(
#     ps,
#     sigr_init,
#     inputs.rotor_panel_centers;
#     label="initial",
#     linestyle=:dash,
#     color=mycolors[1],
# )

#save
savefig(ps, project_dir * "/examples/dfdc_comp/source-strength-check.pdf")

