r1 = [0.25; 0.5; 0.75; 1.0]
Rtip = [1.0, 1.0]
rnondim1 = r1 ./ Rtip[1]
rnondim = [rnondim1 rnondim1]
afparams1 = dt.c4b.DFDCairfoil(
    0.0, 1.5, -1.0, 6.28, 0.5, 0.2, 0.012, 0.1, 0.005, 0.0, 200000.0, 0.35, 0.7
)
rotorzloc = [0.25, 0.75]
r = rnondim
chords = 0.1 * ones(size(rnondim))
twists = 20.0 * pi / 180.0 * ones(size(rnondim))
airfoils = fill(afparams1, 4, 2)
Rhub = [0.25, 0.25]
Rtip = Rtip
tip_gap = [0.0, 0.0]
B = [2, 4]
fliplift = [0.0, 0.0]

rotorstator_parameters = dt.RotorStatorParameters(
    B, rotorzloc, r, Rhub, Rtip, chords, twists, tip_gap, airfoils, fliplift
)

ncenterbody_inlet = 1
nduct_inlet = 1
nwake_sheets = 3
wake_length = 1.0
npanels = [2, 1, 4]
dte_minus_cbte = 0

paneling_constants = dt.PanelingConstants(
    nduct_inlet, ncenterbody_inlet, npanels, dte_minus_cbte, nwake_sheets, wake_length
)

Vinf = [10.0]
rhoinf = [1.226]
muinf = [1.78e-5]
asound = [340.0]
Omega = [5000.0, 0.0] * pi / 30  # convert from RPM to rad/s

operating_point = dt.OperatingPoint(Vinf, rhoinf, muinf, asound, Omega)

Vref = [10.0]
Rref = [Rtip]
reference_parameters = dt.ReferenceParameters(Vref, Rref)

duct_coordinates = [1.0 2.0; 0.5 1.5; 0.0 2.0; 0.5 2.5; 1.0 2.0]
centerbody_coordinates = [0.0 0.0; 0.5 0.5; 1.0 0.0]

propulsor = dt.Propulsor(
    duct_coordinates,
    centerbody_coordinates,
    rotorstator_parameters,
    operating_point,
    paneling_constants,
    reference_parameters,
)
#
#
#
#
#
# plot(
#     duct_coordinates[:, 1],
#     duct_coordinates[:, 2] .- 0.75;
#     aspectratio=1,
#     marker=true,
#     label="",
# )
# plot!(hub_coordinates[:, 1], hub_coordinates[:, 2]; marker=true, label="")

# for i in 1:length(rotorzloc)
#     plot!([rotorzloc[i], rotorzloc[i]], [Rhub, Rtip]; marker=true, label="", color=3)
# end
# plot!()
