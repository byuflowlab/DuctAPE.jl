r1 = [0.25; 0.5; 0.75; 1.0]
Rtip = [1.0, 1.0]
rnondim1 = r1 ./ Rtip[1]
rnondim = [rnondim1 rnondim1]
afparams1 = dt.c4b.DFDCairfoil()
rotor_axial_position = [0.25, 0.75]
r = rnondim
chords = 0.1 * ones(size(rnondim))
twists = 20.0 * pi / 180.0 * ones(size(rnondim))
airfoils = [fill(afparams1, 4), fill(afparams1, 4)]
Rhub = [0.25, 0.25]
Rtip = Rtip
tip_gap = [0.0, 0.0]
B = [2, 4]
is_stator = [0.0, 0.0]

rotor = dt.Rotor(B, rotor_axial_position, r, Rhub, Rtip, chords, twists, tip_gap, airfoils, is_stator)

num_center_body_inlet_panels = 1
num_duct_inlet_panels = 1
num_wake_sheets = 3
wake_length = 1.0
num_panels = [2, 1, 4]
dte_minus_cbte = 0

paneling_constants = dt.PanelingConstants(
    num_duct_inlet_panels, num_center_body_inlet_panels, num_panels, dte_minus_cbte, num_wake_sheets, wake_length
)

Vinf = [10.0]
rhoinf = [1.226]
muinf = [1.78e-5]
asound = [340.0]
Omega = [5000.0, 0.0] * pi / 30  # convert from RPM to rad/s

operating_point = dt.OperatingPoint(Vinf, Omega, rhoinf, muinf, asound)

Vref = [10.0]
Rref = [Rtip]
reference_parameters = dt.ReferenceParameters(Vref, Rref)

duct_coordinates = [1.0 2.0; 0.5 1.5; 0.0 2.0; 0.5 2.5; 1.0 2.0]
center_body_coordinates = [0.0 0.0; 0.5 0.5; 1.0 0.0]

ducted_rotor = dt.DuctedRotor(
    duct_coordinates, center_body_coordinates, rotor, paneling_constants
)
