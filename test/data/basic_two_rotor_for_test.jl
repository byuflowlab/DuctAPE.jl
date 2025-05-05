rotorzloc = [0.25, 0.75]
B = [2, 4]
r = [0.25; 0.5; 0.75; 1.0]
r = [r r]
chords = 0.1 * ones(size(r))
twists = 20.0 * pi / 180.0 * ones(size(r))
afparams1 = dt.c4b.DFDCairfoil(
    0.0, 1.5, -1.0, 6.28, 0.5, 0.2, 0.012, 0.1, 0.005, 0.0, 200000.0, 0.35, 0.7
)
airfoils1 = fill(afparams1, length(r[:, 1]))
airfoils2 = fill(afparams1, length(r[:, 1]))
airfoils = [airfoils1 airfoils2]
Rtip = 1.0
Rhub = 0.25
rhoinf = 1.226
muinf = 1.78e-5
Vinf = 10.0
Vref = 10.0
Omega = [5000.0, 0.0] * pi / 30  # convert from RPM to rad/s
asound = 340.0
nhub_inlet = 1
nduct_inlet = 1
nwake_sheets = 3
wake_length = 1.0
npanels = [2, 1, 4]
rotor = [
    (;
        rotorzloc=rotorzloc[i],
        nwake_sheets,
        r=r[:, i] ./ Rtip, #non-dimensionalize
        chords=chords[:, i],
        twists=twists[:, i],
        airfoils=airfoils[:, i],
        Rtip,
        Rhub,
        tip_gap=0.0,
        B=B[i],
        Omega=Omega[i],
        is_stator=false,
    ) for i in 1:length(Omega)
]
paneling_constants = (; npanels, nhub_inlet, nduct_inlet, wake_length, nwake_sheets)
freestream = (; rhoinf, muinf, asound, Vinf)
reference_parameters = (; Vref, Rref=Rtip)
duct_coordinates = [1.0 2.0; 0.5 1.5; 0.0 2.0; 0.5 2.5; 1.0 2.0]
hub_coordinates = [0.0 0.0; 0.5 0.5; 1.0 0.0]
