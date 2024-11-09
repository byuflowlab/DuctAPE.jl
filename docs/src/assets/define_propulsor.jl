afparams = DuctAPE.c4b.DFDCairfoil(;
    alpha0=0.0,
    clmax=1.5,
    clmin=-1.0,
    dclda=6.28,
    dclda_stall=0.5,
    dcl_stall=0.2,
    cdmin=0.012,
    clcdmin=0.1,
    dcddcl2=0.005,
    cmcon=0.0,
    Re_ref=2e5,
    Re_exp=0.35,
    mcrit=0.7,
)

airfoils = fill(afparams, length(r)) # specify the airfoil array

rotor = dt.Rotor(
    [5],
    [rotorzloc],
    r,
    [Rhub+0.01],
    [Rtip-0.025],
    c,
    t,
    [0.0], # currently only zero tip gaps work.
    airfoils,
    [0.0], # can flip the cl lookups on the fly if desired, say, for stator sections
)

# Freestream
Vinf = 0.0 # hover condition
rhoinf = 1.226
asound = 340.0
muinf = 1.78e-5

# Rotation Rate
RPM = 8000.0
Omega = RPM * pi / 30 # if using RPM, be sure to convert to rad/s

# utilizing the constructor function to put things in vector types
operating_point = dt.OperatingPoint(Vinf, Omega, rhoinf, muinf, asound)

nduct_inlet = 50
ncenterbody_inlet = 30
npanels = [50, 10, 30] # the 1 is due to the fact that the duct and center body trailing edges are not quite aligned.
dte_minus_cbte = 1.0 # the duct trailing edge is ahead of the centerbody trailing edge.
nwake_sheets = 22
wake_length = 0.2

paneling_constants = dt.PanelingConstants(
    nduct_inlet, ncenterbody_inlet, npanels, dte_minus_cbte, nwake_sheets, wake_length
)

Vref = 50.0 #this turns out to be close to the average axial velocity at the rotor in our case
Rref = Rtip

reference_parameters = dt.ReferenceParameters([Vref], [Rref])

dz = [reverse(cz); nz[2:end]]
dr = [reverse(cr); nr[2:end]]

ducted_rotor = dt.DuctedRotor(
    [dz dr],
    [cbz cbr],
    rotor,
    operating_point,
    paneling_constants,
    reference_parameters,
)
