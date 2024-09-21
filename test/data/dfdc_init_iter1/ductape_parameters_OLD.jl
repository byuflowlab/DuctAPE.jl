#=
This file contains parameters (from the DFDC case file set) in the input format for DuctAPE
=#

using DuctAPE #for airfoil parameterization

##### ----- Operation Conditions (OPER) in DFDC Case File ----- #####
Vref = 50.0 # Vref in DFDC case file
RPM = 8000.0 # RPM in DFDC case file
rhoinf = 1.226 # Rho in DFDC case file
asound = 340.0 # Vso in DFDC case file
muinf = 1.78e-5 # Rmu in DFDC case file
wake_length = 0.8 # XDwake in DFDC case file
nwake_sheets = 11 # NRPdef in DFDC case file, also 1 more than the nstations in this case

# - DuctAPE additional parameters - #
Omega = RPM * pi / 30  # convert from RPM to rad/s for DuctAPE

nhub_inlet = 30 # chosen somewhat arbitrarily
nduct_inlet = 30 # chosen somewhat arbitrarily
npanels = [30, 1, 30] # chosen somewhat arbitrarily, the 1 is due to the fact that the duct and center body trailing edges are not quite aligned.

##### ----- Airfoil Parameters (AERO in DFDC case file) ----- #####
afparams = DuctAPE.c4b.DFDCairfoil(;
    alpha0=0.0, # A0deg in DFDC case file
    clmax=1.5, # CLmax in DFDC case file
    clmin=-1.0, # CLmin in DFDC case file
    dclda=6.28, # dCLdA in DFDC case file
    dclda_stall=0.5, # dCLdAstall in DFDC case file
    dcl_stall=0.2, # dCLstall in DFDC case file
    cdmin=0.012, # CDmin in DFDC case file
    clcdmin=0.1, # CLCDmin in DFDC case file
    dcddcl2=0.005, # dCDdCL^2 in DFDC case file
    cmcon=0.0, # Cmconst in DFDC case file
    Re_ref=2e5, # REref in DFDC case file
    Re_exp=0.35, # REexp in DFDC case file
    mcrit=0.7, # Mcrit in DFDC case file
)

##### ----- Rotor Parameters (ROTOR in DFDC case file) ----- #####

rotorzloc = 0.12 # Xdisk in DFDC case file
B = 5 # Nblds in DFDC case file

rct = [
    0.50491E-01 0.89142E-01 69.012
    0.61567E-01 0.79785E-01 59.142
    0.72644E-01 0.71300E-01 51.825
    0.83721E-01 0.63979E-01 46.272
    0.94798E-01 0.57777E-01 41.952
    0.10587 0.52541E-01 38.509
    0.11695 0.48103E-01 35.699
    0.12803 0.44316E-01 33.354
    0.13911 0.41061E-01 31.349
    0.15018 0.38243E-01 29.596
] # copied directly from matrix in DFDC case file

# DuctAPE specific parameters
# Rtip = 0.15572 # comes from DFDC output file
# Rhub = 0.0450 # comes from DFDC output file
Rtip = 0.15572081487373543 # Calculated from Geometry
Rhub = 0.04495252299071941 # Calculated from Geometry

r = rct[:, 1] ./ Rtip # non-dimensionalize the radial stations
chords = rct[:, 2]
twists = rct[:, 3] * pi / 180.0 # convert to radians

airfoils = fill(afparams, length(r)) # specify the airfoil array

##### ----- DuctAPE Input Tuples ----- #####

# - Rotor Parameters: Vector of NTuples - #
rotor = [(;
    nwake_sheets,
    rotorzloc,
    r,
    chords,
    twists,
    airfoils,
    Rtip,
    Rhub,
    tip_gap=0.0, # currently only zero tip gaps work.
    B,
    Omega,
    fliplift=false, # can flip the cl lookups on the fly if desired, say, for stator sections
)]

# - Paneling Constants - #
paneling_constants = (; npanels, nhub_inlet, nduct_inlet, wake_length, nwake_sheets)

# - Freestream (will need to update Vinf when running sweep of advance ratios - #
# this will need to be redefined in the loop

D = 2.0 * Rtip # rotor diameter
n = RPM / 60.0 # rotation rate in revolutions per second
Vinf = 1.0 * n * D
freestream = (; rhoinf, muinf, asound, Vinf)

# - Reference Parameters used for Post-processing - #
reference_parameters = (; Vref, Rref=Rtip)

# ##### ----- Solid Body Coordinates (GEOM in the DFDC case file) ----- #####
# # copied from the DFDC case file and put here at end since they're so long
# # also note that they need to be reversed for DuctAPE

# hub_coordinates = reverse(
#     [
#         0.306379 0.035928
#         0.296701 0.039168
#         0.285768 0.041795
#         0.274184 0.043576
#         0.262095 0.044564
#         0.248578 0.044839
#         0.232723 0.044847
#         0.216580 0.044866
#         0.200407 0.044882
#         0.184232 0.044898
#         0.168060 0.044913
#         0.151899 0.044930
#         0.135773 0.044950
#         0.120000 0.044952
#         0.109361 0.044999
#         0.099680 0.044922
#         0.090039 0.044561
#         0.080630 0.043916
#         0.071451 0.042966
#         0.062540 0.041690
#         0.053933 0.040075
#         0.045705 0.038107
#         0.037900 0.035771
#         0.030604 0.033068
#         0.023901 0.030001
#         0.017880 0.026585
#         0.012632 0.022848
#         0.008231 0.018825
#         0.004736 0.014551
#         0.002179 0.010047
#         0.000586 0.005293
#         0.000000 0.000000
#     ];
#     dims=1,
# )

# duct_coordinates = reverse(
#     [
#         0.304542 0.159526
#         0.296440 0.162842
#         0.285240 0.167270
#         0.273301 0.171796
#         0.261206 0.176188
#         0.249026 0.180416
#         0.236794 0.184470
#         0.224549 0.188335
#         0.212269 0.192017
#         0.200020 0.195490
#         0.187807 0.198746
#         0.175647 0.201771
#         0.163586 0.204547
#         0.151644 0.207051
#         0.139859 0.209267
#         0.128380 0.211153
#         0.117508 0.212618
#         0.106950 0.213649
#         0.096610 0.214255
#         0.086499 0.214421
#         0.076647 0.214136
#         0.067113 0.213390
#         0.057935 0.212174
#         0.049193 0.210487
#         0.040949 0.208332
#         0.033312 0.205736
#         0.026403 0.202735
#         0.020339 0.199393
#         0.015227 0.195790
#         0.011135 0.192004
#         0.008090 0.188086
#         0.006112 0.184067
#         0.005242 0.180005
#         0.005484 0.176154
#         0.006854 0.172546
#         0.009324 0.169289
#         0.012842 0.166404
#         0.017419 0.163862
#         0.023110 0.161648
#         0.029956 0.159771
#         0.037937 0.158256
#         0.046983 0.157103
#         0.057025 0.156294
#         0.067995 0.155792
#         0.079836 0.155546
#         0.092531 0.155498
#         0.106044 0.155585
#         0.120000 0.155721
#         0.134222 0.155902
#         0.148679 0.156177
#         0.163490 0.156523
#         0.178507 0.156897
#         0.193399 0.157258
#         0.208123 0.157586
#         0.222751 0.157864
#         0.237332 0.158088
#         0.251898 0.158254
#         0.266505 0.158365
#         0.281130 0.158423
#         0.294972 0.158441
#         0.304466 0.158439
#     ];
#     dims=1,
# )
