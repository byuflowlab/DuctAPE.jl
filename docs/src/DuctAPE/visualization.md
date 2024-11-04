# Visualization

There are several convenience plotting methods implemented in DuctAPE based on [RecipesBase](https://juliaplots.org/RecipesBase.jl/stable/).
In addition a general function for plotting the suite of available plots or animations is provided in the `generate_plots` function.

```@docs; canonical=false
DuctAPE.generate_plots
```

```@setup visualize
using DuctAPE
include("../assets/plots_default.jl")
gr()

duct_coordinates = [
    0.304466 0.158439
    0.294972 0.158441
    0.28113 0.158423
    0.266505 0.158365
    0.251898 0.158254
    0.237332 0.158088
    0.222751 0.157864
    0.208123 0.157586
    0.193399 0.157258
    0.178507 0.156897
    0.16349 0.156523
    0.148679 0.156177
    0.134222 0.155902
    0.12 0.155721
    0.106044 0.155585
    0.092531 0.155498
    0.079836 0.155546
    0.067995 0.155792
    0.057025 0.156294
    0.046983 0.157103
    0.037937 0.158256
    0.029956 0.159771
    0.02311 0.161648
    0.017419 0.163862
    0.012842 0.166404
    0.009324 0.169289
    0.006854 0.172546
    0.005484 0.176154
    0.005242 0.180005
    0.006112 0.184067
    0.00809 0.188086
    0.011135 0.192004
    0.015227 0.19579
    0.020339 0.199393
    0.026403 0.202735
    0.033312 0.205736
    0.040949 0.208332
    0.049193 0.210487
    0.057935 0.212174
    0.067113 0.21339
    0.076647 0.214136
    0.086499 0.214421
    0.09661 0.214255
    0.10695 0.213649
    0.117508 0.212618
    0.12838 0.211153
    0.139859 0.209267
    0.151644 0.207051
    0.163586 0.204547
    0.175647 0.201771
    0.187807 0.198746
    0.20002 0.19549
    0.212269 0.192017
    0.224549 0.188335
    0.236794 0.18447
    0.249026 0.180416
    0.261206 0.176188
    0.273301 0.171796
    0.28524 0.16727
    0.29644 0.162842
    0.304542 0.159526
]

centerbody_coordinates = [
    0.0 0.0
    0.000586 0.005293
    0.002179 0.010047
    0.004736 0.014551
    0.008231 0.018825
    0.012632 0.022848
    0.01788 0.026585
    0.023901 0.030001
    0.030604 0.033068
    0.0379 0.035771
    0.045705 0.038107
    0.053933 0.040075
    0.06254 0.04169
    0.071451 0.042966
    0.08063 0.043916
    0.090039 0.044561
    0.09968 0.044922
    0.109361 0.044999
    0.12 0.044952
    0.135773 0.04495
    0.151899 0.04493
    0.16806 0.044913
    0.184232 0.044898
    0.200407 0.044882
    0.21658 0.044866
    0.232723 0.044847
    0.248578 0.044839
    0.262095 0.044564
    0.274184 0.043576
    0.285768 0.041795
    0.296701 0.039168
    0.306379 0.035928
]

# number of rotors
B = 5

# rotor axial location
rotorzloc = 0.12

# rotor tip radius
Rtip = 0.15572081487373543

# rotor hub radius
Rhub = 0.04495252299071941

# non-dimensional blade element radial stations
r = [
    0.050491
    0.061567
    0.072644
    0.083721
    0.094798
    0.10587
    0.11695
    0.12803
    0.13911
    0.15018
] ./ Rtip

# dimensional chord lengths
chords = [
    0.089142
    0.079785
    0.0713
    0.063979
    0.057777
    0.052541
    0.048103
    0.044316
    0.041061
    0.038243
]

# twist angles (from plane of rotation) in radians
twists = [
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
] .* pi / 180.0

# DFDC-type airfoil object
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

# all airfoils are the same
airfoils = fill(afparams, length(r)) # specify the airfoil array

# assemble rotor parameters
rotor = DuctAPE.Rotor(
    [B],
    [rotorzloc],
    r,
    [Rhub],
    [Rtip],
    chords,
    twists,
    [0.0], # currently only zero tip gaps work.
    airfoils,
    [0.0], # can flip the cl lookups on the fly if desired, say, for stator sections
)

# number of panels for the duct inlet
nduct_inlet = 30

# number of panels for the center body inlet
ncenterbody_inlet = 30

# number of panels from:
#  - rotor to duct trailing edge
#  - duct trailing edge to center body trailing edge
#  - center body trailing edge to end of wake
npanels = [30, 1, 30]

# the duct trailing edge is ahead of the centerbody trailing edge.
dte_minus_cbte = -1.0

# number of wake sheets (one more than blade elements to use)
nwake_sheets = 11

# non-dimensional wake length aft of rear-most trailing edge
wake_length = 0.8

# assemble paneling constants
paneling_constants = DuctAPE.PanelingConstants(
    nduct_inlet, ncenterbody_inlet, npanels, dte_minus_cbte, nwake_sheets, wake_length
)

# assemble ducted_rotor object
ducted_rotor = DuctAPE.DuctedRotor(
    duct_coordinates, centerbody_coordinates, rotor, paneling_constants
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
operating_point = DuctAPE.OperatingPoint(Vinf, rhoinf, muinf, asound, Omega)

# reference velocity (close to average axial velocity at rotor in this case)
Vref = 50.0

# reference radius (usually tip radius of rotor)
Rref = Rtip

# assemble reference parameters
reference_parameters = DuctAPE.ReferenceParameters([Vref], [Rref])

```

The following generates animations across the given advance ratios.

!!! warning "Plotting Streamlines"
    Currently, plotting streamlines, especially animations, takes an exceptionally long time.

```@example visualize
using Plots

# - Advance Ratio Range - #
advance_ratios = range(0.1, 2.0; step=0.01)

# - Calculate Vinfs - #
D = 2.0 * rotor.Rtip[1] # rotor diameter
n = RPM / 60.0 # rotation rate in revolutions per second
Vinfs = advance_ratios * n * D

# - Set Operating Points - #
operating_points = [deepcopy(operating_point) for i in 1:length(Vinfs)]
for (iv, v) in enumerate(Vinfs)
    operating_points[iv].Vinf[] = v
end

# - Run Multi-point Analysis - #
outs, ins, success_flags = DuctAPE.analyze(
    ducted_rotor,
    operating_points,
    reference_parameters,
    DuctAPE.set_options(
        operating_points;
        boundary_layer_options=DuctAPE.HeadsBoundaryLayerOptions(;
            model_drag=true, n_steps=1000, separation_criteria=3.0
        ),
    );
    return_inputs=true,
)

DuctAPE.generate_plots(
    DuctAPE.animatedPlots(),
    Plots, # Pass in the Plots namespace
    ins,
    outs;
    save_path="../assets/",
    static_file_type=".png",
    (;
        custom_defaults...,
        size=(600, 400),
        markersize=4,
        cp_ylim=(-3, 3), # keyword argument to set ylim for cp plots
        vtan_ylim=(0, 3), # keyword argument to set ylim for vtan plots
        bl_ylim=(0.1, 0.25), # keyword argument to set ylim for boundary layer plots
    )...,
)
```

!!! note "Custom Defaults"
    Additional arguments splatted into `generate_plots` are passed into `Plots.plot` directly as keyword arguments.  In this case, `custom_defaults` happens to be the defaults associated with the plot formatting used in these docs.

![](../assets/geometry.png)
![](../assets/surface_pressure.gif)
![](../assets/surface_velocity.gif)
![](../assets/boundary_layer.gif)
![](../assets/streamlines.gif)
