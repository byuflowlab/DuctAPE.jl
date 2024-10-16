# Getting Started


```@contents
Pages = ["tutorial.md"]
Depth = 5
```

The following is a basic tutorial on how to set up and run an analysis of a ducted fan in DuctAPE.

```@setup tutorial
include("../assets/plots_default.jl")
gr()
```

We begin by loading the package:

```@example tutorial
using DuctAPE
nothing # hide
```

## Assemble Inputs

The next step is to create the input object of type `DuctedRotor`.

```@docs; canonical=false
DuctAPE.DuctedRotor
```

### Body Geometry

We begin by defining a matrix of coordinates for the duct and another for the centerbody geometries.
For example:

```@example tutorial
duct_coordinates = [
    0.304466  0.158439
    0.294972  0.158441
    0.28113   0.158423
    0.266505  0.158365
    0.251898  0.158254
    0.237332  0.158088
    0.222751  0.157864
    0.208123  0.157586
    0.193399  0.157258
    0.178507  0.156897
    0.16349   0.156523
    0.148679  0.156177
    0.134222  0.155902
    0.12      0.155721
    0.106044  0.155585
    0.092531  0.155498
    0.079836  0.155546
    0.067995  0.155792
    0.057025  0.156294
    0.046983  0.157103
    0.037937  0.158256
    0.029956  0.159771
    0.02311   0.161648
    0.017419  0.163862
    0.012842  0.166404
    0.009324  0.169289
    0.006854  0.172546
    0.005484  0.176154
    0.005242  0.180005
    0.006112  0.184067
    0.00809   0.188086
    0.011135  0.192004
    0.015227  0.19579
    0.020339  0.199393
    0.026403  0.202735
    0.033312  0.205736
    0.040949  0.208332
    0.049193  0.210487
    0.057935  0.212174
    0.067113  0.21339
    0.076647  0.214136
    0.086499  0.214421
    0.09661   0.214255
    0.10695   0.213649
    0.117508  0.212618
    0.12838   0.211153
    0.139859  0.209267
    0.151644  0.207051
    0.163586  0.204547
    0.175647  0.201771
    0.187807  0.198746
    0.20002   0.19549
    0.212269  0.192017
    0.224549  0.188335
    0.236794  0.18447
    0.249026  0.180416
    0.261206  0.176188
    0.273301  0.171796
    0.28524   0.16727
    0.29644   0.162842
    0.304542  0.159526
]
nothing # hide
```

```@example tutorial
centerbody_coordinates = [
    0.0       0.0
    0.000586  0.005293
    0.002179  0.010047
    0.004736  0.014551
    0.008231  0.018825
    0.012632  0.022848
    0.01788   0.026585
    0.023901  0.030001
    0.030604  0.033068
    0.0379    0.035771
    0.045705  0.038107
    0.053933  0.040075
    0.06254   0.04169
    0.071451  0.042966
    0.08063   0.043916
    0.090039  0.044561
    0.09968   0.044922
    0.109361  0.044999
    0.12      0.044952
    0.135773  0.04495
    0.151899  0.04493
    0.16806   0.044913
    0.184232  0.044898
    0.200407  0.044882
    0.21658   0.044866
    0.232723  0.044847
    0.248578  0.044839
    0.262095  0.044564
    0.274184  0.043576
    0.285768  0.041795
    0.296701  0.039168
    0.306379  0.035928
]
nothing # hide
```

```@example tutorial
pg = plot( # hide
    duct_coordinates[:, 1], # hide
    duct_coordinates[:, 2]; # hide
    aspectratio=1, # hide
    color=1, # hide
    linewidth=2, # hide
    label="Duct", # hide
    xlabel="z", # hide
    ylabel="r", # hide
    legend=:left, # hide
) # hide
plot!( # hide
    pg, # hide
    centerbody_coordinates[:, 1], # hide
    centerbody_coordinates[:, 2]; # hide
    color=2, # hide
    linewidth=2, # hide
    label="Center Body", # hide
) # hide
```

!!! note
    The body geometry coordinates must be input as columns of z (axial) and r (radial) coordinates, in that order.

### Rotor Geometry

The next step is to assemble an object of type `Rotor` which contains the geometric information required to define the rotor(s) and their respective blade elements.

```@docs; canonical=false
DuctAPE.Rotor
```

In this example, we have a single rotor defined as follows.

```@example tutorial
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
nothing # hide
```

```@example tutorial
plot!( # hide
    pg, # hide
    rotorzloc * ones(length(r)), # hide
    r .* Rtip; # hide
    seriestype=:scatter, # hide
    markersize=3, # hide
    markerstrokewidth=0, # hide
    label="Blade Elements", # hide
) # hide
```

!!! note "Airfoils"
    Airfoil types for DuctAPE are currently contained in the C4Blade (Cascade Compatible [CCBlade](https://flow.byu.edu/CCBlade.jl/stable/)) sub-module of DuctAPE which is exported as `c4b` and also contains the various airfoil evaluation functions used for the blade element lookups.
    The available airfoil types include all the airfoil types from CCBlade, as well as `DFDCairfoil` which is an [XROTOR](https://web.mit.edu/drela/Public/web/xrotor/)-like parametric cascade polar used in DFDC.
    In addition there are untested cascade types with similar structure to CCBlades airfoil types called `DTCascade`.
    Furthermore, there is an experimental actuator disk model implemented via the `ADM` airfoil type in C4Blade.

### Paneling Constants

The `PanelingConstants` object contains the constants required for DuctAPE to re-panel the provided geometry into a format compatible with the solve structure.
Specifically, the DuctAPE solver makes some assumptions on the relative positioning of the body surfaces relative to the wakes and each other; and this is most easily guarenteed by a re-paneling of the provided body surface geometry.
The `PanelingConstants` object is also used to build all of the preallocated caches inside DuctAPE, which can be done up-front if desired.
Note that there is some functionality in place for cases when the user wants to keep their own specified geometry, but this functionality should be used with caution and only by users who are certain their provided geometry is in the compatible format.  See the [Examples](@ref "Circumventing the Automated Geometry Re-paneling") for an example.

```@docs; canonical=false
DuctAPE.PanelingConstants
```

```@example tutorial
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
nothing # hide
```
### Assembling the DuctedRotor

We are now posed to construct the `DuctedRotor` input type.

```@example tutorial
# assemble ducted_rotor object
ducted_rotor = DuctAPE.DuctedRotor(
    duct_coordinates,
    centerbody_coordinates,
    rotor,
    paneling_constants,
)
nothing # hide
```

### Operating Point

Next we will assemble the operating point which contains information about the freestream as well as the rotor rotation rate(s).

```@docs; canonical=false
DuctAPE.OperatingPoint
```

```@example tutorial
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
nothing # hide
```

### Reference Parameters

The reference parameters are used in the post-processing non-dimensionalizations.

```@docs; canonical=false
DuctAPE.ReferenceParameters
```

```@example tutorial
# reference velocity (close to average axial velocity at rotor in this case)
Vref = 50.0

# reference radius (usually tip radius of rotor)
Rref = Rtip

# assemble reference parameters
reference_parameters = DuctAPE.ReferenceParameters([Vref], [Rref])
nothing # hide
```


## Set Options

The default options should be sufficient for just starting out and are set through the `set_options` function.

```@docs; canonical=false
DuctAPE.set_options
```

```@example tutorial
options = DuctAPE.set_options()
```

For more advanced option selection, see the examples and API reference.

## Run a Single Analysis

With the ducted_rotor input build, and the options selected, we are now ready to run an analysis.
This is done simply with the `analyze` function which dispatches the appropriate analysis, solve, and post-processing functions based on the selected options.

```@docs; canonical=false
DuctAPE.analyze(::DuctAPE.DuctedRotor, ::DuctAPE.OperatingPoint, ::DuctAPE.ReferenceParameters, ::DuctAPE.Options)
```

```@example tutorial
outs, success_flag = DuctAPE.analyze(
    ducted_rotor, operating_point, reference_parameters, options
)
nothing # hide
```

### Single Run Outputs

There are many outputs contained in the named tuple output from the `analyze` function (see the [post_process docstring](@ref DuctAPE.post_process)), but some that may be of immediate interest include:

```@example tutorial
# Total Thrust Coefficient
outs.totals.CT
```
```@example tutorial
# Total Torque Coefficient
outs.totals.CQ
```

## Run a Multi-Point Analysis

In the case that one wants to run the same geometry at several different operating points, for example: for a range of advance ratios, there is another dispatch of the `analyze` function that accepts an input, `multipoint`, that is a vector of operating points.

```@docs; canonical=false
DuctAPE.analyze(ducted_rotor::DuctedRotor,operating_point::AbstractVector{TO},reference_parameters::ReferenceParameters,options::Options) where TO<:OperatingPoint
```

Running a multi-point analysis on the example geometry given there, it might look something like this:

```@example tutorial
# - Advance Ratio Range - #
Js = range(0.0, 2.0; step=0.01)

# - Calculate Vinfs - #
D = 2.0 * rotor.Rtip[1] # rotor diameter
n = RPM / 60.0 # rotation rate in revolutions per second
Vinfs = Js * n * D

# - Set Operating Points - #
operating_points = [deepcopy(operating_point) for i in 1:length(Vinfs)]
for (iv, v) in enumerate(Vinfs)
    operating_points[iv].Vinf[] = v
end

# - Run Multi-point Analysis - #
outs_vec, success_flags = DuctAPE.analyze(
    ducted_rotor,
    operating_points,
    reference_parameters,
    DuctAPE.set_options(operating_points),
)
nothing #hide
```

There are a few things to note here.
1. We want to make sure that the operating point objects we put into the input vector are unique instances.
2. We need to use the dispatch of `set_options` that accepts the operating point vector to set up the right number of things in the background (like convergence flags for each operating point).
3. The outputs of the analysis are vectors of the same outputs for a single analysis.

### Multi-point Outputs

For multi-point analysis outputs, which are given as a vector of output objects, we might access and plot things as follows.
We also take the opportunity to present some verification against DFDC, showing that DuctAPE matches remarkably well (within 0.5%) of DFDC.
We therefore first provide data from DFDC analyses of the above example geometry at various advance ratios.

```@example tutorial
# Verification Data From DFDC

dfdc_jept = [
    0.0 0.0 0.64763 0.96692
    0.1 0.1366 0.64716 0.88394
    0.2 0.2506 0.6448 0.80785
    0.3 0.3457 0.64044 0.73801
    0.4 0.4251 0.63401 0.67382
    0.5 0.4915 0.62534 0.61468
    0.6 0.547 0.61428 0.56001
    0.7 0.5935 0.6006 0.50925
    0.8 0.6326 0.58411 0.46187
    0.9 0.6654 0.56452 0.41738
    1.0 0.693 0.54158 0.37531
    1.1 0.716 0.51499 0.33522
    1.2 0.7349 0.48446 0.2967
    1.3 0.7499 0.44966 0.25937
    1.4 0.7606 0.41031 0.2229
    1.5 0.7661 0.36604 0.18694
    1.6 0.7643 0.31654 0.15121
    1.7 0.7506 0.26153 0.11547
    1.8 0.7126 0.20061 0.07941
    1.9 0.61 0.13355 0.04287
    2.0 0.1861 0.05993 0.00558
]

dfdc_J = dfdc_jept[:,1]
dfdc_eta = dfdc_jept[:,2]
dfdc_cp = dfdc_jept[:,3]
dfdc_ct = dfdc_jept[:,4]
nothing #hide
```

We can then access the various multi-point analysis outputs however is convenient, we choose a broadcasting approach here:

```@example tutorial
# - extract efficiency, power, and thrust coefficients - #

# efficiency
eta = (p->p.totals.total_efficiency[1]).(outs_vec)

# power
cp = (p->p.totals.CP[1]).(outs_vec)

# thrust
ct = (p->p.totals.CT[1]).(outs_vec)
nothing #hide
```

And then we can plot the data to compare DFDC and DuctAPE.

```julia
using Plots

# set up efficiency plot
pe = plot(; xlabel="Advance Ratio", ylabel="Efficiency")

# plot DFDC data
plot!(
    pe,
    dfdc_J,
    dfdc_eta;
    seriestype=:scatter,
    markersize=5,
    markercolor=plotsgray, # hide
    markerstrokecolor=plotsgray, # hide
    label="DFDC",
)

# Plot DuctAPE outputs
plot!(
    pe,
    Js,
    eta;
    linewidth=2,
    color=primary, # hide
    label="DuctAPE",
)

# setup cp/ct plot
ppt = plot(; xlabel="Advance Ratio")

# plot DFDC data
plot!(
    ppt,
    dfdc_J,
    dfdc_cp;
    seriestype=:scatter,
    markersize=5,
    markercolor=plotsgray, # hide
    markerstrokecolor=primary, # hide
    markerstrokewidth=2,
    label="DFDC Cp",
)

plot!(
    ppt,
    dfdc_J,
    dfdc_ct;
    seriestype=:scatter,
    markersize=5,
    markercolor=plotsgray, # hide
    markerstrokecolor=secondary, # hide
    markerstrokewidth=2,
    label="DFDC Ct",
)

# plot DuctAPE outputs
plot!(
    ppt,
    Js,
    cp;
    linewidth=1.5,
    color=primary, # hide
    label="DuctAPE Cp",
)

plot!(
    ppt,
    Js,
    ct;
    linewidth=1.5,
    color=secondary, # hide
    label="DuctAPE Ct",
)

plot(
    pe,
    ppt;
    size=(700, 350),
    layout=(1, 2),
    margin=2mm, # hide
)
```
