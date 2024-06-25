# Getting Started


```@contents
Pages = ["tutorial.md"]
Depth = 5
```

The following is a basic tutorial on how to set up the inputs to, and run, an analysis of a ducted fan in DuctAPE.

```@setup dfdc
include("../../assets/plots_default.jl")
gr()
```

We begin by loading the package, and optionally create a shorthand name.

```@example dfdc
using DuctAPE
const dt = DuctAPE
nothing # hide
```

## Build Inputs

The next step is to create the input object of type `Propulsor`.

```@docs; canonical=false
DuctAPE.Propulsor
```

### Body Geometry

We begin by defining a matrix of coordinates for the duct and another for the centerbody geometries, for example:

```@example dfdc
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

```@example dfdc
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

```@example dfdc
pg = plot(duct_coordinates[:,1], duct_coordinates[:,2], aspectratio=1, color=1, linewidth=2, label="Duct", xlabel="z", ylabel="r", legend=:left) # hide
plot!(pg, centerbody_coordinates[:,1], centerbody_coordinates[:,2], color=2, linewidth=2, label="Center Body") # hide
```

!!! note
    The body geometry coordinates must be input as columns of z (axial) and r (radial) coordinates, in that order.

### Rotor Geometry

The next step is to assemble an object of type `RotorStatorParameters` which contains the geometric information required to define the rotor(s) and their respective blade elements.

```@docs; canonical=false
DuctAPE.RotorStatorParameters
```

In this example, we have a single rotor defined as follows.

```@example dfdc
B = 5

rotorzloc = 0.12

Rtip = 0.15572081487373543

Rhub = 0.04495252299071941

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
]./Rtip

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
].*pi/180.0


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

rotorstator_parameters = dt.RotorStatorParameters(
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

```@example dfdc
plot!(pg, rotorzloc*ones(length(r)), r.*Rtip, seriestype=:scatter, markerstrokewidth=0, label="Blade Elements") # hide
```

!!! note "Airfoils"
    Airfoil types for DuctAPE are currently contained in the C4Blade (Cascade Compatible [CCBlade](https://flow.byu.edu/CCBlade.jl/stable/)) sub-module of DuctAPE which is exported as `c4b` and also contains the various airfoil evaluation functions used for the blade element lookups.
    The available airfoil types include all the airfoil types from CCBlade, as well as `DFDCairfoil` which is an [XROTOR](https://web.mit.edu/drela/Public/web/xrotor/)-like parametric cascade polar used in DFDC.
    In addition there are untested cascade types with similar structure to CCBlades airfoil types called `DTCascade`.
    Furthermore, there is an experimental actuator disk model implemented via the `ADM` airfoil type in C4Blade.

### Operating Point

Next we will assemble the operating point which contains information about the freestream as well as the rotor rotation rate(s).

```@docs; canonical=false
DuctAPE.OperatingPoint
```

```@example dfdc
# Freestream
Vinf = 0.0 # hover condition
rhoinf = 1.226
asound = 340.0
muinf = 1.78e-5

# Rotation Rate
RPM = 8000.0
Omega = RPM * pi / 30 # if using RPM, be sure to convert to rad/s

# utilizing the constructor function to put things in vector types
operating_point = dt.OperatingPoint(Vinf, rhoinf, muinf, asound, Omega)

nothing # hide
```

### Paneling Constants

The `PanelingConstants` object contains the constants required for DuctAPE to re-panel the provided geometry into a format compatible with the solve structure.
The `PanelingConstants` object is also used to build all of the preallocated caches inside DuctAPE, which can be done up-front if desired.
Note that there is some functionality in place for cases when the user wants to keep their own specified geometry, but this functionality should be used with caution and only by users who are certain their provided geometry is in the compatible format.  See the [Examples](@ref "Circumventing the Automated Geometry Re-paneling") for an example.

```@docs; canonical=false
DuctAPE.PanelingConstants
```

```@example dfdc
nduct_inlet = 30
ncenterbody_inlet = 30
npanels = [30, 1, 30] # the 1 is due to the fact that the duct and center body trailing edges are not quite aligned.
dte_minus_cbte = -1.0 # the duct trailing edge is ahead of the centerbody trailing edge.
nwake_sheets = 11
wake_length = 0.8

paneling_constants = dt.PanelingConstants(
    nduct_inlet, ncenterbody_inlet, npanels, dte_minus_cbte, nwake_sheets, wake_length
)
nothing # hide
```

### Reference Parameters

The reference parameters are used in the post-processing non-dimensionalizations.

```@docs; canonical=false
DuctAPE.ReferenceParameters
```

```@example dfdc
Vref = 50.0 #this turns out to be close to the average axial velocity at the rotor in our case
Rref = Rtip

reference_parameters = dt.ReferenceParameters([Vref], [Rref])
nothing # hide
```

### All Together

We are now posed to construct the `Propulsor` input type.

```@example dfdc
propulsor = dt.Propulsor(
    duct_coordinates,
    centerbody_coordinates,
    rotorstator_parameters,
    operating_point,
    paneling_constants,
    reference_parameters,
)
nothing # hide
```

## Set Options

The default options should be sufficient for just starting out and are set through the `set_options` function.

```@docs; canonical=false
DuctAPE.set_options
```

```@example dfdc
options = dt.set_options()
```

For more advanced option selection, see the examples and API reference.

## Run Analysis

With the propulsor input build, and the options selected, we are now ready to run an analysis.
This is done simply with the `analyze` function which dispatches the appropriate analysis, solve, and post-processing functions based on the selected options.

```@docs; canonical=false
DuctAPE.analyze(::DuctAPE.Propulsor, ::DuctAPE.Options)
```

```@example dfdc
outs, success_flag = dt.analyze(propulsor, options)
nothing # hide
```

## Outputs

There are many outputs contained in the named tuple output from the `analyze` function (see the [post_process() docstring](@ref DuctAPE.post_process)), but some that may be of immediate interest include:

```@example dfdc
# Total Thrust Coefficient
outs.totals.CT
```
```@example dfdc
# Total Torque Coefficient
outs.totals.CQ
```
