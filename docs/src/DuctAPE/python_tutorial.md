# Calling DuctAPE from Python

In this example we repeat the [quick start](tutorial.md), but with Python.  We won't repeat all the usage details described in the quick start.  We are using the [JuliaCall](https://juliapy.github.io/PythonCall.jl/stable/) packages to enable this.

## Setup

```python
# optional: point python to a specific julia installation
# import os
# os.environ["JULIA_BINDIR"] = path_to_julia

from juliacall import Main as jl

# If needed, install DuctAPE into the python julia environment
jl.seval("using Pkg")
jl.Pkg.add("DuctAPE")

# Load DuctAPE
jl.seval("using DuctAPE")

import numpy as np
```

## Example

```python
# - Same Quick Start as Before, but with some minor changes for Python syntax - #
duct_coordinates = np.array(
    [
        [0.304466, 0.158439],
        [0.294972, 0.158441],
        [0.28113, 0.158423],
        [0.266505, 0.158365],
        [0.251898, 0.158254],
        [0.237332, 0.158088],
        [0.222751, 0.157864],
        [0.208123, 0.157586],
        [0.193399, 0.157258],
        [0.178507, 0.156897],
        [0.16349, 0.156523],
        [0.148679, 0.156177],
        [0.134222, 0.155902],
        [0.12, 0.155721],
        [0.106044, 0.155585],
        [0.092531, 0.155498],
        [0.079836, 0.155546],
        [0.067995, 0.155792],
        [0.057025, 0.156294],
        [0.046983, 0.157103],
        [0.037937, 0.158256],
        [0.029956, 0.159771],
        [0.02311, 0.161648],
        [0.017419, 0.163862],
        [0.012842, 0.166404],
        [0.009324, 0.169289],
        [0.006854, 0.172546],
        [0.005484, 0.176154],
        [0.005242, 0.180005],
        [0.006112, 0.184067],
        [0.00809, 0.188086],
        [0.011135, 0.192004],
        [0.015227, 0.19579],
        [0.020339, 0.199393],
        [0.026403, 0.202735],
        [0.033312, 0.205736],
        [0.040949, 0.208332],
        [0.049193, 0.210487],
        [0.057935, 0.212174],
        [0.067113, 0.21339],
        [0.076647, 0.214136],
        [0.086499, 0.214421],
        [0.09661, 0.214255],
        [0.10695, 0.213649],
        [0.117508, 0.212618],
        [0.12838, 0.211153],
        [0.139859, 0.209267],
        [0.151644, 0.207051],
        [0.163586, 0.204547],
        [0.175647, 0.201771],
        [0.187807, 0.198746],
        [0.20002, 0.19549],
        [0.212269, 0.192017],
        [0.224549, 0.188335],
        [0.236794, 0.18447],
        [0.249026, 0.180416],
        [0.261206, 0.176188],
        [0.273301, 0.171796],
        [0.28524, 0.16727],
        [0.29644, 0.162842],
        [0.304542, 0.159526],
    ]
)

center_body_coordinates = np.array(
    [
        [0.0, 0.0],
        [0.000586, 0.005293],
        [0.002179, 0.010047],
        [0.004736, 0.014551],
        [0.008231, 0.018825],
        [0.012632, 0.022848],
        [0.01788, 0.026585],
        [0.023901, 0.030001],
        [0.030604, 0.033068],
        [0.0379, 0.035771],
        [0.045705, 0.038107],
        [0.053933, 0.040075],
        [0.06254, 0.04169],
        [0.071451, 0.042966],
        [0.08063, 0.043916],
        [0.090039, 0.044561],
        [0.09968, 0.044922],
        [0.109361, 0.044999],
        [0.12, 0.044952],
        [0.135773, 0.04495],
        [0.151899, 0.04493],
        [0.16806, 0.044913],
        [0.184232, 0.044898],
        [0.200407, 0.044882],
        [0.21658, 0.044866],
        [0.232723, 0.044847],
        [0.248578, 0.044839],
        [0.262095, 0.044564],
        [0.274184, 0.043576],
        [0.285768, 0.041795],
        [0.296701, 0.039168],
        [0.306379, 0.035928],
    ]
)

# number of rotors
B = 5

# rotor axial location
rotor_axial_position = 0.12

# rotor tip radius
Rtip = 0.15572081487373543

# rotor hub radius
Rhub = 0.04495252299071941

# non-dimensional blade element radial stations
r = (
    np.array(
        [
            0.050491,
            0.061567,
            0.072644,
            0.083721,
            0.094798,
            0.10587,
            0.11695,
            0.12803,
            0.13911,
            0.15018,
        ]
    )
    / Rtip
)

# dimensional chord lengths
chords = np.array(
    [
        0.089142,
        0.079785,
        0.0713,
        0.063979,
        0.057777,
        0.052541,
        0.048103,
        0.044316,
        0.041061,
        0.038243,
    ]
)

# twist angles (from plane of rotation) in radians
twists = (
    np.array(
        [
            69.012,
            59.142,
            51.825,
            46.272,
            41.952,
            38.509,
            35.699,
            33.354,
            31.349,
            29.596,
        ]
    )
    * np.pi
    / 180.0
)

# DFDC-type airfoil object
afparams = jl.DuctAPE.c4b.DFDCairfoil(
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
# NOTE: airfoils are inputs as a vector of vectors rather than a matrix
airfoils = [[afparams] * len(r)]

# Assemble Rotor object
rotor = jl.DuctAPE.Rotor(
    B,
    rotor_axial_position,
    r,
    Rhub,
    Rtip,
    chords,
    twists,
    [0.0],  # tip gap vector
    airfoils,  # vector-of-vectors of DFDCairfoil objects
    [0.0],  # is_stator flag vector
)

# number of panels for the duct inlet
num_duct_inlet_panels = 30

# number of panels for the center body inlet
num_center_body_inlet_panels = 30

# number of panels:
#  - rotor → duct TE
#  - duct TE → center-body TE
#  - center-body TE → end of wake
num_panels = [30, 1, 30]

# duct TE ahead of center-body TE
dte_minus_cbte = -1.0

# number of wake sheets (one more than blade elements)
num_wake_sheets = 11

# non-dimensional wake length aft of rear-most TE
wake_length = 0.8

# assemble paneling constants
paneling_constants = jl.DuctAPE.PanelingConstants(
    num_duct_inlet_panels,
    num_center_body_inlet_panels,
    num_panels,
    dte_minus_cbte,
    num_wake_sheets,
    wake_length,
)

# assemble ducted rotor object
ducted_rotor = jl.DuctAPE.DuctedRotor(
    duct_coordinates,
    center_body_coordinates,
    rotor,
    paneling_constants,
)

# Freestream
Vinf = 30.0
rhoinf = 1.226
asound = 340.0
muinf = 1.78e-5

# Rotation rate
RPM = 8000.0
Omega = RPM * np.pi / 30.0  # rad/s

# assemble operating point
operating_point = jl.DuctAPE.OperatingPoint(
    Vinf,
    Omega,
    rhoinf,
    muinf,
    asound,
)

# reference velocity (close to average axial velocity at rotor)
Vref = 50.0

# reference radius (usually rotor tip radius)
Rref = Rtip

# assemble reference parameters
reference_parameters = jl.DuctAPE.ReferenceParameters(
    [Vref],
    [Rref],
)

# set options
options = jl.DuctAPE.set_options()

# Run DuctAPE
outs, success_flag = jl.DuctAPE.analyze(
    ducted_rotor,
    operating_point,
    reference_parameters,
    options,
)

# Print some outputs.
print("CT: ", outs.totals.CT)
print("CQ: ", outs.totals.CQ)
```

## Automatic Derivatives

We now continue the example demonstrating how to get derivatives for use in Python (but where the derivative computation occurs in Julia via algorithmic differentiation).  For this functionality we need to load the [PythonCall](https://juliapy.github.io/PythonCall.jl/stable/) package (which enables us to call back into Python from Julia), and we need the [ImplicitAD](https://github.com/byuflowlab/ImplicitAD.jl) packages which provides a convenience function for the derivative computation.  Alternatively, for more advanced users you can just use the Julia differentiation packages directly (ForwardDiff, ReverseDiff, etc.).  Note that for Julia AD to work, all the function calls will need to be Julia function calls.  Even though all the setup is happening here in Python, we are only setting up inputs.  All functions are calls to Julia (jl.somefunction())

The below example demonstrates a forward mode Jacobian and a forward mode Jacobian-vector product.  Other options (reverse Jacobian, reverse vector-Jacobian product) are discussed in [ImplicitAD](https://github.com/byuflowlab/ImplicitAD.jl) docs.

```python
# load PythonCall BEFORE ImplicitAD
jl.seval("using PythonCall")

# Install ImplicitAD as needed
jl.seval("using Pkg")
jl.Pkg.add("ImplicitAD")

# Load ImplicitAD.derivativesetup
jl.seval("using ImplicitAD: derivativesetup")

nc = len(chords)

# ImplicitAD expects a funciton we want to differentiate in the form f = func(x, p)
# where f is output vector, x is input vector, and p are parameters we do not differentiate w.r.t.
def dtwrap(x, p):
    chords = x[:nc]
    twists = x[nc:]

    # Assemble Rotor object
    rotor = jl.DuctAPE.Rotor(
        B,
        rotor_axial_position,
        r,
        Rhub,
        Rtip,
        chords,
        twists,
        [0.0],  # tip gap vector
        airfoils,  # vector-of-vectors of DFDCairfoil objects
        [0.0],  # is_stator flag vector
    )

    # assemble ducted rotor object
    ducted_rotor = jl.DuctAPE.DuctedRotor(
        duct_coordinates,
        center_body_coordinates,
        rotor,
        paneling_constants,
    )

    # Run DuctAPE
    outs, success_flag = jl.DuctAPE.analyze(
        ducted_rotor,
        operating_point,
        reference_parameters,
        options,
    )

    return [outs.totals.CT, outs.totals.CQ]


x = np.concatenate([chords, twists])
p = ()
jacobian = jl.derivativesetup(
    dtwrap, x, p, "fjacobian"
)  # a forward-mode Jacobian is one option

# preallocate Jacobian then evaluate
J = np.zeros((2, len(x)))
jacobian(J, x)
print(J)
# can now change x, and evaluate jacobian(J, x) repeatedly at other points


# demonstrate a Jacobian-vector product
jvp = jl.derivativesetup(dtwrap, x, p, "jvp")
xdot = np.ones(len(x))
fdot = np.zeros(2)
jvp(fdot, x, xdot)
print(fdot)
# can continue to call jvp for different x, xdot pairs
```
