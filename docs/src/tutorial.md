# Quick Start Guide

In the current state, DuctTAPE.jl can only be used in conjunction with FLOWFoil and CCBlade, of which only CCBlade is currently publicly available.

```
using CCBlade
using FLOWFoil
using DuctTAPE
```

We first need to call FLOWFoil to obtain an inviscid panel solution from which we can calculate the velocity field.
```
##### ----- CALL FLOWFOIL

### --- Create Duct Geometry and Run FLOWFoil

## -- Geometry

#use a symmetric naca airfoil for the hub
_, xhub, _, rhub = FLOWFoil.naca4(0.0, 0.0, 10.0; split=true)

#use a flipped naca airfoil for the duct
xduct, rduct = FLOWFoil.naca4(6.0, 4.0, 20.0)
coords = [xduct rduct]
scale = 1.0
angle = 8.0
location = [0.0; 0.5]
xduct, rduct = FLOWFoil.position_coordinates(coords, scale, angle, location; flipped=true)

# generate mesh objects
duct = FLOWFoil.generate_axisym_mesh(xduct, rduct; bodyofrevolution=false)
hub = FLOWFoil.generate_axisym_mesh(xhub, rhub; bodyofrevolution=true)
meshes = [duct; hub]

## -- Set up FLOWFoil Problem
problem = FLOWFoil.Problem(meshes; axisymmetric=true, viscous=false)

## -- Solve FLOWFoil System
ff_solution = FLOWFoil.solve(problem)
```

Next, we need to define the rotor object we would like to use.

```
##### ----- DEFINE ROTOR

#set up parameters
NR = 5 #number of radial stations
xlocation = 0.25
numblades = 5
nref = 25 #number of radial stations to use in refinement
radialstations = range(0.0, 1.0; length=NR)
chords = 0.1 * ones(NR)
twists = range(40.0, 10.0; length=NR)
datapath = "./data/airfoils/"
airfoils = fill(CCBlade.AlphaAF(datapath * "ccb_naca4412.dat"), NR)
rpm = 5000.0

#define rotor geometry object
rotor = DuctTAPE.RotorGeometry(
    xlocation, numblades, nref, radialstations, chords, twists, airfoils, rpm
)
```

The last piece we need is the freestream information we'd like to use.

```
##### ----- DEFINE FREESTREAM
vinf = 5.0
rho = 1.225 #kg/m3
mu = 1.81e-5
asound = 343.0

freestream = DuctTAPE.Freestream(vinf, rho, mu, asound)
```

We are now ready to get the inputs for, and run, the CCBlade solution function.

```
##### ----- COUPLE DUCT -> CCBLADE
ccbrotor, sections, op, rdist = DuctTAPE.ff2ccb(ff_solution, rotor, freestream)

# Call CCBlade
ccb_out = CCBlade.solve.(Ref(ccbrotor), sections, op)
```


