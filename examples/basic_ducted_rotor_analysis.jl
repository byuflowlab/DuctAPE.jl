using FLOWFoil
using DuctTAPE
using CCBlade
using FLOWMath
using PyPlot

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

# rhub2 = 1e-8 * ones(length(rhub))
# xduct = [reverse(xhub); xhub[2:end]]
# rduct = [0.5 * ones(length(rhub)); rhub[2:end] .+ 0.5]

# generate mesh objects
duct = FLOWFoil.generate_axisym_mesh(xduct, rduct; bodyofrevolution=false)
hub = FLOWFoil.generate_axisym_mesh(xhub, rhub; bodyofrevolution=true)
meshes = [duct; hub]

## -- Set up FLOWFoil Problem
problem = FLOWFoil.Problem(meshes; axisymmetric=true, viscous=false)

## -- Solve FLOWFoil System
ff_solution = FLOWFoil.solve(problem)

##### ----- DEFINE ROTOR

#set up parameters
NR = 5
xlocation = 0.25
numblades = 5
nref = 25
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

##### ----- DEFINE FREESTREAM
vinf = 5.0
rho = 1.225 #kg/m3
mu = 1.81e-5
asound = 343.0

freestream = DuctTAPE.Freestream(vinf, rho, mu, asound)

##### ----- COUPLE DUCT -> CCBLADE
ccbrotor, sections, op, rdist = DuctTAPE.ff2ccb(ff_solution, rotor, freestream)

# Call CCBlade
ccb_out = CCBlade.solve.(Ref(ccbrotor), sections, op[1])

ccbrotor, sections, op, rdist = DuctTAPE.ff2ccb(ff_solution, rotor, freestream; debug=true)
ccb_outnoduct = CCBlade.solve.(Ref(ccbrotor), sections, op)

figure(1)
clf()
plot(rdist, ccb_out.Np; label="with duct")
plot(rdist, ccb_outnoduct.Np; label="no duct")
xlabel("Radial Station (m)")
# ylabel("Inflow Angle")
ylabel(L"N_P")
legend()
savefig("ccbladecomp.png"; bbox_inches="tight")

figure(2)
clf()
plot(xduct, rduct)
plot(xhub, rhub)
plot([rotor.xlocation; rotor.xlocation], [rdist[1]; rdist[end]]; linewidth=2)
savefig("ductrotorgeom.png"; bbox_inches="tight")
