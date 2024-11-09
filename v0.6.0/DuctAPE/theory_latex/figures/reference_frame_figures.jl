#=
Figures for reference frame section of duct solver chapter
=#


using FLOWFoil

# Get duct wall geometry
xwall, ywall = naca4(4, 4, 20)
wallcoords = [xwall ywall]
wallangle = -30.0
walllocation = [0.0; 0.1]
FLOWFoil.position_coordinates(wallcoords, 1.5, wallangle, walllocation)
wallx = wallcoords[:, 1]
wallr = wallcoords[:, 2]

f = open("relative-frame-airfoil.dat", "w")
for (z,r) in zip(eachrow(wallx), eachrow(wallr))
    write(f,"$(z[1]) $(r[1])\n")
end
close(f)

xwall, ywall = naca4(6, 4, 20)
wallcoords = [xwall ywall]
wallangle = 55.0
walllocation = [0.0; 0.0]
FLOWFoil.position_coordinates(wallcoords, 1.5, wallangle, walllocation)
wallx = wallcoords[:, 1]
wallr = wallcoords[:, 2]

f = open("bladeelement-angles.dat", "w")
for (z,r) in zip(eachrow(wallx), eachrow(wallr))
    write(f,"$(z[1]) $(r[1])\n")
end
close(f)

xwall, ywall = naca4(10, 5, 12)
wallcoords = [xwall ywall]
wallangle = -45.0
walllocation = [0.0; 0.0]
wallscale = 1.0
FLOWFoil.position_coordinates(wallcoords, wallscale, wallangle, walllocation)
wallx = wallcoords[:, 1]
wallr = wallcoords[:, 2]

f = open("wake-screw-airfoil.dat", "w")
for (z,r) in zip(eachrow(wallx), eachrow(wallr))
    write(f,"$(z[1]) $(r[1])\n")
end
close(f)

