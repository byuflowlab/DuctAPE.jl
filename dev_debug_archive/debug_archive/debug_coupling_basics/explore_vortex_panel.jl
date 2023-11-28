using FLOWFoil
const ff = FLOWFoil
using DuctAPE
const dt = DuctAPE

xcoords = [1.0; 0.0; -1.0; 0.0; 1.0]
rcoords = [1.5; 0.5; 1.5; 2.5; 1.5]
coords = [xcoords rcoords]
npan = length(xcoords) - 1

method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false])
panels = ff.generate_panels(method, coords)
mesh = dt.generate_one_way_mesh(panels, panels)

## -- Make sure to understand the panel angles -- ##
beta = panels.panel_angle*180.0/pi
println("β (in degrees):")
display(beta)

## -- Make sure to understand what sinβ and cosβ are doing -- ##
println("cosβ:")
display(cosd.(beta))
println("sinβ:")
display(sind.(beta))

## -- understand how a panel induces surface velocity on another panel -- ##
println("ξ:")
display(mesh.x)
println("ρ:")
display(mesh.r)
println("k^2:")
display(mesh.m)

vx = zeros(npan, npan)
vr = zeros(npan, npan)
for i in 1:npan #affected panels
    for j in 1:npan #influencing panels
        vx[i, j] = dt.get_vx_ring_vortex_off_body(
            mesh.x[i, j],
            mesh.r[i, j],
            panels.panel_center[j, 2],
            mesh.m[i, j],
        )
        vr[i, j] = dt.get_vr_ring_vortex_off_body(
            mesh.x[i, j],
            mesh.r[i, j],
            panels.panel_center[j, 2],
            mesh.m[i, j],
        )
    end
end
println("v_x:")
display(vx)
println("v_r:")
display(vr)
