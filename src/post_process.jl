"""
    calculate_duct_thrust(ff_solution, Vinf; rho=1.225)

Calculate the thrust of the duct.
TODO: need to test!

**Arguments:**
- `ff_solution::FLOWFoil.InviscidSolution` : Inviscid solution object from FLOWFoil.
- `Vinf::Float` : freestream velocity

**Keyword Arguments:**
- `rho::Float` : air density

**Returns:**
- `thrust::Float` : duct thrust
"""
function calculate_duct_thrust(ff_solution, Vinf; rho=1.225)

    #unpack for convenience
    gammas = ff_solution.panelgammas
    duct_mesh = ff_solution.meshes[1]

    #calculate dynamic pressure
    q = 0.5 * rho * Vinf^2

    #initialize output
    fx = 0.0

    #get gammas specific to duct mesh (assumed to be first mesh)
    #note that FLOWFoil solves things in terms of 1/Vinf, so these surface velocities are actually Vti/Vinf already.
    duct_surface_velocities = get_mesh_gammas(gammas, duct_mesh, 1)

    # loop through panels for this mesh
    for j in 1:length(duct_mesh.panels)

        #get current panel
        panel = duct_mesh.panels[j]

        #calculate pressure on panel
        cp_panel = 1.0 - (duct_surface_velocities[j])^2

        #dimensionalize
        P = cp_panel * q

        #add panel pressure in x-direction to total sectional force
        f_x += P * panel.length * panel.normal[1]
    end

    println("fx: ", fx)
    #return total duct thrust for whole annulus: -fx*2pi
    return -fx * 2.0 * pi
end

