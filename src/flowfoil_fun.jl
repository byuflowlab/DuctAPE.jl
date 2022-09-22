#=
Functions related to FLOWFoil
=#

"""

Obtain initial wall panel vortex strengths using symmetric geometry; also return A matrix and boundary condition vector for axisymmetric geometry.
#TODO requires some modification to FLOWFoil code to allow for intermediate outputs.

**Arguments:**
- `systempanels::PanelSystem`
- `freestream::Freestream`

**Returns:**
- `wall_panel_strengths::Vector{Float64}`
- `Avinf::Matrix{Float64,2}`
- `bvinf::Array{Float64}`
"""
function first_guess_flowfield(systempanels, rotors, freestream)

    # get LHS and RHS from system
    Avinf, bvinf = get_axisymmetric_coefficients(systempanels)

    # get first guess for velocities at rotors
    rotor_velocities = first_guess_rotor_velocities(systempanels, rotors, freestream)

    # return everything.
    return rotor_velocities, Avinf, bvinf
end

"""
"""
function get_axisymmetric_coefficients(systempanels)

    # convert panels to FLOWFoil geometry inputs
    duct_meshes = duct2foil(systempanels)

    # define Problem
    problem = FLOWFoil.Problem(meshes; viscous=false)

    #return LHS and RHS matrix and array
    return FLOWFoil.get_system(problem)
end

"""
"""
function duct2foil(systempanels)
    return duct_meshes
end

"""
"""
function first_guess_rotor_velocities(systempanels, rotors, freestream)

    # set up mirrored system
    mirror_meshes = duct2foilmirror(systempanels)

    #define Problem
    problem = FLOWFoil.Problem(mirror_meshes; viscous=false)

    #solve problem
    solution = FLOWFoil.solve(problem)

    #get flow velocities at rotor stations
    velocities = gamma2vel(solution, rotors, freestream)

    # initialize rotor_velocities object
    rotor_velocities = initialize_rotor_velocities(velocities, rotors, freestream)

    return rotor_velocities
end

"""
"""
function duct2foilmirror(systempanels)
    return mirror_meshes
end

"""
"""
function gamma2vel(panel_solution, rotors, freestream)
    return velocities
end

"""
TODO: move to rotors.jl
"""
function initialize_rotor_velocities(velocities, rotors, freestream)
    return RotorVelocities()
end
