#=
Put functions here until ready for reorganization
=#

"""
    `solve_inviscid_problem`

uses FLOWFoil to solve inviscid problem

**Arguments:**
- duct_geometry

**Returns:**
- inviscid_system : LHS, RHS, and vortex strengths

**Calls:**
- various FLOWFoil Functions
- `mirror_geometry`
"""
function solve_inviscid_problem()
    return nothing
end

"""
    `calculate_AIC_matrices`

Calcualte the aerodynamic influence coefficient matrices of the rotor and wake panels onto the wall panels.

**Arguments:**
- all the panels

**Returns:**
- AIC matrices
"""
function calculate_AIC_matrices()
    return nothing
end

"""
    `initialize_blade_section_data`

Initialize rotor blade section aerodynamic data.

**Arguments:**
- rotors
- freestream (from inviscid panel solution)

**Returns:**
- section_circulation
- section_enthalpy_jump
- average_meridional_velocity
- rotor_velocities

**Calls:**
- CCBlade functions
- `blade_section_circulation`
- `blade_section_enthalpy`
"""
function initialize_blade_section_data()
    return nothing
end

"""
    `blade_section_circulation`

Calculates circulation of rotor blade elements

**Arguments:**
- inflow velocity
- section chord
- section lift coefficient

**Returns:**
- section_circulation
"""
function blade_section_circulation()
    return nothing
end

"""
    `blade_section_enthalpy`

Calculates enthalpy jump across rotor

**Arguments:**
- rotors
- freestream

**Returns:**
- section_enthalpy_jump
"""
function blade_section_enthalpy()
    return nothing
end

"""
    `set_rotor_panel_strengths`

Set the constant source strengths on the rotor source panels

**Arguments:**
- wall_panel_strengths
- wake_vortex_strengths

**Returns:**
- rotor_source_panel_strengths

**Calls:**
- `probe_velocity_field`
- `blade_section_sigma`
"""
function set_rotor_panel_strengths()
    return nothing
end

"""
    `probe_velocity_field`

Use freestream (inviscid panel solution) and wake panel strengths to find the velocity at a give location.

**Arguments:**
- wall panel strengths
- wake panel strengths
- point to probe

**Returns:**
- velocity at point
"""
function probe_velocity_field()
    return nothing
end

"""
    `solve_system`

Solve invsicid system (augmented with wake and rotor panel contributions) for wall vortex strengths

**Arguments:**
- inviscid_system
- AIC matrices
- wake panel strengths
- rotor panel strengths

**Returns:**
- updated wall panel strengths

**Calls:**
- linear solver
"""
function solve_system()
    return nothing
end

##########FOR NEWTON SOLVE


"""
    `update_velocities`

Description

**Arguments:**

**Returns:**

**Calls:**
"""
function update_velocities()
    return nothing
end

"""
    `update_blade_section_data`

Description

**Arguments:**

**Returns:**

**Calls:**
"""
function update_blade_section_data()
    return nothing
end

"""
    `update_wake_grid_data`

Description

**Arguments:**

**Returns:**

**Calls:**
"""
function update_wake_grid_data()
    return nothing
end

"""
    `update_wake_panel_strengths`

Description

**Arguments:**

**Returns:**

**Calls:**

"""
function update_wake_panel_strengths()
    return nothing
end

"""
    `update_source_panel_strengths`

Description

**Arguments:**

**Returns:**

**Calls:**
"""
function update_source_panel_strengths()
    return nothing
end
