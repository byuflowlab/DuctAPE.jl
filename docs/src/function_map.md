# Function Map

Attempt to list functions in order of called along with ins/outs and links to functions called within.

KEY:

### `Main Functions`

#### _`sub fuctions`_



---

### `Template1`

Description

**Arguments:**

**Returns:**

**Calls:**

---

#### _`Template2`_

Description

**Arguments:**

**Returns:**

**Calls:**

---





### `generate_geometry`

Generates duct geometry object

**Arguments:**
- wall coordinates
- hub coordinates

**Returns:**
- duct geometry object

**Calls:**
- nothing

---

### `solve_inviscid_problem`

uses FLOWFoil to solve inviscid problem

**Arguments:**
- duct_geometry

**Returns:**
- inviscid_system : LHS, RHS, and vortex strengths

**Calls:**
- various FLOWFoil Functions

---

### `generate_rotors`

Generates rotor objects

**Arguments:**
- rotor geometry arrays
- rotation speeds

**Returns:**
- array of rotor objects

**Calls:**
- nothing

---

### `generate_wake_grid`

Generates elliptic wake grid

**Arguments:**
- duct_geometry
- rotors
- grid_options (number of radial stations, wake length, wake expansion factor)

**Returns:**
- elliptic grid points

**Calls:**
- [grid_initalization](#gridinitalization) function
- [grid_relaxation](#gridrelaxation) function

---

#### _`grid_initalization`_

Initializes wake grid using conservation of mass, attempting to set aspect ratio to equal while lining up with rotors and trailing edges.  Also sets initial trailing wake to minimum of desired length applying expansion factor.

**Arguments:**
- duct_geometry
- rotors
- grid_options

**Returns:**
- initial_grid_points

**Calls:**
- nothing

---

#### _`grid_relaxation`_

relaxes grid using SLOR

**Arguments:**
- initial_grid_points
- convergence_criteria

**Returns:**
- elliptic_grid_points

**Calls:**
- nothing

---

### `generate_wake_panels`

Define panel edges, centers, and unit normals for wake vortex sheets.

**Arguments:**
- elliptic_grid_points

**Returns:**
- wake panels

**Calls:**
- nothing

---


### `generate_rotor_panels`

Define panel edges, centers, and unit normals for rotor source sheets.

**Arguments:**
- rotors
- elliptic_grid_points

**Returns:**
- rotor panels

**Calls:**
- reinterpolate_rotor!

---

#### _`reinterpolate_rotor!`_

Reinterpolates rotor data using wake grid positions.  Original data is splined, then sampled at the new locations to reinterpolate.

**Arguments:**
- rotors
- elliptic_grid_points

**Returns:**
- updates rotors in place

**Calls:**
- FLOWMath functions

---

### `calculate_AIC_matrices`

Calcualte the aerodynamic influence coefficient matrices of the rotor and wake panels onto the wall panels.

**Arguments:**
- all the panels

**Returns:**
- AIC matrices

**Calls:**
- nothing

---

### `initialize_blade_section_data`

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
- [blade_section_circulation](#rotorcirculation)
- [blade_section_enthalpy](#rotorenthalpy)

---

#### _`blade_section_circulation`_

Calculates circulation of rotor blade elements

**Arguments:**
- inflow velocity
- section chord
- section lift coefficient

**Returns:**
- section_circulation

**Calls:**
- nothing

---

#### _`blade_section_enthalpy`_

Calculates enthalpy jump across rotor

**Arguments:**
- rotors
- freestream

**Returns:**
- section_enthalpy_jump

**Calls:**
- nothing

---

### `initialize_wake_data`

Initializes wake grid aerodynamic properties (flow data)

**Arguments:**
- section_circulation
- section_enthalpy
- average_meridional_velocity
- grid points

**Returns:**
- grid circulation
- grid enthalpy

**Calls:**

---

### `set_wake_panel_strengths`

Sets the wake panel vorticity strengths.

**Arguments:**
- grid circulation
- grid enthalpy
- grid points
- average_meridional_velocity

**Returns:**
- wake_vortex_strengths

**Calls:**
- [calc_gamma_i](#calcgammai)

---

#### _`calc_gamma_i`_

Calculates panel vortex strength

**Arguments:**
- grid circulation
- grid enthalpy
- radial location
- average_meridional_velocity

**Returns:**
- panel_vortex_strength

**Calls:**
- nothing

---

### `set_rotor_panel_strengths`

Set the constant source strengths on the rotor source panels

**Arguments:**
- wall_panel_strengths
- wake_vortex_strengths

**Returns:**
- rotor_source_panel_strengths

**Calls:**
- [probe_velocity_field](#probevelocityfield)
- [blade_section_sigma](#bladesectionsigma)

---

#### _`probe_velocity_field`_

Use freestream (inviscid panel solution) and wake panel strengths to find the velocity at a give location.

**Arguments:**
- wall panel strengths
- wake panel strengths
- point to probe

**Returns:**
- velocity at point

**Calls:**
- nothing

---

#### _`blade_section_sigma`_

Calculates source strength.

**Arguments:**
- inflow velocity
- section chord
- section radial location
- number of rotor blades
- section drag coefficient

**Returns:**

**Calls:**

---

### `solve_system`

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

---

FOR NEWTON SOLVE

### `residual`

Calculate rotor circulation and source panel strengths as a function of wall panel strengths.

**Arguments:**

**Returns**
- blade element circulations
- rotor source panel strengths

Calls:
- all the update functions below.

---

#### _`update_velocities`_

Description

**Arguments:**

**Returns:**

**Calls:**

---

#### _`update_blade_section_data`_

Description

**Arguments:**

**Returns:**

**Calls:**

---

#### _`update_wake_grid_data`_

Description

**Arguments:**

**Returns:**

**Calls:**

---

#### _`update_wake_panel_strengths`_

Description

**Arguments:**

**Returns:**

**Calls:**

#### _`update_source_panel_strengths`_

Description

**Arguments:**

**Returns:**

**Calls:**

---
