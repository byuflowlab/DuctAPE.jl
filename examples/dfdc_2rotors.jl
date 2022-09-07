using DuctTAPE
using Plots
include("../plots_default.jl")

#include dfdc case file data
include("../data/dfdc/dstestr2_case.jl");

#include rotor setup
include("rotor_setup.jl")

#include geometry setup
include("geometry_setup.jl")

#include wake grid setup
include("wakegrid_setup.jl")

#include paneling setup
include("panel_setup.jl")

#include rotor aero setup
include("rotorgrid_aero_setup.jl")

### --- Run Examples --- ###

# Set up Rotors
rotors = setup_rotors()

# Set up Duct Geometery
ductgeometry, ductsplines = setup_geometry(; plotgeometry=true)

# Set up Wake Grid
wakegrid = setup_wakegrid(ductgeometry, ductsplines, rotors; plotgrid=true)

# Get Blade Dimensions
blades = setup_blades(ductgeometry, ductsplines, rotors)

# Set up System Panels
wallpanels, hubpanels, rotor_source_panels, wakepanels = setup_panels(
    ductgeometry, ductsplines, rotors, wakegrid; plotpanels=true
)

# Initialize System Aerodynamics (Rotor and effects on Grid)
system_aero, rotor_velocities, average_axial_velocity = setup_rotorgrid_aero(
    rotors, blades, wakegrid, rotor_source_panels
)

