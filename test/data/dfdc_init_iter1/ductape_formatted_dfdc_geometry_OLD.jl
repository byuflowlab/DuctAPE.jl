#=
This file takes the dfdc geometry written out directly from DFDC (we added that feature in our copy of DFDC) and generates the system_geometry object required by DuctAPE
=#

# - Load Geometry Written from DFDC - #
#= NOTES:
Duct Coordinates are given from trailing edge, counter-clockwise
Hub Coordinates are given from trailing to leading edge (also counter-clockwise)
Wake Coordinates for inner wakes start at the wake nearest the hub and at the first rotor then proceed down the wake sheet axially before moving the wake sheet next closest to the duct, again starting from the foremost rotor
The Hub and Duct wakes are given from the Hub and Duct trailing edges
=#
include("dfdc_geometry.jl")

# - Add the coordinates of the body-wake interface to the Hub and Duct wakes - #
#= NOTES:
DuctAPE indexing requires the full rotor shed wake coordinates to be present, even though the "wake" that touches the body isn't actually a wake numerically.
=#

# Find the indices where the wakes start (based on rotor axial position at hub and tip)
dii = searchsortedlast(dfdc_duct_coordinates[:, 1], dfdc_rotor_coordinates[end, 1])
hii = findfirst(x -> x == dfdc_rotor_coordinates[1, 1], dfdc_hub_coordinates[:, 1])

# Concatenate the splice of the body coordinates and the body wake coordinates, remembering to exclude the repeated value at the trailing edge junction
hub_interface_wake = [
    reverse(dfdc_hub_coordinates[1:hii, :]; dims=1)
    dfdc_hub_wake_coordinates[2:end, :]
]
hubTE_index = length(dfdc_hub_coordinates[1:hii, 1])

duct_interface_wake = [
    dfdc_duct_coordinates[dii:end, :]
    dfdc_duct_wake_coordinates[2:end, :]
]
ductTE_index = length(dfdc_duct_coordinates[dii:end, 1])

# - Combine all the wake nodes together - #
#= NOTES:
DuctAPE takes in the wake coordinates in a 3D Array with dimensions of:
dim 1: z,r
dim 2: axial coordinates for each node
dim 3: radial coordinates for each node
=#

# Initialize the wake grid object
grid = zeros(2, num_wake_z_nodes, num_inner_wake_sheets + 2)

# Reshape the inner wake coordinates and fill in the inner grid nodes
grid[1, :, 2:(end - 1)] .= reshape(
    dfdc_wake_coordinates[:, 1], num_wake_z_nodes, num_inner_wake_sheets
)
grid[2, :, 2:(end - 1)] .= reshape(
    dfdc_wake_coordinates[:, 2], num_wake_z_nodes, num_inner_wake_sheets
)

# Fill in the extended duct and hub wake nodes
grid[:, :, 1] .= hub_interface_wake'
grid[:, :, end] .= duct_interface_wake'

# - Flip the Body Coordinates to go clockwise - #
#= NOTES;
The DuctAPE Panel Method is formulated with the geometry defined in a clockwise manner
=#

duct_coordinates = reverse(dfdc_duct_coordinates; dims=1)' .* 1.0
hub_coordinates = reverse(dfdc_hub_coordinates; dims=1)' .* 1.0

system_geometry = (;
    duct_coordinates,
    hub_coordinates,
    grid,
    Rtips=[dfdc_duct_coordinates[dii, 2]],
    Rhubs=[dfdc_hub_coordinates[hii, 2]],
    rpe=dfdc_rotor_coordinates[:, 2],
    zwake=grid[1, :, 1],
    rwake=dfdc_rotor_coordinates[:, 2],
    rotor_indices_in_wake=[1],
    ductTE_index,
    hubTE_index,
    nohub=false,
    noduct=false,
    rotoronly=false,
)

