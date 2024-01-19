# read in index file
include("dfdc_element_indices.jl")
# read in geometry file to help find indices
include("dfdc_geometry.jl")

# read in data files for which to reformat
include("dfdc_initial_values.jl")
include("dfdc_first_iteration.jl")

# in the case here, element 1 is the hub, 2 is the duct, 3 is the rotor, 4 is the hub wake, 14 is the duct wake, and the rest are the interior wake sheets.

# get ranges of element idices
nidx = dfdc_eid[:, 2:3]
pidx = dfdc_eid[:, 4:5]
nidranges = nidx[:, 2] .- nidx[:, 1]
pidranges = pidx[:, 2] .- pidx[:, 1]

# get wake sheet length
num_wake_z_panels = num_wake_z_nodes - 1

# get number of nodes and panels on bodies interface with wake
nhub_interface_nodes = num_wake_z_nodes - nidranges[4]
nduct_interface_nodes = num_wake_z_nodes - nidranges[end]
nhub_interface_panels = num_wake_z_panels - pidranges[4]
nduct_interface_panels = num_wake_z_panels - pidranges[end]

# fill in first and last wake sheet with zeros?
# these are the average velocities at the wake nodes
Wm_avg_init = [
    zeros(nhub_interface_nodes - 1) #start with zeros on hub wake interface, don't repeat TE node
    VMAVG_init[1:(end - nidranges[end] - 1)] #include hub wake and interior wakes
    zeros(nduct_interface_nodes - 1) # start with zeros in duct wake interface
    VMAVG_init[(end - nidranges[end]):end] # finish with duct wake, don't repeate TE node
]

B = 5 #there are 5 blades
Gamr_init = BGAM / B
sigr_init = zeros(length(BGAM) + 1)

# Format the Wake Strengths
hub_gamw = GAMTH[nidx[1, 1]:nidx[1, 2]]
duct_gamw = GAMTH[nidx[2, 1]:nidx[2, 2]]
hubwake_gamw = GAMTH[nidx[4, 1]:nidx[4, 2]]
interior_gamw = GAMTH[nidx[5, 1]:nidx[13, 2]]
ductwake_gamw = GAMTH[nidx[14, 1]:nidx[14, 2]]

gamw_init = [
    reverse(hub_gamw[2:nhub_interface_nodes]) #don't repeat the TE node
    hubwake_gamw
    interior_gamw
    duct_gamw[(end - nhub_interface_nodes + 3):end] #don't repeat the TE node
    ductwake_gamw
]

