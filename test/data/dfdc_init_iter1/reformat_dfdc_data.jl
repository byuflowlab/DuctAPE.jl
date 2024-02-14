#---------------------------------#
#              SETUP              #
#---------------------------------#
# - Parameters - #
B = 5 #there are 5 blades

# - read in index file - #
# in the case here, element 1 is the hub, 2 is the duct, 3 is the rotor, 4 is the hub wake, 14 is the duct wake, and the rest are the interior wake sheets.
include("element_indices.jl")

# get ranges of element idices
nidx = dfdc_eid[:, 2:3]
pidx = dfdc_eid[:, 4:5]
nidranges = nidx[:, 2] .- nidx[:, 1] .+ 1
pidranges = pidx[:, 2] .- pidx[:, 1] .+ 1

# - read in geometry file to help find indices - #
include("dfdc_geometry.jl")

# get wake sheet length
num_wake_z_panels = num_wake_z_nodes - 1

# get number of nodes and panels on bodies interface with wake
nhub_interface_nodes = num_wake_z_nodes - nidranges[4]
nduct_interface_nodes = num_wake_z_nodes - nidranges[end]
nhub_interface_panels = num_wake_z_panels - pidranges[4]
nduct_interface_panels = num_wake_z_panels - pidranges[end]

# - read in functions to reformat - #
include("reformatting_functions.jl")

#---------------------------------#
#             REFORMAT            #
#---------------------------------#

#### ----- control point velocities ----- #####
# - Right Before Iterations - #
include("initial_controlpoint_vzvr.jl")
iter0_controlpoint_vels = reformat_controlpoint_velocities(
    vzvr0, pidx, nhub_interface_panels, nduct_interface_panels
)
# vel object contains: (; duct_vel, hub_vel, wake_vel, rotor_vel)

# # - Beginning of First Iteration - #
# include("iter1_controlpoint_vzvr.jl")
# iter1_controlpoint_vels = reformat_controlpoint_velocities(
#     vzvr1, pidx, nhub_interface_panels, nduct_interface_panels
# )

# - Beginning of Second Iteration - #
include("iter2_controlpoint_vzvr.jl")
iter2_controlpoint_vels = reformat_controlpoint_velocities(
    vzvr2, pidx, nhub_interface_panels, nduct_interface_panels
)

#### ----- control point pressures ----- #####

# ##### ----- Gamr ----- #####
include("initial_bgam.jl")
Gamr0 = reformat_circulation(BGAM, B)
include("iter1_relaxed_bgam.jl")
Gamr1 = reformat_circulation(bgamr1, B)

##### ----- gamw ----- #####
# include("initial_gamw_in_solve.jl")
include("initial_gamw.jl")
gamw0 = reformat_gamw(GAMTH0, nidx, nhub_interface_nodes, nduct_interface_nodes)

# # gamw_est
# include("iter1_gamw.jl")
# gamw_est = reformat_gamw(GAMTH1, nidx, nhub_interface_nodes, nduct_interface_nodes)

include("iter1_relaxed_gamw.jl")
gamw1 = reformat_gamw(relaxed_GTH, nidx, nhub_interface_nodes, nduct_interface_nodes)

##### ----- sigr ----- #####
include("initial_sigr.jl")
# include("iter1_sigr.jl")
include("iter2_sigr.jl")

##### ----- blade element data ----- #####
#these are the values that give us (and include) the Gamr_est values
#note that DFDC doesn't actually use Wr_rotor values, just Wz and Wtheta
# include("iter1_extended_blade_element_values.jl")
# bev1 = reformat_extended_blade_elements(extended_blade_element_values2, B)

# ##### ----- wake velocities ----- #####
# include("iter1_wm_wake_panels.jl")
# Wm_wake = reformat_wake_panel_vels(VMAV, nhub_interface_panels, nduct_interface_panels)

# include("iter1_wm_wake_wm_avg.jl")
# Wm_avg = reformat_wmavg(wakevels_iter1[:, 2], nhub_interface_nodes, nduct_interface_nodes)

##### ----- Body Pressures ----- #####
include("initial_body_cp.jl")
duct_cpR0, hub_cpR0 = reformat_body_cp(body_cpR0, pidx)
include("initial_body_cp_withdeltas.jl")
duct_cpdel0, hub_cpdel0 = reformat_body_cp(body_cpR0, pidx)
include("iter1_body_cp.jl")
duct_cpR1, hub_cpR1 = reformat_body_cp(body_cpR1, pidx)
include("iter2_body_cp.jl")
duct_cpR2, hub_cpR2 = reformat_body_cp(body_cpR2, pidx)
include("iter2_body_cp_withdeltas.jl")
duct_cpdel2, hub_cpdel2 = reformat_body_cp(body_cpR2, pidx)

##### ----- Wake on body AIC's ----- #####
include("wake_on_body_coeffs.jl")
w2baic, w2baic_pcp = reformat_wake_aic(
    aicgth, nidx, pidx, nidranges, pidranges, nhub_interface_nodes, nduct_interface_nodes
)
wakerhs0 = reformat_rhs(aicgth, GAMTH0, pidx)

# ##### ----- Linear System RHS's ----- #####
include("initial_vinfrhs.jl")
vinfrhs0 = reformat_rhs(vinf0rhs, pidx)

##### ----- Linear System Solution ----- #####
include("initial_gamb.jl")
gamb0 = reformat_sol(res0[:, 1], nidx)

##### ----- Rotor on body AIC's ----- #####
include("rotor_on_body_coeffs.jl")
r2baic, r2baic_pcp = reformat_rotor_aic(aicsigr, nidx, pidx, nidranges, pidranges)

##### ----- Body induced unit velocities ----- #####
include("body_aicz.jl")
include("body_aicr.jl")
b2bvhat, b2rvhat, b2wvhat = reformat_body_vhat(
    bodyaicz, bodyaicr, nidx, pidx, pidranges, nhub_interface_panels, nduct_interface_panels
)

##### ----- Rotor induced unit velocities ----- #####
include("rotor_aicz.jl")
include("rotor_aicr.jl")
r2bvhat, r2rvhat, r2wvhat = reformat_rotor_vhat(
    rotoraicz,
    rotoraicr,
    nidx,
    nidranges,
    pidx,
    pidranges,
    nhub_interface_panels,
    nduct_interface_panels,
)

##### ----- Wake induced unit velocities ----- #####
include("wake_aicz.jl")
include("wake_aicr.jl")
w2bvhat, w2rvhat, w2wvhat = reformat_wake_vhat(
    wakeaicz,
    wakeaicr,
    nidx,
    nidranges,
    pidx,
    pidranges,
    nhub_interface_panels,
    nduct_interface_panels,
)
