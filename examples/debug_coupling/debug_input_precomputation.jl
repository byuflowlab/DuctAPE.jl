# number of rotors
nrotor = length(rotor_parameters)

#------------------------------------#
# Discretize Wake and Repanel Bodies #
#------------------------------------#
rotoronly = false
nohub = false
noduct = false
if hub_coordinates == nothing && duct_coordinates != nothing
    nohub = true
    debug_duct_coordinates = duct_coordinates
    debug_hub_coordinates = [
        minimum(duct_coordinates[:, 1]) 0.0
        maximum(duct_coordinates[:, 1]) 0.0
    ]
elseif hub_coordinates != nothing && duct_coordinates == nothing
    noduct = true
    debug_hub_coordinates = hub_coordinates
    debug_duct_coordinates = [
        minimum(hub_coordinates[:, 1]) 0.0
        maximum(hub_coordinates[:, 1]) 0.0
    ]
elseif hub_coordinates == nothing && duct_coordinates == nothing
    rotoronly = true
end

# - Discretize Wake x-coordinates - #
# also returns indices of rotor locations in the wake
xwake, rotor_indices = dt.discretize_wake(
    debug_duct_coordinates,
    debug_hub_coordinates,
    rotor_parameters.xrotor,
    paneling_constants.wake_length,
    paneling_constants.npanels,
)

# - Repanel Bodies - #
rp_duct_coordinates, rp_hub_coordinates = dt.update_body_geometry(
    debug_duct_coordinates,
    debug_hub_coordinates,
    xwake,
    paneling_constants.nhub_inlet,
    paneling_constants.nduct_inlet;
    finterp=fm.akima,
)

#-----------------------------------#
# Position Duct and Get Rotor Radii #
#-----------------------------------#
# check that tip gap isn't too small
# set vector of tip gaps to zero for initialization (any overrides to zero thus taken c of automatically)
tip_gaps = zeros(eltype(rotor_parameters.tip_gap), nrotor)
if rotor_parameters[1].tip_gap != 0.0
    if rotor_parameters[1].tip_gap < 1e-4
        @warn "You have selected a tip gap for the foremost rotor that is smaller than 1e-4. Overriding to 0.0 to avoid near singularity issues."
    else
        tip_gaps[1] = rotor_parameters[1].tip_gap
    end
end

# can't have non-zero tip gaps for aft rotors
for ir in 2:nrotor
    if rotor_parameters[ir].tip_gap != 0.0
        @warn "DuctTAPE does not currently have capabilities for adding tip gap to any but the foremost rotor. Overriding to 0.0."
    else
        tip_gaps[ir] = rotor_parameters[ir].tip_gap
    end
end

# if hub was nothing, set hub radius to dimensional inner rotor radius
if nohub
    rp_hub_coordinates[:, 2] .= rotor_parameters.Rhub
end
t_duct_coordinates, Rtips, Rhubs = dt.place_duct(
    rp_duct_coordinates,
    rp_hub_coordinates,
    rotor_parameters[1].Rtip,
    tip_gaps,
    rotor_parameters.xrotor,
)

# generate body paneling
if nohub
    body_panels = dt.generate_body_panels(t_duct_coordinates, nothing)
else
    body_panels = dt.generate_body_panels(t_duct_coordinates, rp_hub_coordinates)
end

#----------------------------------#
# Generate Discretized Wake Sheets #
#----------------------------------#
# get discretization of wakes at leading rotor position

for i in 1:nrotor
    @assert Rtips[i] > Rhubs[i] "Rotor #$i Tip Radius is set to be less than its Hub Radius."
end

#rotor panel edges
rpe = range(Rhubs[1], Rtips[1]; length=paneling_constants.nwake_sheets)

# wake sheet starting radius including dummy sheets for tip gap.
if tip_gaps[1] == 0.0
    rwake = rpe
else
    rwake = [rpe; Rtips[1] + tip_gaps[1]]
end

# Initialize wake "grid"
xgrid, rgrid = dt.initialize_wake_grid(t_duct_coordinates, rp_hub_coordinates, xwake, rwake)

# Relax "Grid"
dt.relax_grid!(xgrid, rgrid; max_iterations=100, tol=1e-9, verbose=false)

# generate wake sheet paneling
wake_vortex_panels = dt.generate_wake_panels(
    xgrid[:, 1:length(rpe)], rgrid[:, 1:length(rpe)]
)

#------------------------------------------#
# Generate Rotor Panels and Blade Elements #
#------------------------------------------#

# rotor source panel objects
rotor_source_panels = [
    dt.generate_rotor_panels(
        rotor_parameters[i].xrotor, rgrid[rotor_indices[i], 1:length(rpe)]
    ) for i in 1:nrotor
]

# rotor blade element objects
blade_elements = [
    dt.generate_blade_elements(
        rotor_parameters[i].B,
        rotor_parameters[i].Omega,
        rotor_parameters[i].xrotor,
        rotor_parameters[i].r,
        rotor_parameters[i].chords,
        rotor_parameters[i].twists,
        rotor_parameters[i].airfoils,
        Rtips[i],
        Rhubs[i],
        rotor_source_panels[i].panel_center[:, 2],
    ) for i in 1:nrotor
]

#--------------------------------#
# Generate Relational Geometries #
#--------------------------------#

# body to body
# body_meshff = ff.generate_mesh(body_method, body_panels)
mesh_bb = dt.generate_one_way_mesh(body_panels, body_panels)

# body to rotor
mesh_rb = [
    dt.generate_one_way_mesh(body_panels, rotor_source_panels[i]) for
    i in 1:length(rotor_source_panels), j in 1:1
]

# rotor to body
mesh_br = [
    dt.generate_one_way_mesh(rotor_source_panels[j], body_panels) for i in 1:1,
    j in 1:length(rotor_source_panels)
]

# rotor to rotor
# note: broadcasting like this throws an error. so using comprehension instead
# mesh_rr = generate_one_way_mesh.(rotor_source_panels, rotor_source_panels')
mesh_rr = [
    dt.generate_one_way_mesh(rotor_source_panels[j], rotor_source_panels[i]) for
    i in 1:length(rotor_source_panels), j in 1:length(rotor_source_panels)
]

# wake to body
mesh_bw = [
    dt.generate_one_way_mesh(wake_vortex_panels[j], body_panels) for i in 1:1,
    j in 1:length(wake_vortex_panels)
]

# wake to rotor
# mesh_rw = generate_one_way_mesh.(wake_vortex_panels, rotor_source_panels')
mesh_rw = [
    dt.generate_one_way_mesh(wake_vortex_panels[j], rotor_source_panels[i]) for
    i in 1:length(rotor_source_panels), j in 1:length(wake_vortex_panels)
]

#---------------------------------#
# Calculate Coefficient Matrices  #
#---------------------------------#

##### ----- Induced Velcocities on Bodies ----- #####

# - body to body - #
# A_bbff = ff.assemble_ring_vortex_matrix_raw(ff.Constant(), [false; true], body_panels, body_meshff)
A_bb = dt.assemble_induced_velocity_on_body_matrix(
    mesh_bb, body_panels, body_panels; singularity="vortex"
)

# apply back-diagonal correction to duct portions of coefficient matrix
dt.apply_back_diagonal_correction!(
    A_bb, body_panels[1], mesh_bb.affect_panel_indices[1], mesh_bb.mesh2panel_affect
)

# - freestream to body - #
# b_bfff = ff.assemble_ring_boundary_conditions_raw(
#     ff.Dirichlet(), [false; true], body_panels, body_meshff
# )
b_bf =
    freestream.Vinf .* dt.assemble_body_freestream_boundary_conditions(body_panels, mesh_bb)

# - rotor to body - #
A_br = [
    dt.assemble_induced_velocity_on_body_matrix(
        mesh_br[i, j], [rotor_source_panels[j]], body_panels; singularity="source"
    ) for i in 1:1, j in 1:length(rotor_source_panels)
]

# - wake to body - #
A_bw = [
    dt.assemble_induced_velocity_on_body_matrix(
        mesh_bw[i, j], [wake_vortex_panels[j]], body_panels; singularity="vortex"
    ) for i in 1:1, j in 1:length(wake_vortex_panels)
]

##### ----- Induced Velcocities on Rotors ----- #####
# - rotor to rotor - #
A_rr = [
    dt.assemble_induced_velocity_matrices(
        mesh_rr[i, j], rotor_source_panels[j], rotor_source_panels[i]; singularity="source"
    ) for i in 1:length(rotor_source_panels), j in 1:length(rotor_source_panels)
]

# axial components
vx_rr = [
    A_rr[i, j][1] for i in 1:length(rotor_source_panels), j in 1:length(rotor_source_panels)
]

# radial components
vr_rr = [
    A_rr[i, j][2] for i in 1:length(rotor_source_panels), j in 1:length(rotor_source_panels)
]

# - body to rotor - #
A_rb = [
    dt.assemble_induced_velocity_matrices(
        mesh_rb[i, j], body_panels, rotor_source_panels[i]
    ) for i in 1:length(rotor_source_panels), j in 1:1
]

# axial components
vx_rb = [A_rb[i, j][1] for i in 1:length(rotor_source_panels), j in 1:1]

# radial components
vr_rb = [A_rb[i, j][2] for i in 1:length(rotor_source_panels), j in 1:1]

# - wake to rotor - #
A_rw = [
    dt.assemble_induced_velocity_matrices(
        mesh_rw[i, j], wake_vortex_panels[j], rotor_source_panels[i]
    ) for i in 1:length(rotor_source_panels), j in 1:length(wake_vortex_panels)
]

# axial components
vx_rw = [
    A_rw[i, j][1] for i in 1:length(rotor_source_panels), j in 1:length(wake_vortex_panels)
]

# radial components
vr_rw = [
    A_rw[i, j][2] for i in 1:length(rotor_source_panels), j in 1:length(wake_vortex_panels)
]

## -- Miscellaneous Values for Indexing -- ##

# - Get rotor panel edges and centers - #
rotor_panel_edges = [rgrid[rotor_indices[i], 1:length(rpe)] for i in 1:nrotor]
rotor_panel_centers = [rotor_source_panels[i].panel_center[:, 2] for i in 1:nrotor]

rotor_panel_edges = reduce(hcat, rotor_panel_edges)
rotor_panel_centers = reduce(hcat, rotor_panel_centers)

# get the total number of vortex panels on the bodies
num_body_panels = length(b_bf)

inputs = (;
    converged=[false],
    # body_geometry, # body geometry
    # - rotors
    blade_elements, # blade elements
    num_rotors=nrotor,
    rotor_panel_edges,
    rotor_panel_centers,
    # panels
    rotor_indices,
    num_wake_x_panels=length(xwake) - 1,
    num_body_panels,
    # body_panels, # body paneling
    # rotor_source_panels, # rotor paneling
    # wake_vortex_panels, # wake paneling
    # - unit induced velocities (INCLUDING PANEL LENGTH)
    A_bb, # body to body
    b_bf, # freestream contribution to body boundary conditions
    A_br, # rotor to body (total)
    A_bw, # wake to body (total)
    vx_rb, # body to rotor (x-direction)
    vr_rb, # body to rotor (r-direction)
    vx_rr, # rotor to rotor (x-direction)
    vr_rr, # rotor to rotor ( r-direction)
    vx_rw, # wake to rotor (x-direction)
    vr_rw, # wake to rotor ( r-direction)
    # kutta condition indices
    kutta_idxs=dt.get_kutta_indices([false; true], mesh_bb),
    # operating conditions
    Vinf=freestream.Vinf, # freestream parameters
    # - Debugging/Plotting
    t_duct_coordinates,
    rp_hub_coordinates,
    body_panels,
    rotor_source_panels,
    wake_vortex_panels,
    mesh_bb,
    mesh_rb,
    mesh_br,
    mesh_bw,
    mesh_rw,
)

######################################################################
#                                                                    #
#                             PLOTTING                               #
#                                                                    #
######################################################################

plot(; aspectratio=1)

plot!(
    xrotor * ones(length(rotor_panel_edges)),
    rotor_panel_edges;
    color=mycolors[2],
    linewidth=0.25,
    markersize=0.5,
    markershape=:rect,
    label="",
)
plot!(
    inputs.rotor_source_panels[1].panel_center[:, 1],
    inputs.rotor_source_panels[1].panel_center[:, 2];
    color=mycolors[2],
    seriestype=:scatter,
    markersize=0.5,
    markershape=:circle,
    label="",
)

for iw in 1:nwake_sheets
    plot!(
        xgrid[:, iw],
        rgrid[:, iw];
        linewidth=0.25,
        markersize=0.5,
        markershape=:rect,
        color=mycolors[3],
        label="",
    )

    plot!(
        inputs.wake_vortex_panels[iw].panel_center[:, 1],
        inputs.wake_vortex_panels[iw].panel_center[:, 2];
        seriestype=:scatter,
        markersize=0.5,
        markershape=:circle,
        color=mycolors[3],
        label="",
    )
end
savefig("examples/debug_coupling/precomputed-rotor-wake-geometry.pdf")

# for ib in 1:2
for ib in 1:1
    plot!(
        inputs.body_panels[ib].panel_center[:, 1],
        inputs.body_panels[ib].panel_center[:, 2];
        color=mycolors[1],
        label="",
    )
end

savefig("examples/debug_coupling/precomputed-full-geometry.pdf")
