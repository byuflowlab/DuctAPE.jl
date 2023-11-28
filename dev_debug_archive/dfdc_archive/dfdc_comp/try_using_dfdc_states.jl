project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

include(project_dir * "/plots_default.jl")

savepath = project_dir * "/dev_debug_archive/dfdc_comp/"

include(savepath * "gather_dfdc_data.jl")
include(savepath * "run_dfdc_geom.jl")

using DuctAPE
const dt = DuctAPE

dfdc = get_dfdc()

duct_coordinates, hub_coordinates, wake_coordinates, rotor_parameters, freestream, reference_parameters = init(
    dfdc
)

inputs = dt.manual_precomputed_inputs(
    duct_coordinates,      # panel node locations
    hub_coordinates,       # panel node locations
    wake_coordinates,      # tuple containing xgrid[x,r], rgrid[x,r]
    rotor_parameters,      # tuple with rotor paramters
    freestream,            # tuple containing Vinf, rho, mu, asound
    reference_parameters;  # tuple containing Vref, Rref
    debug=false,
)

# put together all the wake strengths
gamw1 = [reverse(dfdc.gamw_hubinterface); dfdc.gamw_hub[2:end]]
gamwend = [dfdc.gamw_ductinterface; dfdc.gamw_duct[2:end]]
gamw_dfdc = [gamw1'; dfdc.gamw; gamwend']

# average (streamwise) the wake strengths to get things down to length 50
nr, nx = size(gamw_dfdc)
gamw_aprx = zeros(nr, nx - 1)
for i in 1:(length(gamw1) - 1)
    gamw_aprx[:, i] = 0.5 * (gamw_dfdc[:, i] .+ gamw_dfdc[:, i + 1])
end

# calculate wake-on-wake velocities using vx_ww's

vxa_wake = zeros(nr - 1, nx - 1)
for iplane in 1:(nx - 1)
    for jwake in 1:nr
        @views vxa_wake[:, iplane] .+= inputs.vx_ww[iplane, jwake] * gamw_aprx[jwake, :]
    end
end

vxbar = zeros(nr, nx - 1)
for iplane in 1:(nx - 1)
    vxbar[:, iplane] = dt.radially_average_velocity(vxa_wake[:, iplane], 1)
end

# put together all the wake velocities
#TODO:one too long
vxww1 = [reverse(dfdc.vx_hiw); dfdc.vx_hww]
vxwwend = [dfdc.vx_diw; dfdc.vx_dww]
vxww_dfdc = [vxww1'; dfdc.vx_ww; vxwwend']
# put together all the wake locations
#TODO:one too long
cp1 = [reverse(dfdc.hubinterface_ctrlpt_x); dfdc.hubwake_ctrlpt_x]
cpend = [dfdc.ductinterface_ctrlpt_x; dfdc.ductwake_ctrlpt_x]
wake_cp = [cp1'; dfdc.wake_ctrlpt_x; cpend']


# - plot wakes on wakes - #
vxww = plot(; xlabel="wake-induced axial velocity", ylabel="r")
plot!(vxww, dfdc.vx_rw, dfdc.rotor_ctrlpt_r; color=:black, label="DFDC on rotor")

for i in 1:8:length(dfdc.vx_ww[1, :])
    plot!(
        vxww,
        # vxbar[:, i],
        # (p -> p.panel_center[i, 2]).(inputs.wake_vortex_panels);
        vxa_wake[:, i],
        inputs.wake_affect_panels[i].panel_center[:,2];
        label="DT: $i stations behind rotor",
    )
    plot!(
        vxww,
        dfdc.vx_ww[:, i],
        dfdc.wake_ctrlpt_r[:, i];
        linestyle=:dash,
        label="DFDC: $i stations behind rotor",
    )
end

savefig(vxww, datapath * "wake-on-wake-vx_dfdcgamw_ducttapevxww.pdf")
