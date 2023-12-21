project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

include(project_dir * "/plots_default.jl")
datapath = project_dir * "/dev_debug_archive/dfdc_comp/"

include(datapath * "gather_dfdc_data.jl")
dfdc = get_dfdc()

function plot_vx(dfdc)

    # - velocity on rotor - #
    vxr = plot(; xlabel="induced axial velocity on rotor", ylabel="r")

    plot!(vxr, dfdc.vz_rb, dfdc.rotor_ctrlpt_r; label="from bodies")
    plot!(vxr, dfdc.vz_rr, dfdc.rotor_ctrlpt_r; label="from rotor")
    plot!(vxr, dfdc.vz_rw, dfdc.rotor_ctrlpt_r; label="from wake")
    plot!(vxr, dfdc.vz_rotor, dfdc.rotor_ctrlpt_r; label="total")

    savefig(vxr, datapath * "dfdc-vx-on-rotor.pdf")

    # - velocity on bodies - #
    vxb = plot(; ylabel="induced axial velocity", xlabel="x")
    plot!(vxb, dfdc.duct_ctrlpt_x, dfdc.vz_bb_duct; label="from bodies on duct")
    plot!(vxb, dfdc.duct_ctrlpt_x, dfdc.vz_br_duct; label="from rotor on duct")
    plot!(vxb, dfdc.duct_ctrlpt_x, dfdc.vz_bw_duct; label="from wake on duct")
    plot!(vxb, dfdc.duct_ctrlpt_x, dfdc.vz_duct; linewidth=2, label="total on duct")
    plot!(
        vxb, dfdc.hub_ctrlpt_x, dfdc.vz_bb_hub; linestyle=:dash, label="from bodies on hub"
    )
    plot!(
        vxb, dfdc.hub_ctrlpt_x, dfdc.vz_br_hub; linestyle=:dash, label="from rotor on hub"
    )
    plot!(vxb, dfdc.hub_ctrlpt_x, dfdc.vz_bw_hub; linestyle=:dash, label="from wake on hub")
    plot!(
        vxb,
        dfdc.hub_ctrlpt_x,
        dfdc.vz_hub;
        linewidth=2,
        linestyle=:dash,
        label="total on hub",
    )

    savefig(vxb, datapath * "dfdc-vx-on-bodies.pdf")

    # - body on wakes - #
    vxbw = plot(; xlabel="body-induced axial velocity", ylabel="r")
    plot!(vxbw, dfdc.vz_rb, dfdc.rotor_ctrlpt_r; color=:black, label="on rotor")
    for i in 1:8:length(dfdc.vz_wb[1, :])
        plot!(
            vxbw,
            dfdc.vz_wb[:, i],
            dfdc.wake_ctrlpt_r[:, i];
            label="$i stations behind rotor",
        )
    end
    savefig(vxbw, datapath * "dfdc-body-induced-vx-on-wake.pdf")

    # - body on wakes - #
    vxww = plot(; xlabel="wake-induced axial velocity", ylabel="r")
    plot!(vxww, dfdc.vz_rw, dfdc.rotor_ctrlpt_r; color=:black, label="on rotor")
    for i in 1:8:length(dfdc.vz_ww[1, :])
        plot!(
            vxww,
            dfdc.vz_ww[:, i],
            dfdc.wake_ctrlpt_r[:, i];
            label="$i stations behind rotor",
        )
    end
    savefig(vxww, datapath * "dfdc-wake-induced-vx-on-wake.pdf")

    return nothing
end

plot_vx(dfdc)
# plot_vr(dfdc)
# plot_Vm(dfdc)
# plot_Gamr(dfdc)
# plot_sigr(dfdc)
# plot_gamb(dfdc)
# plot_gamw(dfdc)

