
#---------------------------------#
#       Load Packages, etc.       #
#---------------------------------#

# - Get Project Directory - #
project_dir = dirname(dirname(@__FILE__))
if project_dir == ""
    project_dir = "."
end

# - load DuctTAPE - #
using DuctTAPE
const dt = DuctTAPE
using Printf

# - load plotting defaults - #
include(project_dir * "/visualize/plots_default.jl")
gr()
savepath = project_dir * "/test/figures/"

#---------------------------------#
#    Load Nominal Airfoil Data    #
#---------------------------------#
include(project_dir * "/test/data/naca_4412_raw.jl")
aoa = alpha
clo = cl
cdo = cd

#---------------------------------#
#  transonic limiter correction   #
#---------------------------------#

machglob = 0.3
clcdminglob = clo[findmin(cdo)[2]]
clmaxglob, maxid = findmax(clo)
clminglob, minid = findmin(clo)
dcldaglob = (clmaxglob - clminglob) / (aoa[maxid] * pi / 180 - aoa[minid] * pi / 180)

N = 501
cls = range(-2.0, 2.0, N)
machs = range(0.0, 1.25, N)
mcrits = range(0.0, 1.0, N)
clcdmins = range(-2.0, 2.0, N)
dcldas = range(0.01, 4 * pi, N)
dcldastalls = range(0.01, 4 * pi, N)
dclstalls = range(0.0, 5.0, N)
clmaxs = range(0.0, 5.0, N)
clmins = range(0.0, -5.0, N)

nom = zeros(N)
smooth = zeros(N)

#---------------------------------#
#             Out v CL            #
#---------------------------------#
##### ----- Out v CL, Sweep Mach's ----- #####
for (is, mach) in enumerate(machs)

    # output vs cl at sweep of machs
    for (i, cl) in enumerate(cls)
        # println("Mach= ", mach)

        nom[i] = dt.transonicliftlimiter(
            cl, mach, clcdminglob, clmaxglob, clminglob, dcldaglob
        )

        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(cl)], mach, clcdminglob, clmaxglob, clminglob, dcldaglob
        )[1]
    end
    titletext = @sprintf "%1.3f" mach
    plot(;
        xlabel=L"c_\ell",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="Mach = " * titletext * "(#$(is))",
        ylim=(-0.75, 1.0),
    )
    plot!(cls, nom; color=myblue, label="Nominal")
    plot!(cls, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovcl-sweepmach." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v CL, Sweep Mcrits's ----- #####
for (is, mcrit) in enumerate(mcrits)

    # output vs cl at sweep of machs
    for (i, cl) in enumerate(cls)
        # println("Mach= ", mach)

        nom[i] = dt.transonicliftlimiter(
            cl, machglob, clcdminglob, clmaxglob, clminglob, dcldaglob; mcrit=mcrit
        )

        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(cl)], machglob, clcdminglob, clmaxglob, clminglob, dcldaglob; mcrit=mcrit
        )[1]
    end
    titletext = @sprintf "%1.3f" mcrit
    plot(;
        xlabel=L"c_\ell",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="Mcrit = " * titletext * "(#$(is))",
        ylims=(clminglob, clmaxglob),
    )
    plot!(cls, nom; color=myblue, label="Nominal")
    plot!(cls, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovcl-sweepmcrit." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v CL, Sweep clcdmin's ----- #####
for (is, clcdmin) in enumerate(clcdmins)
    # output vs cl at sweep of machs
    for (i, cl) in enumerate(cls)
        nom[i] = dt.transonicliftlimiter(
            cl, machglob, clcdmin, clmaxglob, clminglob, dcldaglob
        )
        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(cl)], machglob, clcdmin, clmaxglob, clminglob, dcldaglob
        )[1]
    end
    titletext = @sprintf "%1.3f" clcdmin
    plot(;
        xlabel=L"c_\ell",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="clcdmin = " * titletext * "(#$(is))",
        ylims=(clminglob, clmaxglob),
    )
    plot!(cls, nom; color=myblue, label="Nominal")
    plot!(cls, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovcl-sweepclcdmin." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v CL, Sweep dclda's ----- #####
for (is, dclda) in enumerate(dcldas)
    # output vs cl at sweep of machs
    for (i, cl) in enumerate(cls)
        nom[i] = dt.transonicliftlimiter(
            cl, machglob, clcdminglob, clmaxglob, clminglob, dclda
        )
        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(cl)], machglob, clcdminglob, clmaxglob, clminglob, dclda
        )[1]
    end
    titletext = @sprintf "%1.3f" dclda
    plot(;
        xlabel=L"c_\ell",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="dclda = " * titletext * "(#$(is))",
        ylims=(-1.0, clmaxglob),
    )
    plot!(cls, nom; color=myblue, label="Nominal")
    plot!(cls, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovcl-sweepdclda." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v CL, Sweep dcldastall's ----- #####
for (is, dcldastall) in enumerate(dcldastalls)
    # output vs cl at sweep of machs
    for (i, cl) in enumerate(cls)
        nom[i] = dt.transonicliftlimiter(
            cl,
            machglob,
            clcdminglob,
            clmaxglob,
            clminglob,
            dcldaglob;
            dclda_stall=dcldastall,
        )
        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(cl)],
            machglob,
            clcdminglob,
            clmaxglob,
            clminglob,
            dcldaglob;
            dclda_stall=dcldastall,
        )[1]
    end
    titletext = @sprintf "%1.3f" dcldastall
    plot(;
        xlabel=L"c_\ell",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="dcldastall = " * titletext * "(#$(is))",
        ylim=(-1.0, clmaxglob)
    )
    plot!(cls, nom; color=myblue, label="Nominal")
    plot!(cls, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovcl-sweepdcldastall." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v CL, Sweep dclstall's ----- #####
for (is, dclstall) in enumerate(dclstalls)
    # output vs cl at sweep of machs
    for (i, cl) in enumerate(cls)
        nom[i] = dt.transonicliftlimiter(
            cl, machglob, clcdminglob, clmaxglob, clminglob, dcldaglob; dcl_stall=dclstall
        )
        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(cl)], machglob, clcdminglob, clmaxglob, clminglob, dcldaglob; dcl_stall=dclstall
        )[1]
    end
    titletext = @sprintf "%1.3f" dclstall
    plot(;
        xlabel=L"c_\ell",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="dclstall = " * titletext * "(#$(is))",
        ylim=(-1.0, clmaxglob),
    )
    plot!(cls, nom; color=myblue, label="Nominal")
    plot!(cls, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovcl-sweepdclstall." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v CL, Sweep clmax's ----- #####
for (is, clmax) in enumerate(clmaxs)
    # output vs cl at sweep of machs
    for (i, cl) in enumerate(cls)
        nom[i] = dt.transonicliftlimiter(
            cl, machglob, clcdminglob, clmax, clminglob, dcldaglob
        )
        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(cl)], machglob, clcdminglob, clmax, clminglob, dcldaglob
        )[1]
    end
    titletext = @sprintf "%1.3f" clmax
    plot(;
        xlabel=L"c_\ell",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="clmax = " * titletext * "(#$(is))",
        ylim=(-1, maximum(cls)),
    )
    plot!(cls, nom; color=myblue, label="Nominal")
    plot!(cls, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovcl-sweepclmax." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v CL, Sweep clmin's ----- #####
for (is, clmin) in enumerate(clmins)
    # output vs cl at sweep of machs
    for (i, cl) in enumerate(cls)
        nom[i] = dt.transonicliftlimiter(
            cl, machglob, clcdminglob, clmaxglob, clmin, dcldaglob
        )
        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(cl)], machglob, clcdminglob, clmaxglob, clmin, dcldaglob
        )[1]
    end
    titletext = @sprintf "%1.3f" clmin
    plot(;
        xlabel=L"c_\ell",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="clmin = " * titletext * "(#$(is))",
        ylim=(minimum(cls), clmaxglob),
    )
    plot!(cls, nom; color=myblue, label="Nominal")
    plot!(cls, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovcl-sweepclmin." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

# ##### ----- Out v Mach, Sweep cl's ----- #####
# for (icl, cl) in enumerate(cls)

#     # output vs cl at sweep of machs
#     for (i, mach) in enumerate(machs)
#         # println("Mach= ", mach)

#         nom[i] = dt.transonicliftlimiter(
#             cl, mach, clcdminglob, clmaxglob, clminglob, dcldaglob
#         )

#         smooth[i] = dt.transonicliftlimitersmooth!(
#             [copy(cl)], mach, clcdminglob, clmaxglob, clminglob, dcldaglob
#         )[1]
#     end
#     titletext = @sprintf "%1.3f" cl
#     plot(;
#         xlabel="Mach Number",
#         ylabel=L"c_{\ell_\mathrm{trans}}",
#         title="cl = " * titletext * "(#$(icl))",
#         ylim=(minimum(cls), maximum(cls)),
#     )
#     plot!(machs, nom; color=myblue, label="Nominal")
#     plot!(machs, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
#     filename = "ovmach-sweepcl." * lpad(icl, 4, "0") * ".png"
#     savefig(savepath * filename)
# end

# ##### ----- Out v Mcrit, Sweep cl's ----- #####
# for (icl, cl) in enumerate(cls)

#     # output vs cl at sweep of machs
#     for (i, mcrit) in enumerate(mcrits)
#         # println("Mach= ", mach)

#         nom[i] = dt.transonicliftlimiter(
#             cl, machglob, clcdminglob, clmaxglob, clminglob, dcldaglob;mcrit=mcrit
#         )

#         smooth[i] = dt.transonicliftlimitersmooth!(
#             [copy(cl)], machglob, clcdminglob, clmaxglob, clminglob, dcldaglob;mcrit=mcrit
#         )[1]
#     end
#     titletext = @sprintf "%1.3f" cl
#     plot(;
#         xlabel="Critical Mach Number",
#         ylabel=L"c_{\ell_\mathrm{trans}}",
#         title="cl = " * titletext * "(#$(icl))",
#         ylim=(minimum(cls), maximum(cls)),
#     )
#     plot!(machs, nom; color=myblue, label="Nominal")
#     plot!(machs, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
#     filename = "ovmcrit-sweepcl." * lpad(icl, 4, "0") * ".png"
#     savefig(savepath * filename)
# end

# # - generate gifs - #
# # convert -delay 2 -loop 0 ovmach-sweepcl_*.png ovmach_sweepcl.gif

# # output vs mach at sweep of mcrits
# # output vs mach at sweep of clcdmin
# # output vs mach at sweep of dclda
# # output vs mach at sweep of dcldastall
# # output vs mach at sweep of dclstall

# # output vs mcrit at sweep of clcdmin
# # output vs mcrit at sweep of dclda
# # output vs mcrit at sweep of dcldastall
# # output vs mcrit at sweep of dclstall

# # output vs dclda at sweep of clcdmin
# # output vs dclda at sweep of dcldastall
# # output vs dclda at sweep of dclstall

# # output vs clcdmin at sweep of dcldastall
# # output vs clcdmin at sweep of dclstall

# # output vs dcldastall at sweep of dclstall