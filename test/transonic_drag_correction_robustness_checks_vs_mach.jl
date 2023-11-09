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

mach = 0.3
cdglob = 0.05
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
#             Out v Mach            #
#---------------------------------#
##### ----- Out v Mach, Sweep mach's ----- #####
for (is, cl) in enumerate(cls)
    # output vs cl at sweep of machs
    for (i, mach) in enumerate(machs)
        clnom = dt.transonicliftlimiter(
            clglob, mach, clcdminglob, clmaxglob, clminglob, dcldaglob
        )
        clsmooth = dt.transonicliftlimitersmooth!(
            [copy(clglob)], mach, clcdminglob, clmaxglob, clminglob, dcldaglob
        )[1]
        nom[i] = dt.transonicdragaddition(cdglob, clnom, clcdminglob, mach;)
        smooth[i] = dt.transonicdragadditionsmooth!(
            [copy(cdglob)], clsmooth, clcdminglob, mach;
        )[1]
    end
    titletext = @sprintf "%1.3f" cl
    plot(;
        xlabel="Mach Number",
        ylabel=L"c_{d_\mathrm{trans}}",
        title="cl = " * titletext * "(#$(is))",
        ylim=(0.0, 0.5),
    )
    plot!(machs, nom; color=myblue, label="Nominal")
    plot!(machs, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovmach-sweepcl_drag." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v Mach, Sweep Mcrits's ----- #####
for (is, mcrit) in enumerate(mcrits)
    # output vs cl at sweep of machs
    for (i, mach) in enumerate(machs)
        clnom = dt.transonicliftlimiter(
            clglob, mach, clcdminglob, clmaxglob, clminglob, dcldaglob; mcrit=mcrit
        )
        clsmooth = dt.transonicliftlimitersmooth!(
            [copy(clglob)], mach, clcdminglob, clmaxglob, clminglob, dcldaglob; mcrit=mcrit
        )[1]
        nom[i] = dt.transonicdragaddition(cdglob, clnom, clcdminglob, mach; mcrit=mcrit)
        smooth[i] = dt.transonicdragadditionsmooth!(
            [copy(cdglob)], clsmooth, clcdminglob, mach; mcrit=mcrit
        )[1]
    end
    titletext = @sprintf "%1.3f" mcrit
    plot(;
        xlabel="Mach Number",
        ylabel=L"c_{d_\mathrm{trans}}",
        title="Mcrit = " * titletext * "(#$(is))",
        ylim=(0.0, 0.5),
    )
    plot!(machs, nom; color=myblue, label="Nominal")
    plot!(machs, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovmach-sweepmcrit_drag." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v Mach, Sweep clcdmin's ----- #####
for (is, clcdmin) in enumerate(clcdmins)
    # output vs cl at sweep of machs
    for (i, mach) in enumerate(machs)
        clnom = dt.transonicliftlimiter(
            clglob, mach, clcdmin, clmaxglob, clminglob, dcldaglob
        )
        clsmooth = dt.transonicliftlimitersmooth!(
            [copy(clglob)], mach, clcdmin, clmaxglob, clminglob, dcldaglob
        )[1]
        nom[i] = dt.transonicdragaddition(cdglob, clnom, clcdmin, mach;)
        smooth[i] = dt.transonicdragadditionsmooth!(
            [copy(cdglob)], clsmooth, clcdmin, mach;
        )[1]
    end
    titletext = @sprintf "%1.3f" clcdmin
    plot(;
        xlabel="Mach Number",
        ylabel=L"c_{d_\mathrm{trans}}",
        title="clcdmin = " * titletext * "(#$(is))",
        ylim=(0.0, 0.5),
    )
    plot!(machs, nom; color=myblue, label="Nominal")
    plot!(machs, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovmach-sweepclcdmin_drag." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v Mach, Sweep dclda's ----- #####
for (is, dclda) in enumerate(dcldas)
    # output vs cl at sweep of machs
    for (i, mach) in enumerate(machs)
        clnom = dt.transonicliftlimiter(
            clglob, mach, clcdminglob, clmaxglob, clminglob, dclda
        )
        clsmooth = dt.transonicliftlimitersmooth!(
            [copy(clglob)], mach, clcdminglob, clmaxglob, clminglob, dclda
        )[1]
        nom[i] = dt.transonicdragaddition(cdglob, clnom, clcdminglob, mach;)
        smooth[i] = dt.transonicdragadditionsmooth!(
            [copy(cdglob)], clsmooth, clcdminglob, mach;
        )[1]
    end
    titletext = @sprintf "%1.3f" dclda
    plot(;
        xlabel="Mach Number",
        ylabel=L"c_{d_\mathrm{trans}}",
        title="dclda = " * titletext * "(#$(is))",
        ylim=(0.0, 0.5),
    )
    plot!(machs, nom; color=myblue, label="Nominal")
    plot!(machs, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovmach-sweepdclda_drag." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v Mach, Sweep dcldastall's ----- #####
for (is, dcldastall) in enumerate(dcldastalls)
    # output vs cl at sweep of machs
    for (i, mach) in enumerate(machs)
        clnom = dt.transonicliftlimiter(
            clglob,
            mach,
            clcdminglob,
            clmaxglob,
            clminglob,
            dcldaglob;
            dclda_stall=dcldastall,
        )
        clsmooth = dt.transonicliftlimitersmooth!(
            [copy(clglob)],
            mach,
            clcdminglob,
            clmaxglob,
            clminglob,
            dcldaglob;
            dclda_stall=dcldastall,
        )[1]
        nom[i] = dt.transonicdragaddition(cdglob, clnom, clcdminglob, mach;)
        smooth[i] = dt.transonicdragadditionsmooth!(
            [copy(cdglob)], clsmooth, clcdminglob, mach;
        )[1]
    end
    titletext = @sprintf "%1.3f" dcldastall
    plot(;
        xlabel="Mach Number",
        ylabel=L"c_{d_\mathrm{trans}}",
        title="dcldastall = " * titletext * "(#$(is))",
        ylim=(0.0, 0.5),
    )
    plot!(machs, nom; color=myblue, label="Nominal")
    plot!(machs, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovmach-sweepdcldastall_drag." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v Mach, Sweep dclstall's ----- #####
for (is, dclstall) in enumerate(dclstalls)
    # output vs cl at sweep of machs
    for (i, mach) in enumerate(machs)
        clnom = dt.transonicliftlimiter(
            clglob, mach, clcdminglob, clmaxglob, clminglob, dcldaglob; dcl_stall=dclstall
        )
        clsmooth = dt.transonicliftlimitersmooth!(
            [copy(clglob)],
            mach,
            clcdminglob,
            clmaxglob,
            clminglob,
            dcldaglob;
            dcl_stall=dclstall,
        )[1]
        nom[i] = dt.transonicdragaddition(cdglob, clnom, clcdminglob, mach;)
        smooth[i] = dt.transonicdragadditionsmooth!(
            [copy(cdglob)], clsmooth, clcdminglob, mach;
        )[1]
    end
    titletext = @sprintf "%1.3f" dclstall
    plot(;
        xlabel="Mach Number",
        ylabel=L"c_{d_\mathrm{trans}}",
        title="dclstall = " * titletext * "(#$(is))",
        ylim=(0.0, 0.5),
    )
    plot!(machs, nom; color=myblue, label="Nominal")
    plot!(machs, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovmach-sweepdclstall_drag." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v Mach, Sweep clmax's ----- #####
for (is, clmax) in enumerate(clmaxs)
    # output vs cl at sweep of machs
    for (i, mach) in enumerate(machs)
        clnom = dt.transonicliftlimiter(
            clglob, mach, clcdminglob, clmax, clminglob, dcldaglob
        )
        clsmooth = dt.transonicliftlimitersmooth!(
            [copy(clglob)], mach, clcdminglob, clmax, clminglob, dcldaglob
        )[1]
        nom[i] = dt.transonicdragaddition(cdglob, clnom, clcdminglob, mach;)
        smooth[i] = dt.transonicdragadditionsmooth!(
            [copy(cdglob)], clsmooth, clcdminglob, mach;
        )[1]
    end
    titletext = @sprintf "%1.3f" clmax
    plot(;
        xlabel="Mach Number",
        ylabel=L"c_{d_\mathrm{trans}}",
        title="clmax = " * titletext * "(#$(is))",
        ylim=(0.0, 0.5),
    )
    plot!(machs, nom; color=myblue, label="Nominal")
    plot!(machs, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovmach-sweepclmax_drag." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v Mach, Sweep clmin's ----- #####
for (is, clmin) in enumerate(clmins)
    # output vs cl at sweep of machs
    for (i, mach) in enumerate(machs)
        clnom = dt.transonicliftlimiter(
            clglob, mach, clcdminglob, clmaxglob, clmin, dcldaglob
        )
        clsmooth = dt.transonicliftlimitersmooth!(
            [copy(clglob)], mach, clcdminglob, clmaxglob, clmin, dcldaglob
        )[1]
        nom[i] = dt.transonicdragaddition(cdglob, clnom, clcdminglob, mach;)
        smooth[i] = dt.transonicdragadditionsmooth!(
            [copy(cdglob)], clsmooth, clcdminglob, mach;
        )[1]
    end
    titletext = @sprintf "%1.3f" clmin
    plot(;
        xlabel="Mach Number",
        ylabel=L"c_{d_\mathrm{trans}}",
        title="clmin = " * titletext * "(#$(is))",
        ylim=(0.0, 0.5),
    )
    plot!(machs, nom; color=myblue, label="Nominal")
    plot!(machs, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovmach-sweepclmin_drag." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end
