
#---------------------------------#
#       Load Packages, etc.       #
#---------------------------------#

# - Get Project Directory - #
project_dir = dirname(dirname(@__FILE__))
if project_dir == ""
    project_dir = "."
end

# - load DuctAPE - #
using DuctAPE
const dt = DuctAPE
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
clcdmin = clo[findmin(cdo)[2]]
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
#             Out v CLCDMIN            #
#---------------------------------#
##### ----- Out v CLCDMIN, Sweep Mach's ----- #####
for (is, mach) in enumerate(machs)

    # output vs cl at sweep of machs
    for (i, clcdmin) in enumerate(clcdmins)
        # println("Mach= ", mach)

        nom[i] = dt.transonicliftlimiter(
            clglob, mach, clcdmin, clmaxglob, clminglob, dcldaglob
        )

        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(clglob)], mach, clcdmin, clmaxglob, clminglob, dcldaglob
        )[1]
    end
    titletext = @sprintf "%1.3f" mach
    plot(;
        xlabel="clcdmin",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="Mach = " * titletext * "(#$(is))",
        ylim=(-0.75, 1.0),
    )
    plot!(clcdmins, nom; color=myblue, label="Nominal")
    plot!(clcdmins, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovclcdmin-sweepmach." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v CLCDMIN, Sweep Mcrits's ----- #####
for (is, mcrit) in enumerate(mcrits)

    # output vs cl at sweep of machs
    for (i, clcdmin) in enumerate(clcdmins)
        # println("Mach= ", mach)

        nom[i] = dt.transonicliftlimiter(
            clglob, machglob, clcdmin, clmaxglob, clminglob, dcldaglob; mcrit=mcrit
        )

        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(clglob)], machglob, clcdmin, clmaxglob, clminglob, dcldaglob; mcrit=mcrit
        )[1]
    end
    titletext = @sprintf "%1.3f" mcrit
    plot(;
        xlabel="clcdmin",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="Mcrit = " * titletext * "(#$(is))",
        ylims=(clminglob, clmaxglob),
    )
    plot!(clcdmins, nom; color=myblue, label="Nominal")
    plot!(clcdmins, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovclcdmin-sweepmcrit." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v CLCDMIN, Sweep cl's ----- #####
for (is, cl) in enumerate(cls)
    # output vs cl at sweep of machs
    for (i, clcdmin) in enumerate(clcdmins)
        nom[i] = dt.transonicliftlimiter(
            clglob, machglob, clcdmin, clmaxglob, clminglob, dcldaglob
        )
        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(clglob)], machglob, clcdmin, clmaxglob, clminglob, dcldaglob
        )[1]
    end
    titletext = @sprintf "%1.3f" cl
    plot(;
        xlabel="clcdmin",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="cl = " * titletext * "(#$(is))",
        ylims=(clminglob, clmaxglob),
    )
    plot!(clcdmins, nom; color=myblue, label="Nominal")
    plot!(clcdmins, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovclcdmin-sweepcl." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v CLCDMIN, Sweep dclda's ----- #####
for (is, dclda) in enumerate(dcldas)
    # output vs cl at sweep of machs
    for (i, clcdmin) in enumerate(clcdmins)
        nom[i] = dt.transonicliftlimiter(
            clglob, machglob, clcdmin, clmaxglob, clminglob, dclda
        )
        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(clglob)], machglob, clcdmin, clmaxglob, clminglob, dclda
        )[1]
    end
    titletext = @sprintf "%1.3f" dclda
    plot(;
        xlabel="clcdmin",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="dclda = " * titletext * "(#$(is))",
        ylims=(-1.0, clmaxglob),
    )
    plot!(clcdmins, nom; color=myblue, label="Nominal")
    plot!(clcdmins, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovclcdmin-sweepdclda." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v CLCDMIN, Sweep dcldastall's ----- #####
for (is, dcldastall) in enumerate(dcldastalls)
    # output vs cl at sweep of machs
    for (i, clcdmin) in enumerate(clcdmins)
        nom[i] = dt.transonicliftlimiter(
            clglob,
            machglob,
            clcdmin,
            clmaxglob,
            clminglob,
            dcldaglob;
            dclda_stall=dcldastall,
        )
        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(clglob)],
            machglob,
            clcdmin,
            clmaxglob,
            clminglob,
            dcldaglob;
            dclda_stall=dcldastall,
        )[1]
    end
    titletext = @sprintf "%1.3f" dcldastall
    plot(;
        xlabel="clcdmin",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="dcldastall = " * titletext * "(#$(is))",
        ylim=(-1.0, clmaxglob),
    )
    plot!(clcdmins, nom; color=myblue, label="Nominal")
    plot!(clcdmins, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovclcdmin-sweepdcldastall." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v CLCDMIN, Sweep dclstall's ----- #####
for (is, dclstall) in enumerate(dclstalls)
    # output vs cl at sweep of machs
    for (i, clcdmin) in enumerate(clcdmins)
        nom[i] = dt.transonicliftlimiter(
            clglob, machglob, clcdmin, clmaxglob, clminglob, dcldaglob; dcl_stall=dclstall
        )
        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(clglob)],
            machglob,
            clcdmin,
            clmaxglob,
            clminglob,
            dcldaglob;
            dcl_stall=dclstall,
        )[1]
    end
    titletext = @sprintf "%1.3f" dclstall
    plot(;
        xlabel="clcdmin",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="dclstall = " * titletext * "(#$(is))",
        ylim=(-1.0, clmaxglob),
    )
    plot!(clcdmins, nom; color=myblue, label="Nominal")
    plot!(clcdmins, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovclcdmin-sweepdclstall." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v CLCDMIN, Sweep clmax's ----- #####
for (is, clmax) in enumerate(clmaxs)
    # output vs cl at sweep of machs
    for (i, clcdmin) in enumerate(clcdmins)
        nom[i] = dt.transonicliftlimiter(
            clglob, machglob, clcdmin, clmax, clminglob, dcldaglob
        )
        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(clglob)], machglob, clcdmin, clmax, clminglob, dcldaglob
        )[1]
    end
    titletext = @sprintf "%1.3f" clmax
    plot(;
        xlabel="clcdmin",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="clmax = " * titletext * "(#$(is))",
        ylim=(-1, maximum(cls)),
    )
    plot!(clcdmins, nom; color=myblue, label="Nominal")
    plot!(clcdmins, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovclcdmin-sweepclmax." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end

##### ----- Out v CLCDMIN, Sweep clmin's ----- #####
for (is, clmin) in enumerate(clmins)
    # output vs cl at sweep of machs
    for (i, clcdmin) in enumerate(clcdmins)
        nom[i] = dt.transonicliftlimiter(
            clglob, machglob, clcdmin, clmaxglob, clmin, dcldaglob
        )
        smooth[i] = dt.transonicliftlimitersmooth!(
            [copy(clglob)], machglob, clcdmin, clmaxglob, clmin, dcldaglob
        )[1]
    end
    titletext = @sprintf "%1.3f" clmin
    plot(;
        xlabel="clcdmin",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="clmin = " * titletext * "(#$(is))",
        ylim=(minimum(cls), clmaxglob),
    )
    plot!(clcdmins, nom; color=myblue, label="Nominal")
    plot!(clcdmins, smooth; coor=myred, linestyle=:dash, linewidth=2, label="Smooth")
    filename = "ovclcdmin-sweepclmin." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end
