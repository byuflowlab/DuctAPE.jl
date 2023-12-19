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
cdglob = 0.05
clglob= 0.5
clcdminglob = clo[findmin(cdo)[2]]
clmaxglob, maxid = findmax(clo)
clminglob, minid = findmin(clo)
dcldaglob = (clmaxglob - clminglob) / (aoa[maxid] * pi / 180 - aoa[minid] * pi / 180)
solidityglob = 1.0
staggerglob = pi / 7.0

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
solidities = range(0.0, 5.0, N)
staggers = range(0.0, 90*pi/180.0, N)

clout = zeros(N)
cdout = zeros(N)

##---------------------------------#
##             Out v Mach            #
##---------------------------------#
###### ----- Out v Mach, Sweep cl's ----- #####
#for (is, cl) in enumerate(cls)
#    for (i, mach) in enumerate(machs)
#        clout[i], cdout[i] = dt.corrected_clcd!(
#            [copy(cl)],
#            [copy(cdglob)],
#            mach,
#            solidityglob,
#            staggerglob,
#            clcdminglob,
#            clmaxglob,
#            clminglob,
#            dcldaglob,
#           )
#    end
#    titletext = @sprintf "%1.3f" cl
#    plot(;
#        xlabel="Mach Number",
#        ylabel=L"c_{\ell_\mathrm{trans}}",
#        title="cl = " * titletext * "(#$(is))",
#        ylim=(-2.0, 2.0),
#    )
#    plot!(machs, cl*ones(N); color=myblue, label="Nominal")
#    plot!(machs, clout; coor=myred, label="Corrected")
#    filename = "ovmach-sweepcl_corrcl." * lpad(is, 4, "0") * ".png"
#    savefig(savepath * filename)

#    plot(;
#        xlabel="Mach Number",
#        ylabel=L"c_{d_\mathrm{trans}}",
#        title="cl = " * titletext * "(#$(is))",
#        ylim=(0.0, 0.5),
#    )
#    plot!(machs, cdglob*ones(N); color=myblue, label="Nominal")
#    plot!(machs, cdout; coor=myred, label="Corrected")
#    filename = "ovmach-sweepcl_corrcd." * lpad(is, 4, "0") * ".png"
#    savefig(savepath * filename)
#end

###### ----- Out v Mach, Sweep cd's ----- #####
#for (is, cd) in enumerate(cds)
#    for (i, mach) in enumerate(machs)
#        clout[i], cdout[i] = dt.corrected_clcd!(
#            [copy(clglob)],
#            [copy(cd)],
#            mach,
#            solidityglob,
#            staggerglob,
#            clcdminglob,
#            clmaxglob,
#            clminglob,
#            dcldaglob,
#           )
#    end
#    titletext = @sprintf "%1.3f" cd
#    # plot(;
#    #     xlabel="Mach Number",
#    #     ylabel=L"c_{\ell_\mathrm{trans}}",
#    #     title="cd = " * titletext * "(#$(is))",
#    #     ylim=(-2.0, 2.0),
#    # )
#    # plot!(machs, cl*ones(N); color=myblue, label="Nominal")
#    # plot!(machs, clout; coor=myred, label="Corrected")
#    # filename = "ovmach-sweepcl_corrcl." * lpad(is, 4, "0") * ".png"
#    # savefig(savepath * filename)

#    plot(;
#        xlabel="Mach Number",
#        ylabel=L"c_{d_\mathrm{trans}}",
#        title="cd = " * titletext * "(#$(is))",
#        ylim=(0.0, 0.5),
#    )
#    plot!(machs, cd*ones(N); color=myblue, label="Nominal")
#    plot!(machs, cdout; coor=myred, label="Corrected")
#    filename = "ovmach-sweepcd_corrcd." * lpad(is, 4, "0") * ".png"
#    savefig(savepath * filename)
#end


## solidity
###### ----- Out v Mach, Sweep solidity's ----- #####
#for (is, solidity) in enumerate(solidities)
#    for (i, mach) in enumerate(machs)
#        clout[i], cdout[i] = dt.corrected_clcd!(
#            [copy(clglob)],
#            [copy(cdglob)],
#            machglob,
#            solidity,
#            staggerglob,
#            clcdminglob,
#            clmaxglob,
#            clminglob,
#            dcldaglob,
#           )
#    end
#    titletext = @sprintf "%1.3f" solidity
#    plot(;
#        xlabel="Mach Number",
#        ylabel=L"c_{\ell_\mathrm{trans}}",
#        title="solidity = " * titletext * "(#$(is))",
#        ylim=(-2.0, 2.0),
#    )
#    plot!(machs, clglob*ones(N); color=myblue, label="Nominal")
#    plot!(machs, clout; coor=myred, label="Corrected")
#    filename = "ovmach-sweepsolidity_corrcl." * lpad(is, 4, "0") * ".png"
#    savefig(savepath * filename)

#    plot(;
#        xlabel="Mach Number",
#        ylabel=L"c_{d_\mathrm{trans}}",
#        title="solidity = " * titletext * "(#$(is))",
#        ylim=(0.0, 0.5),
#    )
#    plot!(machs, cdglob*ones(N); color=myblue, label="Nominal")
#    plot!(machs, cdout; coor=myred, label="Corrected")
#    filename = "ovmach-sweepsolidity_corrcd." * lpad(is, 4, "0") * ".png"
#    savefig(savepath * filename)
#end

# stagger
##### ----- Out v Mach, Sweep stagger's ----- #####
for (is, stagger) in enumerate(staggers)
    for (i, mach) in enumerate(machs)
        clout[i], cdout[i] = dt.corrected_clcd!(
            [copy(clglob)],
            [copy(cdglob)],
            machglob,
            solidityglob,
            stagger,
            clcdminglob,
            clmaxglob,
            clminglob,
            dcldaglob,
           )
    end
    titletext = @sprintf "%1.3f" stagger
    plot(;
        xlabel="Mach Number",
        ylabel=L"c_{\ell_\mathrm{trans}}",
        title="stagger = " * titletext * "(#$(is))",
        ylim=(-2.0, 2.0),
    )
    plot!(machs, clglob*ones(N); color=myblue, label="Nominal")
    plot!(machs, clout; coor=myred, label="Corrected")
    filename = "ovmach-sweepstagger_corrcl." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)

    plot(;
        xlabel="Mach Number",
        ylabel=L"c_{d_\mathrm{trans}}",
        title="stagger = " * titletext * "(#$(is))",
        ylim=(0.0, 0.5),
    )
    plot!(machs, cdglob*ones(N); color=myblue, label="Nominal")
    plot!(machs, cdout; coor=myred, label="Corrected")
    filename = "ovmach-sweepstagger_corrcd." * lpad(is, 4, "0") * ".png"
    savefig(savepath * filename)
end
