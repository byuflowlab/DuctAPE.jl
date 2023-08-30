using Xfoil
using LsqFit
using PrettyTables

include("../plots_default_new.jl")
# include("../geometry_parsing/smooth_af_data_approximator.jl")
include("xrotor_af.jl")
include("rotor_geometry.jl")

savepath = "xrotor_airfoil_parameter_approximations/"

# TODO: write a function to run the xfoil data at the reynolds numbers estimated from the data based on radius and chord.
# then take that data, and the parameters run and estimate the dfdc parameters needed.
# then to double check you can go grab the xrotor function from acoustic analogies and plot at the same reynolds number the range of angles of attack to see if things match well.
# potentially use the xrotor function as a model function and run a least squares fit on the xfoil outputs, with the things known as knowns (e.g. clmax, clmin, etc.), just fitting the lift slope and various factors? can just do all the factors, or will need a function inside the lsq wrapper that gets the various known quantities

#---------------------------------#
#   Global Operating Parameters   #
#---------------------------------#
RPM_rotor = 7808.0 #mechanical RPM in hotwire data readme
Omega_rotor = RPM_rotor * pi / 30.0
r = rotor_rct[:, 1] #from rotor_chord_stagger.jl
chord = rotor_rct[:, 2] #from rotor_chord_stagger.jl
Ma = 0.05 # NASA aero paper claimed all runs were done at Ma=0.05
asound = 341.0
Vinf = Ma * asound
rhoinf = 1.225 # choose values that match asound
muinf = 1.81e-5 # choose values that match asound
Vmags = sqrt.(Vinf^2 .+ (Omega_rotor .* r) .^ 2)
Res = rhoinf * Vmags .* chord / muinf

#---------------------------------#
#            Run XFOIL            #
#---------------------------------#
function run_xfoil(x, y)
    # Angle of attack ranges
    minaoa = -8.0
    maxaoa = 15.0
    da = 0.1
    lenaoa = ceil(Int, (maxaoa - minaoa) / da + 1)
    aoa = range(minaoa, maxaoa, lenaoa)

    cl, cdrag, cdp, cm, conv = Xfoil.alpha_sweep(
        x,
        y,
        aoa,
        1e6;#use an re that will converge, but pass the lower Re in for the xrotor drag stuff
        mach=0.0, #dfdc model will take care of mach corrections.
        reinit=false,
        percussive_maintenance=true,
        printdata=true,
        # clmaxstop=true,
        # clminstop=true,
    )

    alphaconv = []
    clconv = []
    cdconv = []
    for (ic, c) in enumerate(conv)
        if c == 1
            push!(alphaconv, aoa[ic])
            push!(clconv, cl[ic])
            push!(cdconv, cdrag[ic])
        end
    end

    alphacut, clcut, cdcut = poststallcut(alphaconv, clconv, cdconv)

    plot(alphaconv, clconv; label="conv")
    plot!(alphacut, clcut; label="cut")
    savefig(savepath * "clcut.pdf")

    return alphacut * pi / 180.0, clcut, cdcut

    # return alphaconv * pi / 180.0, clconv, cdconv
end

function poststallcut(a, cl, cdrag)
    stallcount = [0]
    poststallcount = [0]
    cutoff = [1]
    for id in 2:length(cl)
        if a[id] > 2.0
            if cl[id] < cl[id - 1]
                stallcount[1] += 1
            end
            if stallcount[1] > 0 && cl[id] > cl[id - 1]
                poststallcount[1] += 1
            end
        end
        if poststallcount[1] > 4
            cutoff[1] = id - poststallcount[1]
            return a[1:cutoff[1]], cl[1:cutoff[1]], cdrag[1:cutoff[1]]
        end
    end

    return a, cl, cdrag
end

# for isec in 1:length(r)
# for isec in 1:1
# for isec in 2:2
# for isec in 3:3
for isec in (3, 23, 24)
    # for isec in 23:23
    # for isec in (1, 5, 10, 15, 25)
    # for isec in length(r):length(r)
    function wrapxrotor(aoa, params)

        # dimensions and types
        TF = promote_type(eltype(aoa), eltype(params))
        na = length(aoa)

        config = XROTORAirfoilConfig(
            params[1],
            params[2],
            params[3],
            params[4],
            params[5],
            params[6],
            params[7],
            params[8],
            params[9],
            params[10],
            params[11],
            params[12],
        )

        clcd_xrotor = zeros(TF, na, 2)
        for (ia, a) in enumerate(aoa)
            clcd_xrotor[ia, 1], clcd_xrotor[ia, 2] = af_xrotor(a, 1e6, 0.0, config)
        end

        return reduce(vcat, clcd_xrotor)
    end

    # load geometry
    include("../geometry_parsing/jls/rotor_section_$(isec)_m_rtheta.jl")

    # run exfoil
    alphadata, cldata, cddata = run_xfoil(mr[:, 1], mr[:, 2])

    #plot xfoil
    pcl = plot(alphadata * 180 / pi, cldata; xlabel="aoa (deg)", ylabel="cl", label="xfoil")
    pcd = plot(alphadata * 180 / pi, cddata; xlabel="aoa (deg)", ylabel="cd", label="xfoil")

    # set up known lsq parameters
    # clmax
    clmax, maxid = findmax(cldata)
    amax = alphadata[maxid]

    # clmin
    clmin, minid = findmin(cldata)
    amin = alphadata[minid]

    # zero lift angle
    #TODO: interpolate this to be more accurate
    cl0, zeroid = findmin(abs.(cldata))
    alpha0 = alphadata[zeroid]

    # cdrag min
    cdmin, minid = findmin(cddata)
    cdamin = alphadata[minid]

    # cl at cdmin
    clcdmin = cldata[minid]

    # reference reynolds
    re_ref = 1e6
    # mcrit
    mcrit = 0.7

    # initial guesses
    dclda = 2.0 * pi
    dcl_stall = 0.1
    dclda_stall = -pi / 2.0
    dcdcl2 = 0.005
    re_exp = -0.4

    # set up lsq
    p0 = [
        alpha0
        dclda
        clmax
        clmin
        dcl_stall
        dclda_stall
        cdmin
        clcdmin
        dcdcl2
        re_ref
        re_exp
        mcrit
    ]

    # lb = [
    #     # alpha0
    #     -pi
    #     pi #lift slope
    #     # clmax
    #     0.1
    #     # clmin
    #     -3.0
    #     0.0001 #dcl_stall
    #     -4 * pi #dclda_stall
    #     # cdmin
    #     0.0 #cdmin
    #     # clcdmin
    #     -3.0 #clcdmin
    #     0.0001 #dcdcd2
    #     re_ref
    #     -1.0 #re_exp
    #     mcrit
    # ]

    lb = [
        alpha0
        pi #lift slope
        # clmax
        eps() # clmax
        # clmin
        -5.0 # clmin
        0.0001 #dcl_stall
        -4 * pi #dclda_stall
        cdmin
        clcdmin
        0.0001 #dcdcd2
        re_ref
        re_exp
        mcrit
    ]

    ub = [
        alpha0
        4 * pi #lift slope
        # clmax
        5.0 # clmax
        # clmin
        5.0 # clmin
        1.0 #dcl_stall
        -pi / 4.0 #dclda_stall
        cdmin
        clcdmin
        0.1 #dcdcd2
        re_ref
        re_exp
        mcrit
    ]

    # ub = [
    #     # alpha0
    #     pi
    #     4 * pi #lift slope
    #     # clmax
    #     3.0
    #     # clmin
    #     50.0
    #     1.0 #dcl_stall
    #     # 4 * pi #dclda_stall
    #     10.0 #dclda_stall
    #     # cdmin
    #     3.0 #cdmin
    #     # clcdmin
    #     3.0 # clcdmin
    #     0.1 #dcdcd2
    #     re_ref
    #     -0.001 #re_exp
    #     mcrit
    # ]

    pretty_table([lb p0 ub]; header=["Lower", "Initial", "Upper"])

    # run lsq
    fit = LsqFit.curve_fit(
        wrapxrotor,
        alphadata,
        [cldata; cddata],
        p0;
        lower=lb,
        upper=ub,
        autodiff=:forwarddiff,
        show_trace=true,
    )

    # check
    alphacheck = range(-pi / 8, pi / 8, 91)
    clcd_fit = wrapxrotor(alphacheck, fit.param)
    clcd = reshape(clcd_fit, length(alphacheck), 2)
    plot!(pcl, alphacheck * 180 / pi, clcd[:, 1]; label="fit")
    plot!(pcd, alphacheck * 180 / pi, clcd[:, 2]; label="fit")
    savefig(pcl, savepath * "sec$(isec)_cl_fit.pdf")
    savefig(pcd, savepath * "sec$(isec)_cd_fit.pdf")

    f = open(savepath * "rotor_section$(isec)_dfdc_params.jl", "w")
    write(f, "alpha0=$(fit.param[1])\n")
    write(f, "dclda=$(fit.param[2])\n")
    write(f, "clmax=$(fit.param[3])\n")
    write(f, "clmin=$(fit.param[4])\n")
    write(f, "dcl_stall=$(fit.param[5])\n")
    write(f, "dclda_stall=$(fit.param[6])\n")
    write(f, "cdmin=$(fit.param[7])\n")
    write(f, "clcdmin=$(fit.param[8])\n")
    write(f, "dcdcl2=$(fit.param[9])\n")
    write(f, "re_ref=$(fit.param[10])\n")
    write(f, "re_exp=$(fit.param[11])\n")
    write(f, "mcrit=$(fit.param[12])\n")
    close(f)
end
