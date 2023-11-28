project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

include(project_dir * "/plots_default.jl")
include(project_dir * "/dev_debug_archive/dfdc_testing/generate_dfdc_case.jl")

using FLOWFoil

function gen_straightwake_case(;savepath="dev_debug_archive/dfdc_testing/")
    filename = "straightwake.case"

    op_data = (;
        rho=1.225, mu=1.81e-5, Vso=343.0, Vinf=20.0, Vref=20.0, Alt=0.0, RPM=5000.0
    )

    wake_data = (; nwake=20, xwake=1.0, rlx_wake="F", nwake_sheets=11)

    #this is the same.
    airfoil_data = [(;
        xisection=0.0,
        alpha0=0.0000,
        dclda=6.2800,
        clmax=1.5000,
        clmin=-1.0000,
        dclda_stall=0.50000,
        dcl_stall=0.20000,
        cmcon=0.0000,
        mcrit=0.70000,
        cdmin=0.12000E-01,
        clcdmin=0.10000,
        dcdcl2=0.50000E-02,
        Re_ref=0.20000E+06,
        Re_exp=0.35000,
    )]

    r, c, t = rotor_geom()

    xrotor = 0.127

    rotor_data = [(; naf=1, xrotor, B=5, r, chord=c, twist=t)]

    hub_coordinates = simple_hub_coordinates()

    duct_coordinates = smooth_duct_coordinates()

    plot_geom(duct_coordinates, hub_coordinates, xrotor, r, c, t)

    gen_dfdc_case(
        filename,
        op_data,
        wake_data,
        airfoil_data,
        rotor_data,
        #reverse coordinates
        reverse(hub_coordinates; dims=1),
        reverse(duct_coordinates; dims=1);
        savepath=savepath,
    )

    # - Generate Matching DuctAPE Parameters - #
    write_ducttape_params(
        filename * ".jl",
        op_data,
        wake_data,
        airfoil_data,
        rotor_data,
        hub_coordinates,
        duct_coordinates;
        savepath=savepath,
        npanels_inlet=round(Int, length(hub_coordinates[:, 1] / 2)),
    )

    return nothing
end

function smooth_duct_coordinates(; Rtip=0.127)
    raw_coords = [
        1.000000 0.000000
        0.950000 0.012800
        0.900000 0.024800
        0.800000 0.049500
        0.700000 0.071500
        0.600000 0.089800
        0.500000 0.104500
        0.400000 0.113700
        0.300000 0.119000
        0.200000 0.117800
        0.150000 0.108500
        0.100000 0.096200
        0.075000 0.088000
        0.050000 0.076200
        0.025000 0.062300
        0.012500 0.052500
        0.000000 0.029300
        0.012500 0.013000
        0.025000 0.007200
        0.050000 0.003700
        0.075000 0.001900
        0.100000 0.000000
        0.150000 0.000000
        0.200000 0.000000
        0.300000 0.000000
        0.400000 0.000000
        0.500000 0.000000
        0.600000 0.000000
        0.700000 0.000000
        0.800000 0.000000
        0.900000 0.000000
        0.950000 0.000000
        1.000000 0.000000
    ]

    smooth_coords = FLOWFoil.repanel_airfoil(raw_coords; N=160)

    xscale = 0.254 #10inches in meters
    yscale = 0.254

    smooth_coords[:, 1] .*= xscale
    smooth_coords[:, 2] .*= yscale
    smooth_coords[:, 2] .+= Rtip

    return smooth_coords
end

function simple_hub_coordinates(; scale=0.0254)

    # create ellipse for front
    a = 3.0 * scale
    b = 1.0 * scale
    th = range(pi, pi / 2.0, 46)

    xel = a * cos.(th) .+ a
    rel = b * sin.(th)

    xflat = range(xel[end], 0.254, 45)
    rflat = rel[end] * ones(45)

    x = [xel; xflat[2:end]]
    r = [rel; rflat[2:end]]
    #for some reason isn't quite zero, but needs to be for DFDC
    r[1] = 0.0

    return [x r]
end

function rotor_geom(; Rtip=0.254 / 2.0, Rhub=0.0254)

    # from DFDC
    rct = [
        0.50491E-01 0.89142E-01 69.012
        0.61567E-01 0.79785E-01 59.142
        0.72644E-01 0.71300E-01 51.825
        0.83721E-01 0.63979E-01 46.272
        0.94798E-01 0.57777E-01 41.952
        0.10587 0.52541E-01 38.509
        0.11695 0.48103E-01 35.699
        0.12803 0.44316E-01 33.354
        0.13911 0.41061E-01 31.349
        0.15018 0.38243E-01 29.596
    ]

    # scale r to hub and tip radii
    r = FLOWFoil.linear_transform([rct[1, 1]; rct[end, 1]], [Rhub; Rtip], rct[:, 1])

    # keep everything else the same
    return r, rct[:, 2], rct[:, 3]
end

function plot_geom(duct_coordinates, hub_coordinates, xrotor, r, c, t)
    plot(; xlabel="x (m)", ylabel="r (m)", aspectratio=1)

    plot!(duct_coordinates[:, 1], duct_coordinates[:, 2]; label="Duct")
    plot!(hub_coordinates[:, 1], hub_coordinates[:, 2]; label="Hub")
    plot!(xrotor * ones(length(r)), r; label="Rotor")

    savefig("dev_debug_archive/dfdc_testing/straightwake-geom.pdf")

    plot(; xlabel="r/Rtip", ylabel="chord (m)")
    plot!(r / r[end], c)
    savefig("dev_debug_archive/dfdc_testing/straightwake-chord.pdf")

    plot(; xlabel="r/Rtip", ylabel="twist (deg)")
    plot!(r / r[end], t)
    savefig("dev_debug_archive/dfdc_testing/straightwake-twist.pdf")

    return nothing
end

gen_straightwake_case()
