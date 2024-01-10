project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

datapath = project_dir * "/dev_debug_archive/dfdc_partial/"
savepath = datapath

include(project_dir * "/convenience_functions/generate_dfdc_case.jl")

# body geometry
scale = 1.0
include(project_dir * "/test/data/naca_662-015.jl")
include(project_dir * "/test/data/bodyofrevolutioncoords.jl")
hub_coordinates = hub_coordinates[1:end-1,:]

op_data = (; rhoinf=1.226, muinf=1.78e-5, Vso=340.0, Vinf=20.0, Vref=20.0, Alt=0.0, RPM=1000.0)

wake_data = (; nwake_sheets=11, nwake=20, xwake=1.0, rlx_wake="F\n")

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
    dcddcl2=0.50000E-02,
    Re_ref=0.20000E+06,
    Re_exp=0.35000,
)]

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

r = range(
    maximum(hub_coordinates[:, 2]),
    minimum(duct_coordinates[:, 2]);
    length=length(rct[:, 1]),
)

rotor_data = [(; naf=1, rotorzloc=0.5/scale, B=5, r, chord=rct[:, 2], twist=rct[:, 3])]

# filename = "lwr"
# gen_dfdc_case(
#     filename,
#     op_data,
#     wake_data,
#     airfoil_data,
#     rotor_data,
#     hub_coordinates./scale,
#     duct_coordinates./scale;
#     savepath=savepath,
#     version=0.70,
#     case_name="Lewis with Rotor",
# )

filename = "lewis_with_rotor.case.jl"
write_ducttape_params(
    filename,
    op_data,
    wake_data,
    airfoil_data,
    rotor_data,
    hub_coordinates,
    duct_coordinates;
    savepath=savepath,
    npanels_inlet=20,
)
