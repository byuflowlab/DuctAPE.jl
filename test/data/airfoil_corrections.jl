using CCBlade

"""
"""
function airfoil_corrections(
    alpha,
    cl,
    cd,
    cr75,
    tsr,
    Re,
    file_name="test/data/naca_4412_extrapolated_rotated_APCshifted.dat",
    file_header="NACA 4412 w/ rotation",
)

    # - Extend Data - #
    alpha_ext, cl_ext, cd_ext = viterna(alpha, cl, cd, cr75)

    # - Rotation Corrections - #
    alpha_rot = alpha_ext
    cl_rot = similar(cl_ext)
    cd_rot = similar(cd_ext)

    rR = 0.75  # r/R = 75%

    for i in 1:length(cl_ext)
        cl_rot[i], cd_rot[i] = rotation_correction(
            DuSeligEggers(), cl_ext[i], cd_ext[i], cr75, rR, tsr, alpha_ext[i]
        )
    end

    # - Write Airfoil - #
    af_final = AlphaAF(alpha_rot, cl_rot, cd_rot, file_header, Re, 0.0)
    write_af(file_name, af_final)

    return nothing
end

include("naca_4412_raw.jl")

cr75 = 0.128
tsr = 6.0
Re = 2e6
file_header = "naca_4412_extrapolated_rotated_APCshifted"
file_name = "test/data/naca_4412_extrapolated_rotated_APCshifted.dat"
airfoil_corrections(
    alpha .- (1.0 * pi / 180.0), cl, cd, cr75, tsr, Re, file_name, file_header
)
