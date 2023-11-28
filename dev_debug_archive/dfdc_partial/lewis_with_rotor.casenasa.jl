
include(
    "../../../NASA_validation/geometry_parsing/parsed_geometry_files/extracted_duct_and_centerbody_geometry_in_julia_arrays.jl",
)
hub_coordinates = hub_coordinates_closed_TE

B = 22

_, leid = findmin(duct_coordinates[:, 1])
inner_duct = reverse(duct_coordinates[1:leid, :]; dims=1)

xrotor = 0.0 # axial position of rotor stacking axis (dimensional)
Rtip = FLOWMath.akima(inner_duct[:, 1], inner_duct[:, 2], xrotor)
Rhub = FLOWMath.akima(hub_coordinates[:, 1], hub_coordinates[:, 2], xrotor)

include("../../../NASA_validation/validation_scripts/rotor_geometry.jl")
r = rotor_rct[:, 1] ./ Rtip
chords = rotor_rct[:, 2]
twists = rotor_rct[:, 3]

airfoils = []
for sidx in 1:25
    include(
        "../../../NASA_validation/validation_scripts/xrotor_airfoil_parameter_approximations/rotor_section$(sidx)_dfdc_params.jl",
    )
    push!(
        airfoils,
        dt.DFDCairfoil(
            alpha0,
            clmax,
            clmin,
            dclda,
            0.1,
            dcl_stall,
            cdmin,
            clcdmin,
            dcdcl2,
            0.0, #cmom
            re_ref,
            re_exp,
            mcrit,
        ),
    )
end

RPM = 1000.0
Omega = RPM * pi / 30  # convert from RPM to rad/s

rhoinf = 1.225
muinf = 1.81e-5
Ma = 0.05 # NASA aero paper claimed all runs were done at Ma=0.05
asound = 341.0
Vinf = Ma * asound
Vref = Vinf

#overwite npanels for now (should put conditions for this in auto generation function)
npi = 40
pref = 1
nhub_inlet = ceil(Int, pref * npi)
nduct_inlet = ceil(Int, pref * npi)
npanels = ceil.(Int, pref .* [30, 20, 40])

nwake_sheets = 11
wake_length = 1.0

rotor_parameters = [(;
    xrotor,
    r,
    chords,
    twists,
    airfoils,
    Rtip,
    Rhub,
    tip_gap=0.0,
    B,
    Omega,
    fliplift=false,
)]

paneling_constants = (; npanels, nhub_inlet, nduct_inlet, wake_length, nwake_sheets)
freestream = (; rhoinf, muinf, asound, Vinf)
reference_parameters = (; Vref=Vinf, Rref=Rtip)

