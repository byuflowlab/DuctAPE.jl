include("1DModel_A.jl")

# Base Values
Vinf_base = 10.0 #m/s
thrust_base = 3.0 #N
exit_area_base = pi * (2.75 * 0.0254)^2 #m
Rtip_base = 2.5 * 0.0254 #m
Rhub_base = 0.25 * Rtip_base #m
radii_base = collect(range(Rhub_base, Rtip_base, 10)) #m
num_blades_base = 2
lift_polars_base = [[5.0 * pi/180.0 0.7; 6.0 * pi/180.0 0.9]]

chords, twists, debug = opt_prelim(
    Vinf_base,
    thrust_base,
    exit_area_base,
    Rtip_base,
    Rhub_base,
    radii_base,
    num_blades_base,
    lift_polars_base;
    rho=1.225, #kg/m3
    flow_coeff=0.4,
    ambient_static_pressure=101325.0, #Pa
    ambient_static_temperature=288.15, #K
    adiabatic_stage_efficiency=0.85,
    c_p=1005.0, #J/kg-K #TODO: putting this in J instead of kJ seems to get correct units
    lift_coefficients=[0.8],
    verbose=true,
)

# using Plots
# function plotstuff(p, x, y, lab, labprefix, units="")
#     return plot!(p, x, y; label=labprefix * "$lab" * units)
# end

# # - thrust Sweeps - #
# thrust_sweep = range(1.0, 5.0, 5)
# labprefix = "Thrust = "
# units = " N"
# pv = plot(; xlabel="chords", ylabel="radial positions")

# for var in thrust_sweep
#     chords, twists, debug = opt_prelim(
#         Vinf_base,
#         var,
#         exit_area_base,
#         Rtip_base,
#         Rhub_base,
#         radii_base,
#         num_blades_base,
#         lift_polars_base;
#         rho=1.225,
#         flow_coeff=0.4,
#         ambient_static_pressure=101325.0,
#         ambient_static_temperature=288.15,
#         adiabatic_stage_efficiency=0.85,
#         c_p=1.005,
#         lift_coefficients=[0.8],
#         verbose=true,
#     )

#     plotstuff(pv, chords, radii_base, var, labprefix, units)
# end
# savefig(pv, "chord-vs-thrust.pdf")
