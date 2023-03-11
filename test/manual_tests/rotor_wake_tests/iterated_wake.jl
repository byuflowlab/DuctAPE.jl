#=

Look at rotor-wake model in an iterated scheme rather than a prescribed circulation scheme

=#
include("../../../plots_default.jl")
using CCBlade
using DuctTAPE
const dt = DuctTAPE
using FLOWFoil
const ff = FLOWFoil

######################################################################
#                                                                    #
#                        CCBLADE Example stuff                       #
#                                                                    #
######################################################################

Rtip = 10 / 2.0 * 0.0254  # inches to meters
Rhub = 0.10 * Rtip
B = 2  # number of blades

rotor = Rotor(Rhub, Rtip, B; tip=nothing)

propgeom = [
    0.15 0.130 32.76
    0.20 0.149 37.19
    0.25 0.173 33.54
    0.30 0.189 29.25
    0.35 0.197 25.64
    0.40 0.201 22.54
    0.45 0.200 20.27
    0.50 0.194 18.46
    0.55 0.186 17.05
    0.60 0.174 15.97
    0.65 0.160 14.87
    0.70 0.145 14.09
    0.75 0.128 13.39
    0.80 0.112 12.84
    0.85 0.096 12.25
    0.90 0.081 11.37
    0.95 0.061 10.19
    0.999 0.041 8.99
]

r = propgeom[:, 1] * Rtip
chord = propgeom[:, 2] * Rtip
theta = propgeom[:, 3] * pi / 180

af = AlphaAF("test/data/naca4412.dat")

sections = Section.(r, chord, theta, Ref(af))

Vinf = 5.0
Omega = 5400 * pi / 30  # convert to rad/s
rho = 1.225

op = simple_op.(Vinf, Omega, r, rho)

out = CCBlade.solve.(Ref(rotor), sections, op)

## -- Get Circulation from outputs -- ##
function get_gamma(ccbouts, chord)
    cl = ccbouts.cl
    W = ccbouts.W

    return 0.5 .* W .* cl .* chord
end

ccbGamma = get_gamma(out, chord)

######################################################################
#                                                                    #
#                        DuctTAPE Version                            #
#                                                                    #
######################################################################

"""
here Gamma_tilde will just be B*Gamma, since we are only looking at a single rotor.
"""
function vm_from_vinf(Vinf, Gamma_tilde, H_tilde, radial_stations)

    #rename for convenience
    nr = length(radial_stations)

    #initialize
    Vm = zeros(nr)

    #march from just outside the tip to the hub
    for i in nr:-1:1
        if i == nr

            #this is the change at the tip, so we need to set Vm2 to Vinf, and set Gamma dn H 2's to zeros
            radical =
                Vinf^2 + (1.0 / (2.0 * pi * radial_stations[i]))^2 * (-Gamma_tilde[i]^2) -
                2.0 * (-H_tilde[i])
            # if radical > 0.0
            Vm[i] = sqrt(radical)
            # else
            #     Vm[i] = 0.0
            # end
        else

            #otherwise, we just take the differences inside as-is
            radical =
                Vm[i + 1]^2 +
                (1.0 / (2.0 * pi * radial_stations[i]))^2 *
                (Gamma_tilde[i + 1]^2 - Gamma_tilde[i]^2) -
                2.0 * (H_tilde[i + 1] - H_tilde[i])
            # if radical > 0.0
            Vm[i] = sqrt(radical)
            # else
            #     Vm[i] = 0.0
            # end
        end
    end

    #Vm should now be in the order of hub to tip
    return Vm
end

"""
"""
function gamma_theta_open_rotor(Vinf, Vm)
    gamma_theta = zeros(length(Vm))

    for i in 1:length(Vm)
        if i == 1
            #this is at the hub, where the internal velocity is zero,
            #so we should see a non-zero vorticity here
            gamma_theta[i] = Vm[i] - Vinf #assumes hollow rotor? 0.0
        elseif i == length(Vm)
            #this is the rotor tip, where we again have a difference from the freestream
            #so we should also see a non-zero vorticity here
            gamma_theta[i] = Vinf - Vm[i]
        else
            # this is all the other cases, where circulation was constant, so the velocities should all be the same and we should see zero circulation
            gamma_theta[i] = Vm[i + 1] - Vm[i]
        end
    end
    return gamma_theta
end

#---------------------------------#
#     Generate Wake Geometry      #
#---------------------------------#

nbe = length(r)
# Just make it a rectangle that is 2D long.
D = 2.0 * Rtip
#make the grid square using the same spacing in x as was done in r
xrange = range(0.0, 2.0 * D; step=0.05 * Rtip)

grid_points = [[x r] for x in xrange, r in r]

# grid ponts for making panels are in x rows, and r columns matrices
# all the x's are in one matrix, and all the r's in another.
x_grid_points = repeat(xrange; outer=(1, nbe))
r_grid_points = repeat(r'; outer=(length(xrange), 1))

#returns the wake panels as vectors for each "row"
wake_panels = dt.generate_wake_panels(
    x_grid_points,
    r_grid_points,
    nbe;
    method=ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true]),
)

# Rotor points of interest
x_sample_points = zeros(nbe + 1)
r_sample_points = zeros(nbe + 1)
r_sample_points[2:(end - 1)] = (r[1:(end - 1)] .+ r[2:end]) / 2.0

r_sample_points[1] = 2.0 * r[1] - r_sample_points[2]
r_sample_points[end] = 2.0 * r[end] - r_sample_points[end - 1]

# use body of revolution method so there aren't any rogue kutta conditions used later
method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true])

# - Generate panels to sample velocity at rotor blade element locations - #
rotor_affect_panels = ff.generate_panels(method, [x_sample_points r_sample_points])

## - PLOT GEOMETRY -- ##
# - Plot grid to max sure you made it correctly - #
plot(
    getindex.(grid_points, 1),
    getindex.(grid_points, 2);
    color=mycolors[1],
    aspectratio=:equal,
    seriestype=:scatter,
    markersize=1,
    xlabel="x",
    ylabel="r",
    label="",
)

# - Plot the panels to make sure you made them correctly - #
# plot(; aspectratio=:equal)
for i in 1:nbe
    plot!(
        wake_panels[i].panel_center[:, 1],
        wake_panels[i].panel_center[:, 2];
        color=mycolors[2],
        seriestype=scatter,
        markersize=1,
        label="",
    )
end

# - Plot affect panel centers to make sure they're in the right place - #
plot!(
    rotor_affect_panels.panel_center[:, 1],
    rotor_affect_panels.panel_center[:, 2];
    color=mycolors[3],
    seriestype=scatter,
    markersize=1,
    label="",
)
savefig("./test/manual_tests/rotor_wake_tests/iter_wake_geometry.pdf")

## -- Generate the geometry meshes for finding the velocity at all the points based on the wake panels -- ##
wake_to_rotor_mesh = [
    dt.generate_one_way_mesh(wake_panels[i], rotor_affect_panels) for
    i in 1:length(wake_panels)
]

## -- Generate the unit induced velocities -- ##
#TODO: I think the assemble_induced_velocity_matrices needs a new unit test after the update.
A_wake_to_rotor = [
    dt.assemble_induced_velocity_matrices(
        wake_to_rotor_mesh[i], wake_panels[i], rotor_affect_panels
    ) for i in 1:length(wake_panels)
]
vxd_wake_to_rotor = [A_wake_to_rotor[i][1] for i in 1:length(wake_panels)]
vrd_wake_to_rotor = [A_wake_to_rotor[i][2] for i in 1:length(wake_panels)]

#---------------------------------#
#            INITIALIZE           #
#---------------------------------#

# - Define some basic blade element named tuples as place holders
blade_elements = [(omega=Omega, radial_positions=r, chords=chord)]

# initialize induced velocities at rotors to zeros
vm_init = zeros(nbe)
vtheta_init = zeros(nbe)

# calculate inflow from freestream and rotation
inflow_init = dt.calculate_inflow_velocities(
    blade_elements, Ref(Vinf), vm_init, vtheta_init
)

# calculate angle of attack
alpha_init = dt.calculate_angle_of_attack(theta, inflow_init.Wm, inflow_init.Wtheta)

# - Look up lift and drag data - #
cl_init = zeros(nbe)
cd_init = zeros(nbe)
for a in 1:nbe
    cl_init[a], cd_init[a] = dt.search_polars(af, alpha_init[a])
end

# - Calculate Circulation from lift and inflow_init - #
Gamma = 0.5 .* inflow_init.Wmag .* chord .* cl_init
# Gamma = deepcopy(ccbGamma)

# get net circulation and B*Gamma
BGamma_init, Gamma_tilde_init = dt.calculate_net_circulation(Gamma, B)

# - Calculate Enthalpy Jumps - #
H_tilde_init = dt.calculate_enthalpy_jumps(Gamma, Omega, B)

# - Calculate Vm in wake - #
Vm_init = vm_from_vinf(Vinf, Gamma_tilde_init, H_tilde_init, r)

# - Calculate gamma_thetas - #
wake_vortex_strengths = gamma_theta_open_rotor(Vinf, Vm_init)

# - Populate Wake velocity and strength matrices - #
# wake_Vms = repeat(Vm_init; inner=(1, length(xrange) - 1))
wake_gammas = repeat(wake_vortex_strengths; inner=(1, length(xrange) - 1))

## -- Set up some plots for iteration -- ##
pcirc = plot(Gamma, r; xlabel=L"\Gamma = 0.5 W c c_\ell", ylabel="r", label="Init")
palpha = plot(
    alpha_init * 180.0 / pi, r; xlabel=L"\alpha (degrees)", ylabel="r", label="Init"
)
pvm = plot(
    Vm_init .- Vinf,
    r;
    xlabel="farfield induced axial (meridional) velocity",
    ylabel="r",
    label="Init",
)
pvt = plot(
    vtheta_init, r; xlabel="nearfield induced tangential velocity", ylabel="r", label="Init"
)
pw = plot(inflow_init.Wmag, r; xlabel="inflow velocity magnitude", ylabel="r", label="Init")
pgammatheta = plot(
    wake_vortex_strengths, r; xlabel="wake vortex strengths", ylabel="r", label="Init"
)
pvmwake = plot(; xlabel="wake-induced meridional velocity at rotor plane", ylabel="r")

plot!(pcirc, ccbGamma, r; label="CCBlade")
plot!(palpha, out.alpha * 180.0 / pi, r; label="CCBlade")
plot!(pvm, out.u * 2.0, r; label="CCBlade")
plot!(pvt, out.v, r; label="CCBlade")
plot!(pw, out.W, r; label="CCBlade")

#---------------------------------#
#             Iterate             #
#---------------------------------#
Gamma_temp = 99 * ones(length(Gamma))
iter = [0]
Vm = deepcopy(Vm_init)
while abs(Gamma_temp[5] - Gamma[5]) > 1e-3
    iter[1] += 1
    println("iter $(iter[1])")

    # - Calculate induced velocities at rotor - #

    # vm comes from wake affect on rotor
    vm_wake_on_rotor = zeros(nbe)
    for i in length(wake_to_rotor_mesh)
        vx = vxd_wake_to_rotor[i] * wake_gammas[i, :]
        vr = vrd_wake_to_rotor[i] * wake_gammas[i, :]
        vm_wake_on_rotor .+= sqrt.(vx .^ 2 .+ vr .^ 2)
    end

    # vtheta comes from rotor self-induction
    vtheta = 1.0 ./ (2.0 .* pi .* r) .* (0.5 .* B .* Gamma)

    # calculate inflow from freestream and rotation
    inflow = dt.calculate_inflow_velocities(
        blade_elements, Ref(Vinf), vm_wake_on_rotor, vtheta
    )

    # calculate angle of attack
    alpha = dt.calculate_angle_of_attack(theta, inflow.Wm, inflow.Wtheta)

    # - Look up lfit and drag data - #
    cl = zeros(nbe)
    cd = zeros(nbe)
    for a in 1:length(alpha)
        cl[a], cd[a] = dt.search_polars(af, alpha[a])
    end

    # - swap old to new - #
    Gamma_temp .= Gamma

    # - Calculate Circulation from lift and inflow - #
    @. Gamma = 0.5 * inflow.Wmag * chord * cl

    # get net circulation and B*Gamma
    BGamma = Gamma * B
    Gamma_tilde = B * Gamma #only because we only have one rotor

    # - Calculate Enthalpy Jumps - #
    H_tilde = dt.calculate_enthalpy_jumps(Gamma, Omega, B)

    # - Calculate Vm in wake - #
    Vm .= vm_from_vinf(Vinf, Gamma_tilde, H_tilde, r)

    # - Calculate gamma_thetas - #
    wake_vortex_strengths .= gamma_theta_open_rotor(Vinf, Vm)

    # - Populate Wake velocity and strength matrices - #
    # wake_Vms .= repeat(Vm; inner=(1, length(xrange) - 1))
    wake_gammas .= repeat(wake_vortex_strengths; inner=(1, length(xrange) - 1))

    # - PLOT - #
    plot!(pcirc, Gamma, r; label="iter #$(iter[1])")
    plot!(palpha, alpha * 180.0 / pi, r; label="iter #$(iter[1])")
    plot!(pvm, Vm .- Vinf, r; label="iter #$(iter[1])")
    plot!(pvt, vtheta, r; label="iter #$(iter[1])")
    plot!(pw, inflow.Wmag, r; label="iter #$(iter[1])")
    plot!(pgammatheta, wake_vortex_strengths, r; label="iter #$(iter[1])")
    plot!(pvmwake, vm_wake_on_rotor, r; label="iter #$(iter[1])")
end

savefig(pcirc, "test/manual_tests/rotor_wake_tests/circulation_iteration.pdf")
savefig(palpha, "test/manual_tests/rotor_wake_tests/alpha_iteration.pdf")
savefig(pvm, "test/manual_tests/rotor_wake_tests/vm_iteration.pdf")
savefig(pvt, "test/manual_tests/rotor_wake_tests/vtheta_iteration.pdf")
savefig(pw, "test/manual_tests/rotor_wake_tests/Wmag_iteration.pdf")
savefig(
    pgammatheta, "test/manual_tests/rotor_wake_tests/wake_vortex_strength_iteration.pdf"
)
savefig(pvmwake, "test/manual_tests/rotor_wake_tests/vm_wake_on_rotor_iteration.pdf")

## -- Plots -- ##
plot(ccbGamma, r; xlabel=L"\Gamma", ylabel="r", label="CCBlade")
plot!(Gamma, r; label="DuctTAPE")
savefig("test/manual_tests/rotor_wake_tests/iterwake_circulation.pdf")

plot(out.u, r; xlabel="induced axial velocity", ylabel="r", label="CCBlade")
plot!((Vm .- Vinf), r; label="DuctTAPE")
savefig("test/manual_tests/rotor_wake_tests/iterwake_Vx.pdf")

# # f = open("test/manual_tests/rotor_wake_tests/save_gamma_ccb_init.jl", "w+")
# f = open("test/manual_tests/rotor_wake_tests/save_gamma_my_init.jl", "a+")

# # write(f, "myinitgam = [\n")
# write(f, "ccbinitgam = [\n")

# for i in 1:length(Gamma)
#     write(f, "$(Gamma[i])\n")
# end
# write(f, "]\n")
# close(f)
