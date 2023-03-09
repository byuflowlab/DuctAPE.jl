#=

Compare rotor/wake model with ccblade to see if you get anything remotely close

=#
include("../../../plots_default.jl")

#### ---- WAKE MODEL ---- ####
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
            Vm[i] = sqrt(
                Vinf^2 + (1.0 / (2.0 * pi * radial_stations[i]))^2 * (-Gamma_tilde[i]^2) -
                2.0 * (-H_tilde[i]),
            )
        else

            #otherwise, we just take the differences inside as-is
            Vm[i] = sqrt(
                Vm[i + 1]^2 +
                (1.0 / (2.0 * pi * radial_stations[i]))^2 *
                (Gamma_tilde[i + 1]^2 - Gamma_tilde[i]^2) -
                2.0 * (H_tilde[i + 1] - H_tilde[i]),
            )
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

######################################################################
#                                                                    #
#                        CCBLADE Example stuff                       #
#                                                                    #
######################################################################
using CCBlade

Rtip = 10 / 2.0 * 0.0254  # inches to meters
Rhub = 0.15 * Rtip
B = 2  # number of blades

rotor = Rotor(Rhub, Rtip, B; tip=nothing)

propgeom = [
    0.1501 0.130 32.76
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

out = solve.(Ref(rotor), sections, op)

## -- Get Circulation from outputs -- ##
function get_gamma(ccbouts, chord)
    cl = ccbouts.cl
    W = ccbouts.W

    return 0.5 .* W .* cl .* chord
end

Gamma = get_gamma(out, chord)

# - Calculate Vms and gamma_thetas - #
Gamma_tilde = B .* Gamma
H_tilde = Omega * B * Gamma / (2.0 * pi)
Vms = vm_from_vinf(Vinf, Gamma_tilde, H_tilde, r)
gamma_thetas = gamma_theta_open_rotor(Vinf, Vms)

## -- Plots -- ##
#Circulation
plot(Gamma, r; xlabel=L"\Gamma", ylabel="r")
savefig("test/manual_tests/rotor_wake_tests/ccblade_circulation.pdf")

# axial induced velocity
plot(out.u, r; xlabel="induced axial velocity", ylabel="r", label="CCBlade")
plot!((Vms .- Vinf), r; label="DuctTAPE")
savefig("test/manual_tests/rotor_wake_tests/ccblade_Vx.pdf")

# tangential induced velocity
plot(out.v, r; xlabel="induced tangential velocity", ylabel="r", label="CCBlade")
plot!(Gamma * B ./ (2.0 * pi * r), r; label="DuctTAPE")
savefig("test/manual_tests/rotor_wake_tests/ccblade_Vtheta.pdf")

# plot(gamma_thetas, r; xlabel=L"\gamma_\theta", ylabel="r")
# savefig("test/manual_tests/rotor_wake_tests/gamma_theta_fyi.pdf")

######################################################################
#                                                                    #
#                        varying Vx inputs                           #
#                                                                    #
######################################################################

Vy = Omega * r
Vx = Vinf * range(0.1, 1.0; length=length(r))
rho = 1.225
mu = 1.81e-5
pitch = 0.0
asound = 341.0

opvar = CCBlade.OperatingPoint.(collect(Vx), Vy, rho, pitch, mu, asound)

out = solve.(Ref(rotor), sections, opvar)

Gamma = get_gamma(out, chord)

## -- Plots -- ##
plot(Gamma, r; xlabel=L"\Gamma", ylabel="r")
savefig("test/manual_tests/rotor_wake_tests/ccblade_circulation_vxvar.pdf")

Gamma_tilde = B .* Gamma
H_tilde = Omega * B * Gamma / (2.0 * pi)
Vms = vm_from_vinf(Vinf, Gamma_tilde, H_tilde, r)
gamma_thetas = gamma_theta_open_rotor(Vinf, Vms)

plot(out.u, r; xlabel="induced axial velocity", ylabel="r", label="CCBlade")
plot!((Vms .- Vx), r; label="DuctTAPE")
savefig("test/manual_tests/rotor_wake_tests/ccblade_Vxvar.pdf")

plot(out.v, r; xlabel="induced tangential velocity", ylabel="r", label="CCBlade")
plot!(Gamma * B ./ (4.0 * pi * r), r; label="DuctTAPE")
savefig("test/manual_tests/rotor_wake_tests/ccblade_Vthetavar.pdf")
