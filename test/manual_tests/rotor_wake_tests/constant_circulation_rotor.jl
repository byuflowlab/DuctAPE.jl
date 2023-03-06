#=

Looking at rotor only with constant circulation.

Need to use first expression to get Vm in wake since we don't know Vm_avg.

We can set a known freestream, then march that in to the hub.

There will only be vorticity at the tip and hub, where the changes in circulation take place.

=#

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

"""
"""
function gamma_m_for_completeness(Gamma_tilde, radial_stations)
    gamma_m = zeros(length(Gamma_tilde))
    for i in 1:length(radial_stations)
        if i == length(radial_stations)
            #if at tip, Gamma_tilde2 is zero since it's outside the rotor wake.
            gamma_m[i] = -(-Gamma_tilde[i]) / (2.0 * pi * radial_stations[i])
        elseif i == 1
            #if at the hub, Gamma_tilde1 is zero since it's in the hub
            gamma_m[i] = -(Gamma_tilde[i]) / (2.0 * pi * radial_stations[i])
        else
            #otherwise we have delta Gamma_tilde in the numerator
            gamma_m[i] =
                -(Gamma_tilde[i + 1] - Gamma_tilde[i]) / (2.0 * pi * radial_stations[i])
        end
    end
    return gamma_m
end

##### ----- TEST ----- #####
# let's use simple nubmers, ones all around before moving on to something else
Vinf = 1.0
B = 1.0
radial_stations = range(0.1, 1.0; step=0.1)
Gamma = ones(length(radial_stations))
Gamma_tilde = B .* Gamma
Omega = 1.0 # 2.0 * pi
H_tilde = Omega * B * Gamma / (2.0 * pi)

######
#TODO: NEEDD TO FIND A SOURCE THAT CAN TELL YOU WHAT THESE VALUES SHOULD BE!!!
#the math all works out, as in I'm getting the values I expect from the math,
#but is the original expression correct?
#
#Also need to figure out signs.  should tips have negative vorticity or positive? (see theory and think about it)
######

Vms = vm_from_vinf(Vinf, Gamma_tilde, H_tilde, radial_stations)
gamma_thetas = gamma_theta_open_rotor(Vinf, Vms)
gamma_ms = gamma_m_for_completeness(Gamma_tilde, radial_stations)

println()
println("NOTE: The following vectors are from Hub to Tip\n")
println("Vms should all be the same")
display(Vms)

println()
println(
    "gamma_thetas should be zero except for ends.  Ends should be equal and opposite if assuming Vinf at both tip and hub",
)
display(gamma_thetas)

println()
println(
    "gamma_ms should also be zero except for ends, but these aren't ever used in practice., also, the hub should be 1/rhub times the magnitude with opposite sign.",
)
display(gamma_ms)

