"""
    calculate_enthalpy_jumps(Γr, Ωr, num_blades)

Calculate enthalpy jump across each blade.
"""
function calculate_enthalpy_jumps(Γr, Ωr, num_blades)

    # number of blade elements (or radial positions) and rotors
    nr, nrotor = size(Γr)

    # compute cumulative sum of enthalpy jumps
    H_tilde = similar(Γr) .= 0
    for ir in 1:nr
        for irotor in 1:nrotor
            Ω = Ωr[irotor]
            B = num_blades[irotor]

            Γ = Γr[ir, irotor]
            H_tilde[ir, irotor] += Ω * B * Γ / (2.0 * pi)

            # add the upstream contributions
            for jrotor in 1:(irotor - 1)
                H_tilde[ir, irotor] += H_tilde[ir, jrotor]
            end
        end
    end

    return H_tilde
end

"""
    calculate_net_circulation(Γr, num_blades)

Calculate net circulation from upstream rotors.
"""
function calculate_net_circulation(Γr, num_blades)

    # number of blade elements (or radial positions) and rotors
    nr, nrotor = size(Γr)

    # calculate net circulations
    Γ_tilde = similar(Γr) .= 0
    for ir in 1:nr
        for irotor in 1:nrotor
            B = num_blades[irotor]

            Γ = Γr[ir, irotor]
            Γ_tilde[ir, irotor] += B * Γ

            # add the upstream contributions
            for jrotor in 1:(irotor - 1)
                Γ_tilde[ir, irotor] += Γ_tilde[ir, jrotor]
            end
        end
    end

    return Γ_tilde
end

"""
    calculate_wake_vortex_strengths!(Γw, Rr_wake, Vmr, Γ_tilde, H_tilde)

Calculate wake vortex strengths

`Rr_wake::Matrix{Float}` : radial locations of the wake sheets at each rotor plane. (synonymous with the rotor source panel edge locations)
`Vmr::Matrix{Float}` : meridional velocity at each radial position on the rotor planes (the meridional velocity at the rotor source panel centers, i.e., at the blade element radial stations)
"""
function calculate_wake_vortex_strengths!(Γw, Rr_wake, Vmr, Γ_tilde, H_tilde)

    # number of streamwise and radial panels
    nr, nrotor = size(Γw)

    #loop  through rotors where wake strengths are generated
    for irotor in 1:nrotor

        # loop through each radial position
        for ir in 1:nr

            # compute average meridional velocities
            if ir == 1
                # use panel velocity as average velocity
                Vm_avg = Vmr[ir, irotor]
            elseif ir == nr
                # use panel velocity as average velocity
                Vm_avg = Vmr[ir - 1, irotor]
            else
                # average velocities above and below the current panel
                Vm_avg = 1 / 2 * (Vmr[ir - 1, irotor] + Vmr[ir, irotor])
            end

            # calculate the wake vortex strength
            if Vm_avg <= 0.0
                # avoid division by zero
                Γw[ir, irotor] = 0.0
            else

                # radial positions of wake sheets at rotor plane
                r = Rr_wake[ir, irotor]

                # constant in expression for convenience
                K = -1 / (8 * pi^2 * r^2)

                # net circulation squared and enthalpy jumps across sheet
                if ir == 1
                    ΔΓ2 = Γ_tilde[ir, irotor]^2
                    ΔH = H_tilde[ir, irotor]
                elseif ir == nr
                    ΔΓ2 = -Γ_tilde[ir - 1, irotor]^2
                    ΔH = -H_tilde[ir - 1, irotor]
                else
                    ΔΓ2 = Γ_tilde[ir, irotor]^2 - Γ_tilde[ir - 1, irotor]^2
                    ΔH = H_tilde[ir, irotor] - H_tilde[ir - 1, irotor]
                end

                # wake strength density taken from rotor to next rotor constant along streamlines
                Γw[ir, irotor] = (K * ΔΓ2 + ΔH) / Vm_avg
            end
        end
    end

    return Γw
end

#TODO: probably delete this, unless it ends up being needed
# function initialize_wake_vortex_strengths(Γr, Ωr, freestream)

#     # set initial Vm values
#     Vm = similar(Γr, nr - 1, nx) .= freestream.Vinf

#     # calculate enthalpy jumps
#     H_tilde = calculate_enthalpy_jumps(Γr, Ωr, B)

#     # calculate net circulation
#     Γ_tilde = calculate_net_circulation(Γr, B)

#     # calculate wake vortex strengths
#     calculate_wake_vortex_strengths!(Γw, wake_panels, Vm, Γ_tilde, H_tilde, rotor_indices)

#     # println(wake_vortex_strengths)
#     return wake_vortex_strengths
# end
