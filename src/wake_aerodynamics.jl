"""
    calculate_enthalpy_jumps(Gamr, Ωr, num_blades)

Calculate enthalpy jump across each blade.
"""
function calculate_enthalpy_jumps(Gamr, Ωr, num_blades)

    # number of blade elements (or radial positions) and rotors
    nr, nrotor = size(Gamr)

    # compute cumulative sum of enthalpy jumps
    H_tilde = similar(Gamr) .= 0
    for ir in 1:nr
        for irotor in 1:nrotor
            Ω = Ωr[irotor]
            B = num_blades[irotor]

            Γ = Gamr[ir, irotor]
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
    calculate_net_circulation(Gamr, num_blades)

Calculate net circulation from upstream rotors.
"""
function calculate_net_circulation(Gamr, num_blades)

    # number of blade elements (or radial positions) and rotors
    nr, nrotor = size(Gamr)

    # calculate net circulations
    Γ_tilde = similar(Gamr) .= 0
    for ir in 1:nr
        for irotor in 1:nrotor
            B = num_blades[irotor]

            Γ = Gamr[ir, irotor]
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
    calculate_wake_vortex_strengths!(gamw, Rr_wake, Vmr, Γ_tilde, H_tilde)

Calculate wake vortex strengths

`Rr_wake::Matrix{Float}` : radial locations of the wake sheets at each rotor plane. (synonymous with the rotor source panel edge locations)
`Vmr::Matrix{Float}` : meridional velocity at each radial position on the rotor planes (the meridional velocity at the rotor source panel centers, i.e., at the blade element radial stations)
"""
function calculate_wake_vortex_strengths!(gamw, Rr_wake, Vmr, Γ_tilde, H_tilde)

    # number of streamwise and radial panels
    nr, nrotor = size(gamw)

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
                gamw[ir, irotor] = 0.0
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
                gamw[ir, irotor] = (K * ΔΓ2 + ΔH) / Vm_avg
            end
        end
    end

    return gamw
end

"""
Fill in values for wake strengths behind rotors
"""
function fill_out_wake_strengths(gamw, rotor_indices, num_wake_x_panels)

    # do things the easy way if there's only one rotor.
    if length(rotor_indices) == 1
        return repeat(gamw; inner=(1, num_wake_x_panels))
    else
        # Initialize Output
        TF = eltype(gamw)
        wake_vortex_strengths = zeros(TF, length(gamw[:, 1]), num_wake_x_panels)

        # set up ranges
        # note that the plus one here and the minus one in the loop account for the difference in panel edge (rotor index) vs panel center (strength location) indexing
        ridx = [rotor_indices; num_wake_x_panels + 1]

        # loop through rotors
        for i in 1:length(rotor_indices)
            wake_vortex_strengths[:, ridx[i]:(ridx[i + 1] - 1)] .= repeat(
                gamw[:, i]; inner=(1, ridx[i + 1] - ridx[i])
            )
        end

        return wake_vortex_strengths
    end
end

function initialize_wake_vortex_strengths(Vinf, Gamr, Omega, B, rotor_panel_edges)

    # get net circulations
    Gamma_tilde = calculate_net_circulation(Gamr, B)

    # get enthalpy jumps
    H_tilde = calculate_enthalpy_jumps(Gamr, Omega, B)

    # initialize vm with freestream and rotation assuming open rotor
    vm = marched_vm(Vinf, Gamma_tilde, H_tilde, rotor_panel_edges)

    # return initialized wake vortex strengths
    gwhub = vm[1, :] .- Vinf
    gw = vm[2:end, :] .- vm[1:(end - 1), :]
    gwduct = Vinf .- vm[end, :]
    return [gwhub'; gw; gwduct']
end

"""
get meridional velocities using the marching method for when Vinf outside of wake is known.
"""
function marched_vm(Vinf, Gamma_tilde, H_tilde, rotor_panel_edges)

    #rename for convenience
    nr, nrotor = size(Gamma_tilde)

    #initialize
    Vm = zeros(nr, nrotor)

    # Loop through rotors
    for irotor in 1:nrotor

        # Loop through radial stations
        #march from just outside the tip to the hub
        for ir in nr:-1:1
            if ir == nr

                #this is the change at the tip, so we need to set Vm2 to Vinf, and set Gamma dn H 2's to zeros
                radical =
                    Vinf^2 +
                    (1.0 / (2.0 * pi * rotor_panel_edges[ir, irotor]))^2 *
                    (-Gamma_tilde[ir, irotor]^2) - 2.0 * (-H_tilde[ir, irotor])
            else

                #otherwise, we just take the differences inside as-is
                radical =
                    Vm[ir + 1, irotor]^2 +
                    (1.0 / (2.0 * pi * rotor_panel_edges[ir, irotor]))^2 *
                    (Gamma_tilde[ir + 1, irotor]^2 - Gamma_tilde[ir, irotor]^2) -
                    2.0 * (H_tilde[ir + 1, irotor] - H_tilde[ir, irotor])
            end

            # compute the meridional velocity value, avoiding negative square roots and divisions by zero later
            if radical > 0.0
                Vm[ir, irotor] = sqrt(radical)
            else
                Vm[ir, irotor] = 0.1 * Vinf
            end
        end
    end

    #Vm should now be in the order of hub to tip
    return Vm
end
