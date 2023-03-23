"""
    calculate_enthalpy_jumps(Γr, Ωr, num_blades)

Calculate enthalpy jump across each blade.
"""
function calculate_enthalpy_jumps(Γr, Ωr, num_blades)

    # number of blade elements (or radial positions) and rotors
    nr, nrotor = size(Γr)

    # compute cumulative sum of enthalpy jumps
    H_tilde = similar(Γr) .= 0
    for irotor in 1:nrotor
        Ω = Ωr[irotor]
        B = num_blades[irotor]
        for ir = 1:nr
            Γ = Γr[ir, irotor]
            H_tilde[ir, irotor] += Ω*B*Γ/(2.0*pi)
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
    for irotor in 1:nrotor
        B = num_blades[irotor]
        for ir in 1:nr
            Γ = Γr[ir, irotor]
            Γ_tilde[ir, irotor] += B*Γ
        end
    end

    return Γ_tilde
end

"""
    calculate_wake_velocities(Ax_bw, Ar_bw, Γb, Ax_ww, Ar_ww, Γw)

Calculate the magnitude of the meridional velocity.  This velocity is defined tangent to a
streamline in the x-r plane such that its magnitude is given by: `Vm = Vx + Vr`
"""
function calculate_wake_velocities(Ax_bw, Ar_bw, Γb, Ax_ww, Ar_ww, Γw, Ax_rw, Ar_rw, Σr)

    # number of streamwise and radial panels
    nx, nr = size(Γw)

    # initialize meridional velocity magnitude
    Vm = similar(Γw) .= 0

    # loop through each radial set of wake panels
    for ir in 1:nr

        # add body induced velocities
        Vm[:,ir] .+= Ax_bw[ir] * Γb
        Vm[:,ir] .+= Ar_bw[ir] * Γb

        # add wake induced velocities
        for iwake in 1:size(Γw, 1)
            Vm[:, ir] .+= Ax_ww[iwake, ir] * Γw[:, iwake]
            Vm[:, ir] .+= Ar_ww[iwake, ir] * Γw[:, iwake]
        end

        # add rotor induced velocities
        for irotor in 1:size(Σr, 1)
            Vm[:, ir] .+= Ax_rw[irotor, ir] * Σr[:, irotor]
            Vm[:, ir] .+= Ar_rw[irotor, ir] * Σr[:, irotor]
        end

    end

    return Vm
end

"""
    calculate_wake_vortex_strengths!(Γw, wake_panels, Vm, Γ_tilde, H_tilde, rotor_indices)

Calculate wake vortex strengths
"""
function calculate_wake_vortex_strengths!(Γw, wake_panels, Vm, Γ_tilde, H_tilde, rotor_indices)

    # number of streamwise and radial panels
    nx, nr = size(Γw)

    # loop through each radial position
    for ir in 1:nr
        # loop through each streamwise position
        for ix in 1:nx

            # find upstream rotor
            irotor = findlast(idx -> idx <= x, rotor_indices)

            # compute average meridional velocities
            if ir == 1 || ir == nr
                # use panel velocity as average velocity
                Vm_avg = Vm[ir, ix]
            else
                # average velocities above and below the current panel
                Vm_avg = 1/2 * (Vm[ir-1, ix] + Vm[ir+1, ix])
            end

            # calculate the wake vortex strength
            if Vm_avg <= 0.0
                Γw[ir,ix] = 0.0
            else
                r = wake_panels[ir].panel_center[ix, 2]
                K = -1/(8*pi^2*r^2)
                ΔΓ2 = Γ_tilde[ir+1, irotor]^2 - Γ_tilde[ir, irotor]^2
                ΔH = H_tilde[ir+1, irotor] - H_tilde[ir, irotor]
                Γw[ir,ix] = (K*ΔΓ2 + ΔH)/Vm_avg
            end
        end
    end

    return Γw
end

function initialize_wake_vortex_strengths(Γr, Ωr, freestream)

    # set initial Vm values
    Vm = similar(Γr, nr - 1, nx) .= freestream.Vinf

    # calculate enthalpy jumps
    H_tilde = calculate_enthalpy_jumps(Γr, Ωr, B)

    # calculate net circulation
    Γ_tilde = calculate_net_circulation(Γr, B)

    # calculate wake vortex strengths
    calculate_wake_vortex_strengths!(Γw, wake_panels, Vm, Γ_tilde, H_tilde, rotor_indices)

    # println(wake_vortex_strengths)
    return wake_vortex_strengths
end