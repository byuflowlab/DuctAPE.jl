
"""
"""
function setup_rotorgrid_aero(rotors, blades, wakegrid, rotor_source_panels)

    #set freestream from values in dfdc case file
    freestream = DuctTAPE.Freestream(vinf, vref, rho, vso, rmu)

    #initialize system aerodynamics
    return initialize_system_aerodynamics(
        rotors, blades, wakegrid, rotor_source_panels, freestream; niter=10, rlx=0.5
    )
end

