"""
"""
function setup_rotors(; datapath="data/dfdc/airfoils/")

    # generate airfoils for rotor
    af1 = DuctTAPE.ccb.AlphaReAF([
        datapath * "disk1_re5e5.dat",
        datapath * "disk1_re1e6.dat",
        datapath * "disk1_re1.5e6.dat",
        datapath * "disk1_re2e6.dat",
    ])

    # generate airfoils for stator
    af2 = DuctTAPE.ccb.AlphaReAF([
        datapath * "disk2_re5e5.dat",
        datapath * "disk2_re1e6.dat",
        datapath * "disk2_re1.5e6.dat",
        datapath * "disk2_re2e6.dat",
    ])

    #generate rotor object using initialization function
    rotor1 = initialize_rotor_geometry(
        xdisk1, #x position of rotor
        nblade1, # number of blades
        length(rnondim1), #number of radial stations
        chord1, #chords
        beta1, #twists
        fill(af1, length(rnondim1)), #airfoils
        rpm; #RPM
        radialstations=rnondim1, #radial station locations
    )

    #generate stator object (rpm is zero for stator) using direct definition
    rotor2 = DuctTAPE.RotorGeometry(
        xdisk2, #x position of rotor
        nblade2, #number of blades
        rnondim2, #radial stations
        0.0, #tip gap
        chord2, #chords
        beta2, #twists
        nothing, #skews
        nothing, #rakes
        fill(af2, length(rnondim2)), #airfoils
        -1 .* ones(Int, length(rnondim2)), #solidities
        0.0, #RPM
    )

    #return array of rotor objects
    return [rotor1; rotor2]
end

"""
"""
function setup_blades(ductgeometry, ductsplines, rotors)

    #get blade objects
    blade1 = DuctTAPE.initialize_blade_dimensions(ductgeometry, ductsplines, rotors[1])
    blade2 = DuctTAPE.initialize_blade_dimensions(ductgeometry, ductsplines, rotors[2])

    #return blades
    return [blade1; blade2]
end
