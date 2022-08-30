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

    #generate rotor object
    rotor1 = DuctTAPE.RotorGeometry(
        xdisk1, #x position of rotor
        nblade1, #number of blades
        rnondim1, #radial stations
        0.0, #tip gap
        chord1, #chords
        beta1, #twists
        nothing, #skews
        nothing, #rakes
        fill(af1,length(rnondim1)), #airfoils
        nothing, #solidities
        rpm, #RPM
    )

    #generate stator object (rpm is zero for stator)
    rotor2 = DuctTAPE.RotorGeometry(
        xdisk2, #x position of rotor
        nblade2, #number of blades
        rnondim2, #radial stations
        0.0, #tip gap
        chord2, #chords
        beta2, #twists
        nothing, #skews
        nothing, #rakes
        fill(af2,length(rnondim2)), #airfoils
        nothing, #solidities
        0.0, #RPM
    )

    #return array of rotor objects
    return [rotor1; rotor2]
end
