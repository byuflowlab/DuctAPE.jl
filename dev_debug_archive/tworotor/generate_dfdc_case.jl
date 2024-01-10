using DuctAPE
const dt = DuctAPE

using FLOWMath

"""
function for generating DFDC case files
"""
function gen_dfdc_case(
    filename,
    op_data,
    wake_data,
    airfoil_data,
    rotor_data,
    hub_coordinates,
    duct_coordinates;
    savepath="dev_debug_archive/tworotor/",
    version=0.7,
    case_name="Custom Test",
)

    #open file
    f = open(savepath * filename, "w")

    # write header
    write(f, "DFDC Version $(version)\n")
    write(f, case_name * "\n")
    write(f, "\n")

    # write operating conditions
    write(f, "OPER\n")
    write(f, "!        Vinf         Vref          RPM\n")
    write(f, "    $(op_data.Vinf)     $(op_data.Vref)     $(op_data.RPM[1])\n")
    write(f, "!         Rho          Vso          Rmu           Alt\n")
    write(f, "    $(op_data.rho)     $(op_data.Vso)     $(op_data.mu)     $(op_data.Alt)\n")
    #NOTE: Nwake is the number of wake x-stations aft of the duct trailing edge, it is used in calculating the wake expansion spacing in DFDC, the DuctAPE equivalent is the last entry in the npanels vector
    write(f, "!       XDwake        Nwake\n")
    write(f, "    $(wake_data.xwake)     $(wake_data.nwake)\n")
    write(f, "!       Lwkrlx\n")
    write(f, "    " * wake_data.rlx_wake * "\n")
    write(f, "ENDOPER\n")
    write(f, "\n")

    for irotor in 1:length(rotor_data)
        # write airfoil properties
        write(f, "AERO\n")
        write(f, "!  #sections\n")
        write(f, "    $(rotor_data[irotor].naf)\n")
        if rotor_data[irotor].naf == 1
            write_af_section(f, airfoil_data[irotor], 1)
        else
            for iaf in 1:(rotor_data[irotor].naf)
                write_af_section(f, airfoil_data[irotor], iaf)
            end
        end
        write(f, "ENDAERO\n")
        write(f, "\n")

        # write rotor section
        write(f, "ROTOR\n")
        write(f, "!       Xdisk        Nblds       NRPdef\n")
        write(
            f,
            "    $(rotor_data[irotor].rotorzloc)     $(rotor_data[irotor].B)     $(wake_data.nwake_sheets)\n",
        )
        write(f, "!  #stations\n")
        write(f, "    $(length(rotor_data[irotor].r))\n")
        write(f, "!           r        Chord         Beta\n")
        for ir in 1:length(rotor_data[irotor].r)
            write(
                f,
                "    $(rotor_data[irotor].r[ir])     $(rotor_data[irotor].chord[ir])     $(rotor_data[irotor].twist[ir])\n",
            )
        end
        write(f, "ENDROTOR\n")
        write(f, "\n")
    end # for number of rotors

    # write body coordintes
    write(f, "GEOM\n")
    write(f, case_name * "\n")
    for ih in 1:length(hub_coordinates[:, 1])
        write(f, "    $(hub_coordinates[ih,1])     $(hub_coordinates[ih,2])\n")
    end
    write(f, "  999.0 999.0\n")
    for id in 1:length(duct_coordinates[:, 1])
        write(f, "    $(duct_coordinates[id,1])     $(duct_coordinates[id,2])\n")
    end
    write(f, "ENDGEOM\n")
    close(f)

    return nothing
end

function write_af_section(f, airfoil_data, idx)
    write(f, "!   Xisection\n")
    write(f, "    $(airfoil_data.xisection)\n")
    write(f, "!       A0deg        dCLdA        CLmax         CLmin\n")
    write(
        f,
        "    $(airfoil_data.alpha0)     $(airfoil_data.dclda)     $(airfoil_data.clmax)     $(airfoil_data.clmin)\n",
    )
    write(f, "!  dCLdAstall     dCLstall      Cmconst         Mcrit\n")
    write(
        f,
        "    $(airfoil_data.dclda_stall)     $(airfoil_data.dcl_stall)     $(airfoil_data.cmcon)     $(airfoil_data.mcrit)\n",
    )
    write(f, "!       CDmin      CLCDmin     dCDdCL^2\n")
    write(
        f,
        "    $(airfoil_data.cdmin)     $(airfoil_data.clcdmin)     $(airfoil_data.dcddcl2)\n",
    )
    write(f, "!       REref        REexp\n")
    write(f, "    $(airfoil_data.Re_ref)     $(airfoil_data.Re_exp)\n")

    return nothing
end

function test_gen_dfdc_case()
    filename = "test.case"

    op_data = (; rho=1.226, mu=1.78e-5, Vso=340.0, Vinf=0.0, Vref=50.0, Alt=0.0, RPM=[8000.0; 0.0])

    wake_data = (; nwake_sheets=11, nwake=20, xwake=0.8, rlx_wake="F\n")

    airfoil_data = [(;
        xisection=0.0,
        alpha0=0.0000,
        dclda=6.2800,
        clmax=1.5000,
        clmin=-1.0000,
        dclda_stall=0.50000,
        dcl_stall=0.20000,
        cmcon=0.0000,
        mcrit=0.70000,
        cdmin=0.12000E-01,
        clcdmin=0.10000,
        dcddcl2=0.50000E-02,
        Re_ref=0.20000E+06,
        Re_exp=0.35000,
    ),
    (;
        xisection=0.0,
        alpha0=0.0000,
        dclda=6.2800,
        clmax=1.000,
        clmin=-1.5000,
        dclda_stall=0.50000,
        dcl_stall=0.20000,
        cmcon=0.0000,
        mcrit=0.70000,
        cdmin=0.12000E-01,
        clcdmin=0.10000,
        dcddcl2=0.50000E-02,
        Re_ref=0.20000E+06,
        Re_exp=0.35000,
    )]


    rct1 = [
  0.50494E-01  0.56234E-01   70.828
  0.61571E-01  0.53641E-01   59.758
  0.72648E-01  0.51101E-01   51.707
  0.83724E-01  0.48825E-01   46.429
  0.94801E-01  0.46895E-01   42.209
  0.10588      0.45316E-01   38.766
  0.11695      0.44065E-01   35.897
  0.12803      0.43111E-01   33.457
  0.13911      0.42424E-01   31.337
  0.15019      0.41981E-01   29.453
    ]

    rct2 =[
  0.50488E-01  0.58025E-01   92.711
  0.61734E-01  0.57606E-01   91.739
  0.72981E-01  0.57052E-01   91.228
  0.84242E-01  0.56419E-01   90.868
  0.95517E-01  0.55738E-01   90.638
  0.10681      0.55032E-01   90.441
  0.11811      0.54317E-01   90.124
  0.12943      0.53608E-01   89.583
  0.14077      0.52918E-01   88.955
  0.15213      0.52256E-01   88.388
 ]

    rotor_data = [(;
        naf=1, rotorzloc=0.12, B=6, r=rct1[:, 1], chord=rct1[:, 2], twist=rct1[:, 3]
    ),(;
        naf=1, rotorzloc=0.22, B=9, r=rct2[:, 1], chord=rct2[:, 2], twist=rct2[:, 3]
       )]

    hub_coordinates = [
        0.306379 0.035928
        0.296701 0.039168
        0.285768 0.041795
        0.274184 0.043576
        0.262095 0.044564
        0.248578 0.044839
        0.232723 0.044847
        0.216580 0.044866
        0.200407 0.044882
        0.184232 0.044898
        0.168060 0.044913
        0.151899 0.044930
        0.135773 0.044950
        0.120000 0.044952
        0.109361 0.044999
        0.099680 0.044922
        0.090039 0.044561
        0.080630 0.043916
        0.071451 0.042966
        0.062540 0.041690
        0.053933 0.040075
        0.045705 0.038107
        0.037900 0.035771
        0.030604 0.033068
        0.023901 0.030001
        0.017880 0.026585
        0.012632 0.022848
        0.008231 0.018825
        0.004736 0.014551
        0.002179 0.010047
        0.000586 0.005293
        0.000000 0.000000
    ]
    duct_coordinates = [
        0.304542 0.159526
        0.296440 0.162842
        0.285240 0.167270
        0.273301 0.171796
        0.261206 0.176188
        0.249026 0.180416
        0.236794 0.184470
        0.224549 0.188335
        0.212269 0.192017
        0.200020 0.195490
        0.187807 0.198746
        0.175647 0.201771
        0.163586 0.204547
        0.151644 0.207051
        0.139859 0.209267
        0.128380 0.211153
        0.117508 0.212618
        0.106950 0.213649
        0.096610 0.214255
        0.086499 0.214421
        0.076647 0.214136
        0.067113 0.213390
        0.057935 0.212174
        0.049193 0.210487
        0.040949 0.208332
        0.033312 0.205736
        0.026403 0.202735
        0.020339 0.199393
        0.015227 0.195790
        0.011135 0.192004
        0.008090 0.188086
        0.006112 0.184067
        0.005242 0.180005
        0.005484 0.176154
        0.006854 0.172546
        0.009324 0.169289
        0.012842 0.166404
        0.017419 0.163862
        0.023110 0.161648
        0.029956 0.159771
        0.037937 0.158256
        0.046983 0.157103
        0.057025 0.156294
        0.067995 0.155792
        0.079836 0.155546
        0.092531 0.155498
        0.106044 0.155585
        0.120000 0.155721
        0.134222 0.155902
        0.148679 0.156177
        0.163490 0.156523
        0.178507 0.156897
        0.193399 0.157258
        0.208123 0.157586
        0.222751 0.157864
        0.237332 0.158088
        0.251898 0.158254
        0.266505 0.158365
        0.281130 0.158423
        0.294972 0.158441
        0.304466 0.158439
    ]

    gen_dfdc_case(
        filename,
        op_data,
        wake_data,
        airfoil_data,
        rotor_data,
        hub_coordinates,
        duct_coordinates;
        savepath=savepath,
    )

 write_ducttape_params(
    filename*".jl",
    op_data,
    wake_data,
    airfoil_data,
    rotor_data,
    reverse(hub_coordinates,dims=1),
    reverse(duct_coordinates,dims=1);
    savepath=savepath,
    npanels_interior=80,
)

    return nothing
end

"""
function for generating DuctAPE input parameters
"""
function write_ducttape_params(
    filename,
    op_data,
    wake_data,
    airfoil_data,
    rotor_data,
    hub_coordinates,
    duct_coordinates;
    savepath="dev_debug_archive/tworotor/",
    npanels_interior=80,
)
    f = open(savepath * filename, "w")

    rotorzloc = rotor_data.rotorzloc
    write(f, "rotorzloc = $rotorzloc\n")
    B = rotor_data.B
    write(f, "B = $B\n")

    #dimensional radius
    write(f, "r = $(reshape(reduce(vcat,rotor_data.r),(length(rotor_data[1].r),length(rotor_data.r))))\n")

    #dimensional chord
    write(f, "chords = $(reshape(reduce(vcat,rotor_data.chord),(length(rotor_data[1].chord),length(rotor_data.chord))))\n")

    #twist in degrees converted to radians
    write(f, "twists = $(pi/180.0*reshape(reduce(vcat,rotor_data.twist),(length(rotor_data[1].twist),length(rotor_data.twist))))\n")

    for ir in 1:length(rotorzloc)
    #Airfoil Data:
    (; alpha0, clmax, clmin, dclda, dclda_stall, dcl_stall, cdmin, clcdmin, dcddcl2, cmcon, Re_ref, Re_exp, mcrit) = airfoil_data[ir]

    write(
        f,
        "afparams$ir = dt.DFDCairfoil(
   $alpha0,
   $clmax,
   $clmin,
   $dclda,
   $dclda_stall,
   $dcl_stall,
   $cdmin,
   $clcdmin,
   $dcddcl2,
   $cmcon,
   $Re_ref,
   $Re_exp,
   $mcrit,
   )\n",
    )

    write(f, "airfoils$ir = fill(afparams$ir, length(r[:,$ir]))\n")

    end

    write(f, "airfoils = [")
    for ir in 1:length(rotorzloc)
        write(f, "airfoils$ir ")
    end
    write(f, "]\n")

    #---------------------------------#
    #       Duct and Hub Geometry     #
    #---------------------------------#

    _, duct_leidx = findmin(duct_coordinates[:, 1])
    ductxin = reverse(duct_coordinates[1:duct_leidx, 1])
    ductrin = reverse(duct_coordinates[1:duct_leidx, 2])

    # load in duct and hub geometry, spline, and find out what the duct and hub radii are at the rotor positions to figure out what Rtip and Rhub are.
    Rhub = FLOWMath.akima(hub_coordinates[:, 1], hub_coordinates[:, 2], rotorzloc)
    Rtip = FLOWMath.akima(ductxin, ductrin, rotorzloc)

    write(f, "Rtip=$(Rtip[1])\n")
    write(f, "Rhub=$(Rhub[1])\n")

    #---------------------------------#
    #      Operation Conditions       #
    #---------------------------------#

    (; rho, mu, Vso, Vinf, Vref, Alt, RPM) = op_data

    write(f, "rho=$rho\n")
    write(f, "mu=$mu\n")
    write(f, "Vinf=$Vinf\n")
    write(f, "Vref=$Vref\n")
    write(f, "Omega = $RPM * pi / 30  # convert from RPM to rad/s\n")
    write(f, "asound = $Vso\n")

    #---------------------------------#
    #        Paneling Options         #
    #---------------------------------#
    (; xwake, nwake_sheets) = wake_data

    wake_length = xwake #times duct chord

    discscale = 1.0

    ductle = minimum(duct_coordinates[:, 1])
    ductte = maximum(duct_coordinates[:, 1])
    huble = minimum(hub_coordinates[:,1])
    hubte = maximum(hub_coordinates[:,1])
    ductchord = maximum(max(hubte, ductte) - min(huble,ductle))

    if isapprox(ductte, hubte)
    xstations = [ductle; rotorzloc; ductte]
    else
    xstations = [ductle; rotorzloc; min(hubte,ductte); max(hubte,ductte)]
    end

    nhub_inlet = round(Int, npanels_interior * discscale * (rotorzloc[1]-ductle)/ductchord)

    nduct_inlet = nhub_inlet

    npanels_aft = [ceil(Int, npanels_interior * discscale * (xstations[i+1]-xstations[i])/ductchord) for i in 2:length(xstations)-1]

    nwake = round(Int, npanels_interior * discscale)

    append!(npanels_aft,nwake)

    write(f, "nhub_inlet = $nhub_inlet\n")
    write(f, "nduct_inlet = $nduct_inlet\n")
    write(f, "nwake_sheets = $nwake_sheets\n")
    write(f, "wake_length = $xwake\n")
    write(f, "npanels = $npanels_aft\n")

    #--------------------------------#
    #      Assemble Named Tuples      #
    #---------------------------------#

    # Rotor Parameters
    write(
        f,
        "rotor_parameters = [(;
        rotorzloc=rotorzloc[i],
   nwake_sheets,
   r=r[:,i] ./ Rtip, #non-dimensionalize
   chords=chords[:,i],
   twists=twists[:,i],
   airfoils = airfoils[:,i],
   Rtip,
   Rhub,
   tip_gap=0.0,
   B=B[i],
   Omega=Omega[i],
   ) for i in 1:length(Omega)]\n",
    )

    # Paneling Parameters
    write(
        f,
        "paneling_constants = (; npanels, nhub_inlet, nduct_inlet, wake_length, nwake_sheets)\n",
    )

    # Freestream Parameters
    write(f, "freestream = (; rho, mu, asound, Vinf)\n")

    write(f, "reference_parameters = (; Vref, Rref=Rtip)\n")

    write(f, "duct_coordinates = $duct_coordinates\n")
    write(f, "hub_coordinates = $hub_coordinates\n")

    close(f)
    return nothing
end
