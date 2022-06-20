"""
"""
function generate_input_file(;
    # File Data
    file_name="test.case",
    version="DFDC Version 0.70E+03",
    case_name="10-2 Case with MA10x5 3-Blade",
    # Oper Section Data
    vinf=0.0, #hover
    vref=50.0, #TODO: pick a better vref
    rpm=8000.0, #TODO: pick a better RPM
    rho=1.225, # air density
    vso=340.0, # speed of sound, m/s
    rmu=0.178e-4, #mu
    alt=0.0, #altitude
    xdwake=1.0, #wake distance (in duct diameters) downstream
    nwake=20, #number of points to use in wake
    lwkrlx="F", #flag for whether to update wake grid as part of solver, "F" or "T"
    # Aero Section Data
    nsections=1, # number of rotor sections to define airfoil data for
    xisection=[0.0], #r/R coordinate along blade
    a0deg=[0.0], #zero lift angle of attack
    dclda=[2 * pi], #lift slope
    clmax=[1.5], #max lift coeff
    clmin=[-1.0], #min lift coeff
    dcldastall=[0.5], #TODO what is this?
    dcldstall=[0.2], #TODO what is this?
    cmconst=[0.0], #TODO what is this?
    mcrit=[0.7], #TODO what is this?
    cdmin=[0.012], #minimum drag coeff
    clcdmin=[0.1], # minimum cl/cd
    dcddcl2=[0.005], #cl vs cd squared slope
    reref=[0.2e6], #reference reynolds number
    reexp=[0.35], #TODO what is this
    # Rotor Section Data
    xdisc=0.35, #axial location of rotor #TODO what are the units and datum?
    nblds=3, #number of blades
    nrsta=14, #number of radial stations to use in analysis
    nstations=14, #number of stations with defined data
    r=[], #radial station locations in meters
    chord=[], #chord lengths at radial stations in meters
    beta=[], #local twist angle in degrees at radial stations
    # Drag Section Data
    pts=2, #number of points defining drag object
    x=[], #x cooridnates of points #TODO units and datum?
    r=[], #r coordinates of points #TODO units and datum?
    cda=[], #drag area TODO what does that mean??
    # Geometry Section Data
    #
    tab="    ",
)

    # Open file for writing
    f = open(file_name, "w")

    # Version
    write(f, version * "\n")

    # Case Name
    write(f, case_name * "\n\n")

    opersection!(
        f;
        vinf=vinf,
        vref=vref,
        rpm=rpm,
        rho=rho,
        vso=vso,
        rmu=rmu,
        alt=alt,
        xdwake=xdwake,
        nwake=nwake,
        lwkrlx=lwkrlx,
    )

    aerosection!(
        f;
        nsections=nsections,
        xisection=xisection,
        a0deg=a0deg,
        dclda=dclda,
        clmax=clmax,
        clmin=clmin,
        dcldastall=dcldastall,
        dcldstall=dcldstall,
        cmconst=cmconst,
        mcrit=mcrit,
        cdmin=cdmin,
        clcdmin=clcdmin,
        dcddcl2=dcddcl2,
        reref=reref,
        reexp=reexp,
        tab=tab,
    )

    rotorsection!(
        f;
        xdisc=xdisc,
        nblds=nblds,
        nrsta=nrsta,
        nstations=nstations,
        r=r,
        chord=chord,
        beta=beta,
        tab=tab,
    )

    dragsection!(f; pts=pts, x=x, r=r, cda=cda, tab=tab)

    close(f)

    return nothing
end

"""
    Operating Conditions

    This section contains the flow condition and operating point data for the
    case.  The first line contains the freestream velocity, reference velocity
    and RPM (used for actuator disk to estimate swirl effects, or a rotor for
    operating condition).  The second line contains the fluid properties and/or
    the altitude (used to set atmosphere properties if they are otherwise set to
    0.0).  The next line contains the downstream vortex wake length to be used
    for the flow domain (in duct diameters, use 0.75 to 2 diameters for most
    cases).   The number of points to use in the wake is set here (20 works well).
    The last input is a flag to set automatic wake realignment and regridding,
    normally not used as it slows the case convergence).
"""
function opersection!(
    f;
    vinf=0.0, #hover
    vref=50.0, #TODO: pick a better vref
    rpm=8000.0, #TODO: pick a better RPM
    rho=1.225, # air density
    vso=340.0, # speed of sound, m/s
    rmu=0.178e-4, #mu
    alt=0.0, #altitude
    xdwake=1.0, #wake distance (in duct diameters) downstream
    nwake=20, #number of points to use in wake
    lwkrlx="F", #flag for whether to update wake grid as part of solver, "F" or "T"
    tab="    ",
)
    write(f, "OPER\n")
    write(f, "!$(tab)Vinf$(tab*tab)Vref$(tab*tab)RPM\n")
    write(f, "$tab$vinf$(tab*tab)$vref$(tab*tab)$rpm\n")
    write(f, "!$(tab)Rho$(tab*tab)Vso$(tab*tab)Rmu$(tab*tab)Alt\n")
    write(f, "$tab$rho$(tab*tab)$vso$(tab*tab)$rmu$(tab*tab)$alt\n")
    write(f, "!$(tab)XDwake$(tab*tab)Nwake\n")
    write(f, "$tab$xdwake$(tab*tab)$nwake\n")
    write(f, "!$(tab*tab)Lwkrlx\n")
    write(f, "$(tab*tab)$lwkrlx\n")
    write(f, "ENDOPER\n\n")

    return nothing
end

"""
    Aero Section

    This section starts with the AERO keyword and ends with ENDAERO.
    This section contains aerodynamic data that is used for each blade element in
    a rotor analysis or design.  The data is used for a parametric model of the
    lifting and drag properties of an airfoil section (2D section).  The details
    of the model are presented in the aero.f code sections but are derived from
    the XROTOR aero model for blade elements.

    The aero properties are used for bladed disk analysis for the blade elements.
    They are also used when a bladed rotor or stator is designed.

    Any number (one or more) aero "sections" are used, located on the blade by
    a XIsection coordinate (XI = r/R). If two or more aero "sections are specifed
    the aero properties for any blade element at XI=xi_i will be obtained by
    interpolation of aero properties of aero "sections" whose XIsection values
    bound the desired blade element radial station.  If only one section is present
    its aero properties will be used for the whole blade from hub to tip.

    Note that if one rotor disk is used only one set of aero properties are
    needed (or none at all, in which case a default set is used). The AERO
    data should precede the associated rotor or actuator disk definition.

    If a case has two disks (rotor/stator or ?) each AERO/disk definition is
    treated as a composite definition, i.e. the first aero data is assumed
    to refer to the first disk defined (ACTDISK or ROTOR), the second aero data
    is assumed to refer to the second disk defined (ACTDISK or ROTOR).  The AERO
    data should normally precede its intended disk definition.

    Note that, for the current release, the aero properties must be altered for
    a disk with negative circulation (BGAM).  DFDC does not automatically "flip"
    the aero properties for a blade with negative circulation.  The negative
    circulation means that negative CL's will generate the negative circulation
    (as the chord is never allowed to go negative!!).  This means the user should
    input the aero properties with CLmax, CLmin, CLCDmin reversed.  This means that
    -CLmax is input as CLmin, -CLmin is input as CLmax and CLCDmin=-CLCDmin.

    The plotted or output print of disk properties will show negative CL's for
    disks with negative circulation.
    """
function aerosection!(
    f;
    nsections=1, # number of rotor sections to define airfoil data for
    xisection=[0.0], #r/R coordinate along blade
    a0deg=[0.0], #zero lift angle of attack
    dclda=[2 * pi], #lift slope
    clmax=[1.5], #max lift coeff
    clmin=[-1.0], #min lift coeff
    dcldastall=[0.5], #TODO what is this?
    dcldstall=[0.2], #TODO what is this?
    cmconst=[0.0], #TODO what is this?
    mcrit=[0.7], #TODO what is this?
    cdmin=[0.012], #minimum drag coeff
    clcdmin=[0.1], # minimum cl/cd
    dcddcl2=[0.005], #cl vs cd squared slope
    reref=[0.2e6], #reference reynolds number
    reexp=[0.35], #TODO what is this
    tab="    ",
)
    write(f, "AERO\n")
    write(f, "!$(tab)#sections\n")
    write(f, "$tab$nsections\n")

    for i in 1:nsections
        write(f, "!$(tab)Xisection\n")
        write(f, "$tab$(xisection[i])\n")
        write(f, "!$(tab)A0deg$(tab*tab)dCLdA$(tab*tab)CLmax$(tab*tab)CLmin\n")
        write(
            f, "$tab$(a0deg[i])$tab$tab$(dclda[i])$tab$tab$(clmax[i])$tab$tab$(clmin[i])\n"
        )
        write(f, "!$(tab)dCLdAstall$(tab*tab)dCLstall$(tab*tab)Cmconst$(tab*tab)Mcrit\n")
        write(
            f,
            "$tab$(dcldastall[i])$tab$tab$(dcldstall[i])$tab$tab$(cmconst[i])$tab$tab$(mcrit[i])\n",
        )
        write(f, "!$(tab)CDmin$(tab*tab)CLCDmin$(tab*tab)dCDdCL^2\n")
        write(f, "$tab$(cdmin[i])$tab$tab$(clcdmin[i])$tab$tab$(dcddcl2[i])\n")
        write(f, "!$(tab)REref$(tab*tab)REexp\n")
        write(f, "$tab$(reref[i])$tab$tab$(reexp)\n")
    end

    write(f, "ENDAERO\n\n")

    return nothing
end

"""
    Rotor Section

    This section contains a specification for a bladed disk in the ducted fan.
    The section starts with the ROTOR keyword and ends with ENDROTOR.

    The rotor is specified by an axial location for the rotor (Xdisk) which must
    be located between the duct and centerbody leading and trailing edges.  The
    number of rotor blades (Nblds) is specified here.
    The number of points to use to represent the rotor blade elements (NRPdef)
    is also required (normallly 11-15 points is enough).

    Adding rotor/actuator disk analysis points slows the method as each rotor
    radial station emits a vortex wake (which has Nwake+ points on each).

    The rotor is specified by a set of radial stations with specified
    blade chord and blade angle.  The blade angle is measured from the plane
    of rotor rotation (0 for flat pitch, 90deg for axial inflow direction). The
    radial station and chord are assumed to be dimensional in meters.

    Note that the user can specify any number of radial stations that characterize
    the rotor, these will be interpolated to the rotor analysis stations.  Rotors,
    once defined, take precedence over actuator disks. Normally you don't use both.
"""
function rotorsection!(
    f;
    xdisc=0.35, #axial location of rotor #TODO what are the units and datum?
    nblds=3, #number of blades
    nrsta=14, #number of radial stations to use in analysis
    nstations=14, #number of stations with defined data
    r=[], #radial station locations in meters
    chord=[], #chord lengths at radial stations in meters
    beta=[], #local twist angle in degrees at radial stations
    tab="    ",
)
    write(f, "ROTOR\n")
    write(f, "!$(tab)Xdisk$(tab*tab)Nblds$(tab*tab)NRsta\n")
    write(f, "$tab$xdisk$tab$tab$nblds$tab$tab$nrsta\n")
    write(f, "!$(tab)#stations\n")
    write(f, "$tab$nstations\n")
    write(f, "!$(tab)r$(tab*tab)Chord$(tab*tab)Beta\n")
    for i in 1:nstations
        write(f, "$tab$(r[i])$tab$tab$(chord[i])$tab$tab$(beta[i])\n")
    end

    write(f, "ENDROTOR\n\n")

    return nothing
end

"""
    Drag Object Section

    This section starts with the DRAGOBJ keyword and ends with ENDDRAGOBJ.
    This section contains a CDA (drag area) and X,Y (x,r in axisymmetric system)
    coordinates of the drag line.  This is used to set up a line source that will
    represent the blockage and loss of the drag object on the flow field.

    The DRAG menu in DFDC can be used to edit, change or add drag objects.  Drag
    objects are saved in the case file when it is written.
"""
function dragsection!(
    f;
    pts=2, #number of points defining drag object
    x=[], #x cooridnates of points #TODO units and datum?
    r=[], #r coordinates of points #TODO units and datum?
    cda=[], #drag area TODO what does that mean??
    tab="    ",
)
    write(f, "DRAGOBJ\n")
    write(f, "!$(tab)#pts\n")
    write(f, "pts\n")
    write(f, "!$(tab)x$(tab*tab)r$(tab*tab)CDA\n")
    for i in 1:pts
        write(f, "$tab$(x[i])$tab$tab$(r[i])$tab$tab$(cda[i])\n")
    end
    write(f, "ENDDRAGOBJ")

    return nothing
end

generate_input_file()
