#=
Coupling of FLOWFoil to CCBlade
=#

"""
    ff2ccb(ff_solution, rotor, freestream)

Takes the FLOWFoil solution object, rotor geometry object, and freestream object and genereates the inputs to the CCBlade solve function.

**Arguments:**
- `ff_solution::FLOWFoil.InviscidSolution` : Inviscid solution object from FLOWFoil.
- `rotor::RotorGeometry` : rotor geometry object
- `freestream::Freestream` : freestream object

**returns:**
- `ccbrotor::ccblade.rotor` : ccblade rotor object
- `ccbsections::array{ccblade.section}` : array of ccblade section objects
- `ccbop::array{ccblade.operatingpoint}` : array of ccblade operation points
- `rdist::array{float}` : array of refined, dimensional rotot radial station r locations.
"""
function ff2ccb(ff_solution, rotor, freestream; debug=false)

    # need to get wall mesh geometry in order to define rotor blade dimensions
    xduct, yduct, xhub, yhub = extract_ff_geom(ff_solution)

    # duct wall geometry needs to be split before being splined
    xductinner, _, yductinner, _ = split_wall(xduct, yduct)

    #spline duct and hub geometries
    ductspline = FLOWMath.akima(reverse(xductinner), reverse(yductinner))
    hubspline = FLOWMath.akima(xhub, yhub)

    #find rhub from point on hub spline corresponding with rotor xlocation
    rhub = hubspline(rotor.xlocation)

    #find rtip from point on duct spline corresponding with rotor xlocation
    rtip = ductspline(rotor.xlocation)

    #get dimensional rotor radial station locations
    dim_rad_stash = lintran(
        rhub, #range a start
        rtip, #range a end
        rotor.radialstations[1], #range b start
        rotor.radialstations[end], #range b end
        rotor.radialstations, #range to transform
    )

    #assemble field points
    field_points = [[rotor.xlocation; dim_rad_stash[i]] for i in 1:length(dim_rad_stash)]

    # probe ff_solution velocity field at rotor x location
    duct_induced_velocities = FLOWFoil.probe_velocity_axisym(ff_solution, field_points)

    #get full magnitude velocities
    velocities = freestream.vinf .* duct_induced_velocities
    for i in 1:length(velocities)
        velocities[i][1] += freestream.vinf
    end

    #define rotor object
    ccbrotor = ccb.rotor(rhub, rtip, rotor.numblades; turbine=false)

    #define section objects
    ccbsections, rdist = generate_ccb_sections(rotor, dim_rad_stash)

    #convert velocities
    vx, vy = ff2ccb_velocity(rotor, dim_rad_stash, velocities, rdist)

    #define operating point
    #note: for single rotor, v_x = vinf and v_y = omega*r
    pitch = 0.0
    op =
        ccb.operatingpoint.(vx, vy, freestream.rho, pitch, freestream.mu, freestream.asound)

    if debug
        ccbrotor,
        ccbsections,
        ccb.simple_op.(freestream.vinf, get_omega(rotor.rpm), rdist, freestream.rho),
        rdist
    else
        return ccbrotor, sections, op, rdist
    end
    #fyi this is how to call ccblade
    # ccb_out = ccb.solve.(ref(ccbrotor), sections, op)
end

"""
    extract_ff_geom(ff_solution)

extract x and r coordiantes of duct and hub geometries from the FLOWFoil solution object.

**arguments:**
- `ff_solution::FLOWFoil.inviscidsolution` : inviscid solution object from FLOWFoil.

**returns:**
- `duct_x::array{float}` : array of x-coordinates of duct wall control points.
- `duct_r::array{float}` : array of r-coordinates of duct wall control points.
- `hub_x::Array{Float}` : Array of x-coordinates of hub wall control points.
- `hub_r::Array{Float}` : Array of r-coordinates of hub wall control points.
"""
function extract_ff_geom(ff_solution)

    #rename for convenience
    mesh1 = ff_solution.meshes[1]
    mesh2 = ff_solution.meshes[2]

    if mesh1.panels[1].controlpoint[2] > mesh2.panels[1].controlpoint[2]

        #duct is mesh1 since TE point is above hub LE point
        cpduct = (panel -> panel.controlpoint).(mesh1.panels)
        cphub = (panel -> panel.controlpoint).(mesh2.panels)

    else

        #otherwise the duct is the second mesh
        cpduct = (panel -> panel.controlpoint).(mesh2.panels)
        cphub = (panel -> panel.controlpoint).(mesh1.panels)
    end

    return getindex.(cpduct, 1),
    getindex.(cpduct, 2), getindex.(cphub, 1),
    getindex.(cphub, 2)
end

"""
    generate_ccb_sections(rotor, dim_rad_stash)

Generates CCBlade section objects from rotor information and dimensional radial stations.

**Arguments:**
- `rotor::RotorGeometry` : Rotor geometry object
- `dim_rad_stash::Array{Float}` : Dimensional Radial Station locations

**Returns:**
- `ccbsections::Array{CCBlade.Section}` : Array of CCBlade section objects
- `rdist::array{float}` : array of refined, dimensional rotot radial station r locations.
"""
function generate_ccb_sections(rotor, dim_rad_stash)

    #rename for convenience
    chords = rotor.chords * dim_rad_stash[end] #dimensionalize
    twists = rotor.twists * pi / 180.0 #convert to radians
    airfoils = rotor.airfoils

    ## -- Refine Blade
    #spline chord and twist
    csp = FLOWMath.Akima(dim_rad_stash, chords)
    tsp = FLOWMath.Akima(dim_rad_stash, twists)

    #get refined radial stations
    rdist = range(dim_rad_stash[1], dim_rad_stash[end]; length=rotor.nref)

    #get refined chord and twist
    cdist = csp.(rdist)
    tdist = tsp.(rdist)

    ## -- assign airfoils
    #initialize array
    afdist = Array{typeof(airfoils[1])}(undef, rotor.nref)

    #get average values of radial stations to better place airfoils (center defined airfoils in refined blade)
    mean_rad_stash = (dim_rad_stash[1:(end - 1)] .+ dim_rad_stash[2:end]) ./ 2.0

    #loop through refined radial stations and apply appropriate airfoil
    for i in 1:(rotor.nref)
        ridx = findfirst(x -> x > rdist[i], mean_rad_stash)
        if ridx != nothing
            afdist[i] = airfoils[ridx]
        else
            afdist[i] = airfoils[end]
        end
    end

    return ccb.Section.(rdist, cdist, tdist, afdist), rdist
end

"""
    ff2ccb_velocity(rotor, dim_rad_stash, velocities, rdist)

Convert velocities from FLOWFoil to the CCBlade reference frame.

**Arguments:**
- `rotor::RotorGeometry` : Rotor geometry object
- `dim_rad_stash::Array{Float}` : Dimensional Radial Station locations
- `velocities::Array{Array{Float}}` : Array of [x; r] velocities from FLOWFoil.
- `rdist::array{float}` : array of refined, dimensional rotot radial station r locations.

**Returns:**
- `Vx::Array{Float}` : x-velocities in the CCBlade reference frame
- `Vy::Array{Float}` : y-velocities in the CCBlade reference frame
"""
function ff2ccb_velocity(rotor, dim_rad_stash, velocities, rdist)

    #spline velocities
    vxsp = FLOWMath.Akima(dim_rad_stash, getindex.(velocities, 1))
    vrsp = FLOWMath.Akima(dim_rad_stash, getindex.(velocities, 2))
    vtsp = FLOWMath.Akima(dim_rad_stash, get_omega(rotor.RPM) .* dim_rad_stash)

    #get refined velocities
    Vx = vxsp.(rdist)
    Vy = vtsp.(rdist)

    return Vx, Vy
end
