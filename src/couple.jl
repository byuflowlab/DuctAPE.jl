#=
Coupling of FLOWFoil to CCBlade
=#

"""
    ff2ccb(ff_solution, rotor,

Takes the FLOWFoil solution object, rotor geometry object, and object and genereates the inputs to the CCBlade solve function.

**Arguments:**
- `ff_solution::FLOWFoil.InviscidSolution` : Inviscid solution object from FLOWFoil.
- `dtrotor::RotorGeometry` : rotor geometry object
- `: : object

**returns:**
- `ccbrotor::CCBlade.rotor` : CCBlade rotor object
- `ccbsections::Array{CCBlade.section}` : Array of CCBlade section objects
- `ccbops:Array{Array{CCBlade.operatingpoint}}` : Array of Arrays of CCBlade operation points for each rotor RPM at which to analyze
- `v_tip::Array{Float}` : Array of tip velocities for each advance ratio
"""
function ff2ccb(dtrotor, vinf, area_ratio; debug=false)

    ## need to get wall mesh geometry in order to define dtrotor blade dimensions
    #xduct, yduct, xhub, yhub = extract_ff_geom(ff_solution)

    ## duct wall geometry needs to be split before being splined
    #xductinner, _, yductinner, _ = split_wall(xduct, yduct)

    ## spline inner duct and hub wall geometries front to back
    #ductspline = FLOWMath.Akima(reverse(xductinner), reverse(yductinner))
    #hubspline = FLOWMath.Akima(xhub, yhub)

    ## find rhub from point on hub spline corresponding with dtrotor xlocation
    #rhub = hubspline(dtrotor.xlocation)

    ## find rtip from point on duct spline corresponding with dtrotor xlocation
    #rtip = ductspline(dtrotor.xlocation)

    ###sanity check on dtrotor tip radius
    ##@assert isapprox(rtip, dtrotor.Rtip, atol=1e-5)

    ## get dimensional dtrotor radial station locations using utility function in this package
    #dtrotor.radialstations = lintran(
    #    rhub, # range a start
    #    rtip, # range a end
    #    dtrotor.radialstations[1], # range b start
    #    dtrotor.radialstations[end], # range b end
    #    dtrotor.radialstations, # range to transform
    #)

    # assemble field points, i.e. blade element locations
    field_points = [
        [dtrotor.xlocation; dtrotor.radialstations[i]] for
        i in 1:length(dtrotor.radialstations)
    ]

    TR = eltype(dtrotor.RPM)

    # define CCBlade rotor object with dummy tip correction (manually added to CCBlade source code)
    ccbrotor = ccb.Rotor(
        dtrotor.radialstations[1],
        dtrotor.radialstations[end],
        dtrotor.numblades;
        turbine=false,
        tip=CCBlade.DuctTip(),
        # mach=CCBlade.PrandtlGlauert(),
    )

    # define CCBlade Section objects
    ccbsections, rdist = generate_ccb_sections(dtrotor)

    # - Define Operation Points
    # rename for convenience
    # nops = length(vinf)
    nops = length(dtrotor.RPM)

    # initialize operating point Arrays
    # ccbops = [Vector{CCBlade.OperatingPoint{TR}}(undef, length(rdist)) for i in 1:nops]

    #initialize axial velocity arays
    v_tip = Vector{TR}(undef, nops)
    # v_tip = Vector{eltype(freestream.vinf)}(undef, nops)

    #use conservation of mass to get average rotor velocity
    avg_vx = vinf * area_ratio

    # loop through RPM's to analyze
    # for i in 1:nops

        # # probe ff_solution velocity field at dtrotor x location
        # velocities = probe_ff_velocity(
        #     ff_solution,
        #     field_points,
        #     vinf[i];
        #     rho=rho,
        #     mu=mu,
        #     treat_singularity=true,
        #     core_size=dtrotor.Rtip * 0.1,
        # )

        # Put velocities in ccblade reference frame
        vx, vy = ff2ccb_velocity(dtrotor, avg_vx, rdist, 1)

        #save tip velocity for mach calculation
        v_tip[1] = sqrt(vx[end]^2 + vy[end]^2)

        # define operating point
        # note: for single open rotor, v_x = vinf and v_y = omega*r
        pitch = 0.0
        rho = 1.225
        mu = 1.81e-5
        asound = 341.0
        ccbops =
            ccb.OperatingPoint.(
                vx, vy, rho, pitch, mu, asound
            )
    # end

    if debug
        #if debugging, return ccblade defaults (tip corrections, simple ops)
        return ccb.Rotor(rhub, rtip, dtrotor.numblades; turbine=false),#, tip=CCBlade.DuctTip()
        ccbsections,
        ccbops,
        v_tip,
        avg_vx
    else
        return ccbrotor, ccbsections, ccbops, v_tip, avg_vx
    end
    # fyi this is how to call CCBlade
    # ccb_out = ccb.solve.(ref(ccbrotor), sections, op)
end

"""
    extract_ff_geom(ff_solution)

Extract x and r coordiantes of duct and hub geometries from the FLOWFoil solution object.

**arguments:**
- `ff_solution::FLOWFoil.InviscidSolution` : inviscid solution object from FLOWFoil.

**returns:**
- `duct_x::Array{float}` : Array of x-coordinates of duct wall control points.
- `duct_r::Array{float}` : Array of r-coordinates of duct wall control points.
- `hub_x::Array{Float}` : Array of x-coordinates of hub wall control points.
- `hub_r::Array{Float}` : Array of r-coordinates of hub wall control points.
"""
function extract_ff_geom(ff_solution)

    # rename for convenience
    mesh1 = ff_solution.meshes[1]
    mesh2 = ff_solution.meshes[2]

    if mesh1.panels[1].controlpoint[2] > mesh2.panels[1].controlpoint[2]

        # duct is mesh1 since TE point is above hub LE point
        cpduct = (panel -> panel.controlpoint).(mesh1.panels)
        cphub = (panel -> panel.controlpoint).(mesh2.panels)

    else

        # otherwise the duct is the second mesh
        cpduct = (panel -> panel.controlpoint).(mesh2.panels)
        cphub = (panel -> panel.controlpoint).(mesh1.panels)
    end

    return getindex.(cpduct, 1),
    getindex.(cpduct, 2), getindex.(cphub, 1),
    getindex.(cphub, 2)
end

"""
    generate_ccb_sections(dtrotor, dtrotor.radialstations)

Generates CCBlade section objects from dtrotor information and dimensional radial stations.

**Arguments:**
- `dtrotor::RotorGeometry` : dtrotor geometry object
- `dtrotor.radialstations::Array{Float}` : Dimensional Radial Station locations

**Returns:**
- `ccbsections::Array{CCBlade.Section}` : Array of CCBlade section objects
- `rdist::Array{float}` : Array of refined, dimensional rotot radial station r locations.
"""
function generate_ccb_sections(dtrotor)

    # rename for convenience
    chords = dtrotor.chords * dtrotor.radialstations[end] # dimensionalize

    twists = dtrotor.twists * pi / 180.0 # convert to radians

    airfoils = dtrotor.airfoils

    # -- Refine Blade
    # spline chord and twist
    csp = FLOWMath.Akima(dtrotor.radialstations, chords)
    tsp = FLOWMath.Akima(dtrotor.radialstations, twists)

    # get dimensiona, refined radial stations
    Rhub = dtrotor.radialstations[1]
    Rtip = dtrotor.radialstations[end]
    rdist = range(Rhub, Rtip; length=dtrotor.nref)

    # get refined chord and twist
    cdist = csp.(rdist)
    tdist = tsp.(rdist)

    # -- assign airfoils
    # initialize airfoil array
    if length(airfoils) > 1
        afdist = Array{typeof(airfoils[1])}(undef, dtrotor.nref)

        # get average values of radial stations to better place airfoils (center defined airfoils in refined blade)
        mean_rad_stash =
            (dtrotor.radialstations[1:(end - 1)] .+ dtrotor.radialstations[2:end]) ./ 2.0

        # loop through refined radial stations and apply appropriate airfoil
        for i in 1:(dtrotor.nref)
            ridx = findfirst(x -> x > rdist[i], mean_rad_stash)
            if ridx != nothing
                afdist[i] = airfoils[ridx]
            else
                afdist[i] = airfoils[end]
            end
        end

        #Return CCBlade sections
        return ccb.Section.(rdist, cdist, tdist, afdist), rdist
    else
        return ccb.Section.(rdist, cdist, tdist, Ref(airfoils[1])), rdist
    end
end

"""
    probe_velocity_axisym(ff_solution, field_points, vinf)

Probe the velocity field for the axisymmetric solution at the given field points.

**Arguements:**
- `ff_solution::FLOWFoil.InviscidSolution` : Inviscid Solution for the axisymmetric problem
- `field_points::Array{Array{Float}}` : Array of field point location arrays.
- `vinf::Float` : Velocity

**Keyword Arguments:**
- `rho::Float` : air density, default: 1.225 kg/m3
- `mu::Float` : air viscosity, default: 1.81e-5 Pa-s
- `treat_singularity::Bool` : flag whether to add treatment for near-wall singularity behavior
- `core_size::Float` : scale factor for weibull function ends up being roughly the radius of a "core size."

**Returns:**
- `velocities::Array{Array{Float}}` : Array of velocities, [u;v], at each field point.
"""
function probe_ff_velocity(
    ff_solution,
    field_points,
    vinf;
    rho=1.225,
    mu=1.81e-5,
    treat_singularity=false,
    core_size=0.025,
)
    T = eltype(vinf)
    #initialize output velocities
    # velocities = [[0.0; 0.0] for i in 1:length(field_points)]
    velocities = [[T(0.0); T(0.0)] for i in 1:length(field_points)]
    if treat_singularity
        edge_velocities = [[T(0.0); T(0.0)] for i in 1:2]
    end

    #loop through each mesh
    for i in 1:length(ff_solution.meshes)

        #get gammas specific to this mesh
        gammas = get_mesh_gammas(ff_solution.panelgammas, ff_solution.meshes, i)

        # loop through panels for this mesh
        for j in 1:length(ff_solution.meshes[i].panels)

            #get current panel
            panel = ff_solution.meshes[i].panels[j]

            #loop through field points
            for k in 1:length(field_points)
                #get relative geometries needed for velocity calculation
                x, r, cpr, dmagj, m = FLOWFoil.get_relative_geometry_axisym(
                    panel, field_points[k]
                )

                ujk = FLOWFoil.get_u_ring(x, r, cpr, dmagj, m; probe=true)
                vjk = FLOWFoil.get_v_ring(x, r, cpr, m; probe=true)

                #add to overall velocity at field point
                velocities[k][1] += ujk * abs(gammas[j]) * dmagj * vinf
                velocities[k][2] += vjk * abs(gammas[j]) * dmagj * vinf
            end

            if treat_singularity

                #get edge of field of affect
                ductcore = [field_points[end][1]; field_points[end][2] - core_size]
                hubcore = [field_points[1][1]; field_points[1][2] + core_size]

                #DUCT
                x, r, cpr, dmagj, m = FLOWFoil.get_relative_geometry_axisym(panel, ductcore)

                ujk = FLOWFoil.get_u_ring(x, r, cpr, dmagj, m; probe=true)
                vjk = FLOWFoil.get_v_ring(x, r, cpr, m; probe=true)

                #add to overall velocity at field point
                edge_velocities[1][1] += ujk * abs(gammas[j]) * dmagj * vinf
                edge_velocities[1][2] += vjk * abs(gammas[j]) * dmagj * vinf

                #HUB
                x, r, cpr, dmagj, m = FLOWFoil.get_relative_geometry_axisym(panel, hubcore)

                ujk = FLOWFoil.get_u_ring(x, r, cpr, dmagj, m; probe=true)
                vjk = FLOWFoil.get_v_ring(x, r, cpr, m; probe=true)

                #add to overall velocity at field point
                edge_velocities[2][1] += ujk * abs(gammas[j]) * dmagj * vinf
                edge_velocities[2][2] += vjk * abs(gammas[j]) * dmagj * vinf
            end
        end
    end

    #get absolute velocity (not just induced)
    for i in 1:length(velocities)
        velocities[i][1] += vinf
    end

    if treat_singularity
        edge_velocities[1][1] += vinf
        edge_velocities[2][1] += vinf

        for k in 1:length(field_points)
            #Duct
            if field_points[end][2] - field_points[k][2] < core_size
                R = (field_points[end][2] - field_points[k][2]) / core_size
                velocities[k][1] = weibull(
                    abs(R); lambda=core_size / 2.0, k=1.0, scale=edge_velocities[1][1]
                )
                velocities[k][2] = weibull(
                    abs(R); lambda=core_size / 2.0, k=1.0, scale=edge_velocities[1][2]
                )
            elseif field_points[k][2] - field_points[1][2] < core_size
                R = (field_points[k][2] - field_points[1][2]) / core_size
                velocities[k][1] = weibull(
                    abs(R); lambda=core_size / 2.0, k=1.0, scale=edge_velocities[2][1]
                )
                velocities[k][2] = weibull(
                    abs(R); lambda=core_size / 2.0, k=1.0, scale=edge_velocities[2][2]
                )
            end
        end
    end

    return velocities
end

"""
    get_mesh_gammas(gammas, meshes, meshidx)

Get the gamma values only for the mesh at index meshidx in meshes.

**Arguments:**
- `gammas::FLOWFoil.InviscidSolution.panelgammas` : vortex strengths at each panel in the system.
- `meshes::Array{FLOWFoil.AxiSymMesh}` : Array of meshes in system
- `meshidx::Int` : index of which mesh in the meshes array for which to obtain the associated gammas.

**Returns:**
- `mesh_gammas::Array{Float}` : panel gamma values for input mesh
"""
function get_mesh_gammas(gammas, meshes, meshidx)

    #initialize offset
    offset = 0

    #if we're interested in values on mesh greater than 1, add to offset
    if meshidx > 1
        for i in 1:(meshidx - 1)
            offset += length(meshes[i].panels)
        end
    end

    #grab the gammas for just the body we want.
    mesh_gammas = gammas[(1 + offset):(offset + length(meshes[meshidx].panels))]

    return mesh_gammas
end

"""
    weibull(x; lambda, k, scale, offset)

Produce a weibull distribution specifically for driving near-wall velocities to zero.

**Arguments:**
- `x::Float` : distance from wall

**Keyword Arguments:**
- `lambda::Float` : parameter of denominator in exponential, roughly associated with halfway point of distribution.
- `k::Float` : parameter of power in exponential, probably want to keep at default, 1.0.
- `scale::Float` : amount to scale distribution output by (velocity where you want to start driving things to zero)
- `offset::Float` : value to offset the output by (if you aren't actually driving things to zero, but some other value)
"""
function weibull(x; lambda=1.0, k=1.0, scale=1.0, offset=0.0)

    #weibull distribution
    w = 1.0 - exp(-(x / lambda)^k)

    #scale into the range from the offset to the max value
    w *= (scale - offset)

    #shift by the offset
    w += offset

    return w
end

"""
    ff2ccb_velocity(rotor, dtrotor.radialstations, velocities, rdist)

Convert velocities from FLOWFoil to the CCBlade reference frame.

**Arguments:**
- `rotor::RotorGeometry` : Rotor geometry object
- `dtrotor.radialstations::Array{Float}` : Dimensional Radial Station locations
- `velocities::Array{Array{Float}}` : Array of [x; r] velocities from FLOWFoil.
- `rdist::Array{float}` : Array of refined, dimensional rotot radial station r locations.

**Returns:**
- `Vx::Array{Float}` : x-velocities in the CCBlade reference frame
- `Vy::Array{Float}` : y-velocities in the CCBlade reference frame
"""
function ff2ccb_velocity(dtrotor, velocities, rdist, jidx)
    # spline velocities
    # vxsp = FLOWMath.Akima(dtrotor.radialstations, getindex.(velocities, 1))
    # vrsp = FLOWMath.Akima(dtrotor.radialstations, getindex.(velocities, 2))

    #use simple conservation of mass value for now
    vxsp = velocities

    vtsp = FLOWMath.Akima(
        dtrotor.radialstations, get_omega(dtrotor.RPM[jidx]) .* dtrotor.radialstations
    )

    # get refined velocities
    Vx = [vxsp for i in 1:length(rdist)] #.(rdist)
    Vy = vtsp.(rdist)

    Vx, Vy = promote(Vx, Vy)

    return Vx, Vy
end
