#=

Funtions related to generation of the wake "grid"

Authors: Judd Mehr,

=#

"""
    initilaize_grid_points(body_geometry, blade_elements, grid_options, debug=false)

Get grid boundary and initial interior points.

**Arguments:**
 - `body_geometry::DuctTAPE.body_geometry` : Duct Geometry Object.
 - `blade_elements::Array{DuctTAPE.Rotor}` : Array of Rotor objects
 - `grid_options::DuctTAPE.GridOptions` : GridOptions object

**Returns:**
 - `x_grid_points::Matrix{Float64,2}` : 2D Array of x grid points
 - `r_grid_points::Matrix{Float64,2}` : 2D Array of r grid points
"""
function initilaize_grid_points(body_geometry, blade_elements; wake_length=1.0, debug=false)

    #check that rotor/stator are in the correct order
    for i in 1:(length(blade_elements) - 1)
        @assert blade_elements[i].rotor_x_location < blade_elements[i + 1].rotor_x_location
    end

    # --- RENAME THINGS FOR CONVENIENCE
    TF = eltype(blade_elements[1].chords)

    duct_range = body_geometry.duct_range
    hub_range = body_geometry.hub_range
    hub_spline = body_geometry.hub_spline
    duct_inner_spline = body_geometry.wallinnerspline

    # Rename options for convenience
    nr = blade_elements[1].num_radial_stations
    # check that both rotors have the same number of elements
    for i in 1:length(blade_elements)
        @assert nr == length(blade_elements[i].radialstations)
    end

    # --- SETUP/FIND BOUNDARY GEOMETRY

    # -- Get rear boundary
    # calculate duct trailing edge
    duct_te = max(duct_range[2], hub_range[2])
    # calculate duct chord
    duct_chord = duct_te - min(duct_range[1], hub_range[1])

    #wake starts at body_geometry trailing edge and extends the input wake_length relative to body_geometry chord.
    wake_te = duct_te * (1.0 + wake_length)

    # --- DEFINE X GRID SPACING PIECES

    # choose x-step to be same length as radial step length, assumes that radial positions are linearly spaced
    x_step = blade_elements[1].radial_positions[2] - radial_positions[1]

    #need to save the x-index of the rotor positions for later use
    rotoridxs = [0 for i in 1:length(blade_elements)]

    # if there are more than 1 rotors, need to account for following rotor locations so that grid aligns with rotors
    if length(blade_elements) > 1

        # find the lengths between rotors using an approximate step size of the blade element spacing
        x_length = zeros(Int, length(blade_elements) - 1)
        for i in 1:(length(blade_elements) - 1)
            x_length[i] = round(
                Int,
                (
                    blade_elements[i + 1].rotor_x_location -
                    blade_elements[i].rotor_x_location
                ) / x_step,
            )
        end

        # save rotor x-index
        rotoridxs = cumsum(x_length)

        #get the length from the last rotor to the end of the wake
        last_length = round(Int, (wake_te - blade_elements[end].rotor_x_location) / x_step)

        # put all the x stations together
        wake_x_stations = reduce(
            vcat,
            [
                [
                    range(
                        blade_elements[i].rotor_x_location,
                        blade_elements[i + 1].rotor_x_location;
                        length=x_length[i],
                    ) for i in length(blade_elements) - 1
                ]
                range(blade_elements[end].rotor_x_location, wake_te; length=last_length)
            ],
        )
    else

        # rotor index is only the first index if there is only one rotor.
        rotoridxs = [1]

        # - If only one rotor, no need for complicated stuff.
        wake_x_stations = range(
            blade_elements[1].rotor_x_location, wake_te; length=x_length
        )
    end

    #get number of x stations
    nx = length(wake_x_stations)

    # --- DEFINE R STATION LOCATIONS ALONG INNER/OUTER BOUNDARIES
    # -- Get radial locations of outer boundary
    router = zeros(TF, nx)
    for i in 1:nx
        if wake_x_stations[i] < duct_range[2]
            router[i] = duct_inner_spline(wake_x_stations[i])
        else
            router[i] = duct_inner_spline(duct_range[2])
        end
    end

    # -- Get radial locations of inner boundary
    rinner = zeros(TF, nx)
    for i in 1:nx
        if wake_x_stations[i] < hub_range[2]
            rinner[i] = hub_spline(wake_x_stations[i])
        else
            rinner[i] = hub_spline(hub_range[2])
        end
    end

    # --- INITIALIZE GRID POINTS
    #initialize grid matrix
    x_grid_points = zeros(TF, nx, nr)
    r_grid_points = zeros(TF, nx, nr)

    #first column of radial grid points are radial stations of foremost rotor
    r_grid_points[1, :] = blade_elements[1].radial_positions

    #first column of x grid points are all the x position of the foremost rotor
    x_grid_points[1, :] .= blade_elements[1].rotor_x_location

    #first and last rows of points from setup above
    x_grid_points[:, 1] = wake_x_stations
    x_grid_points[:, end] = wake_x_stations
    r_grid_points[:, 1] = rinner
    r_grid_points[:, end] = router

    # --- DEFINE INNER GRID POINTS

    # -- Set up grid definition Loop
    #Get area of foremost rotor annulus
    rotorarea = pi * (router[1]^2 - rinner[1]^2)

    #get the annulus area of the radial stations
    drotorarea = [
        pi * (r_grid_points[1, j + 1]^2 - r_grid_points[1, j]^2) for j in 1:(nr - 1)
    ]

    #get areas of subsequent annuli
    localarea = [pi * (r_grid_points[i, end]^2 - r_grid_points[i, 1]^2) for i in 1:nx]

    # March conservation of mass along x-direction to get grid points.
    for i in 2:nx #for each x station after the rotor plane

        #get the x,r coordinates at the x stations
        xwall = x_grid_points[i, end]
        xhub = x_grid_points[i, 1]
        rwall = r_grid_points[i, end]
        rhub = r_grid_points[i, 1]

        # get the area of the annulus at the current x station
        # localarea = pi * (rwall^2 - rhub^2)

        for j in 2:(nr - 1) #for each radial station between the hub and wall

            #get local change in area using conservation of mass
            dlocalarea = (drotorarea[j - 1] * localarea[i]) / rotorarea

            #define radial grid point based on conservation of mass
            r_grid_points[i, j] = sqrt(dlocalarea / pi + r_grid_points[i, j - 1]^2)

            #NOTE that the rest of the loop only does something if xs aren't the same like they are now.
            #calculate ratio for weighted average
            frac = (r_grid_points[i, j] - rhub) / (rwall - rhub)

            # place x grid point based on weighted ratio
            x_grid_points[i, j] = xhub + frac * (xwall - xhub)
        end # for radial stations
    end # for x stations

    return x_grid_points, r_grid_points, nx, nr, rotoridxs
end
