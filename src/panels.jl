#=
Functions for paneling system

Authors: Judd Mehr,
=#

###############################
##### ----- EXPORTS ----- #####
###############################

## -- TYPES

export PanelSystem, Panels

## -- FUNCTIONS

export generate_paneling, generate_panel_system

"""
    Panels{TPEx,TPEr,TPC,TPN,TPT}

**Fields:**
 - `panel_edges_x::Array{Tuple{Float, Float}}` : Array of sets of x locations for panel edges
 - `panel_edges_r::Array{Tuple{Float, Float}}` : Array of sets of r locations for panel edges
 - `panel_centers::Array{Tuple{Float, Float}}` : Array of sets of x,r locations for panel centers
 - `panel_normals::Array{Tuple{Float, Float}}` : Array of panel unit normal vectors
 - `panel_tyes::Array{String}` : Array of panel types (for use in assembling linear system)
"""
struct Panels{TPEx,TPEr,TPC,TPN,TPT}
    panel_edges_x::TPEx
    panel_edges_r::TPEr
    panel_centers::TPC
    panel_normals::TPN
    panel_types::TPT
end

"""
    PanelSystem{TD,TH,TW,TR}

**Fields:**
 - `wall_panels::DuctTAPE.Panels` : panels defining duct wall airfoil
 - `hub_panels::DuctTAPE.Panels` : panels defining hub
 - `wake_panels::DuctTAPE.Panels` : panels defining rotor wake vortex sheets
 - `rotor_source_panels::DuctTAPE.Panels` : panels defining rotor drag source panels
"""
struct PanelSystem{TD,TH,TW,TR}
    wall_panels::TD
    hub_panels::TH
    wake_panels::TW
    rotor_source_panels::TR
    #drag_panels::TD
end

"""
    generate_paneling(ductgeometry, ductsplines, rotors, wakegrid)

Generate panel edges, centers, and unit normals.

**Arguments:**
 - `ductgeometry::DuctTAPE.DuctGeometry` : Duct Geometry object
 - `ductsplines::DuctTAPE.DuctSplines` : Duct Splines object
 - `rotors::Array{DuctTAPE.Rotor}` : Array of rotor objects
 - `wakegrid::DuctTAPE.WakeGridGeometry` : Wake Grid object

**Returns:**
 - `wall_panels::DuctTAPE.Panels` : Panels object for duct wall
 - `hub_panels::DuctTAPE.Panels` : Panels object for hub
 - `wake_panels::DuctTAPE.Panels` : Panels object for vortex wake sheets
 - `rotor_source_panels::Array{DuctTAPE.Panels}` : Array of Panels objects for each rotor

# NOTES:
- The paneling for the rotor sources and the vortex wake sheets is based directly on the wake grid.
- The paneling of the duct wall and hub are set such that the panels aft of the foremost rotor also align perfectly with the wake grid.
- The wall panels in front of the foremost rotor are set using cosine spacing such that the last panel before the foremost rotor is roughly similar in length to the average of the panel lengths in the remainder of the duct.
(Note that the estimation process for this is not particularly robust at this point.)
"""
function generate_paneling(ductgeometry, ductsplines, rotors, wakegrid)

    ## -- INITIALIZE ARRAYS -- ##
    # i.e. Count number of panels

    # - number of hub geometry panels
    #take number of xcoordinates for the hub and subtract the last one to get the number of panels
    numhubpan = length(wakegrid.hub_xstations) - 1

    #if the hub has a blunt trailing edge. need to add another panel.
    if ductgeometry.hubbluntTE
        numhubpan += 1
    end

    # - number of wall geometry panels
    #take the number of x coordinates for the inner and outer wall geometry and subtract the last point from each to get the number of panels for the whole airfoil
    numwallpan = 2 * length(wakegrid.wall_xstations) - 2

    #add another panel if the duct wall has a blunt trailing edge.
    if ductgeometry.wallbluntTE
        numwallpan += 1
    end

    # - number of rotor source panels
    #rename num_radial_stations for convenience
    numrs = length(rotors[1].radialstations)

    #number of rotor panels is going to be the number of rotors multiplied by one fewer than the number of radial stations
    numrotorsourcepan = numrs - 1

    # - number of wake panels
    # TODO: need to make sure that panels don't double up for multiple rotors
    #each internal (not hub/tip) rotor radial station sheds a vortex sheet that extends to the end of the grid.
    numrotorwakepan = (numrs - 2) * (wakegrid.nx - 1)

    #the hub sheds a vortex sheet from its trailing edge to the end of the grid.
    numhubTEwakepan = wakegrid.nx - wakegrid.hubTEidx

    #the duct wall also sheds a vortex sheet from its trailing edge to the end of the grid.
    numwallTEwakepan = wakegrid.nx - wakegrid.wallTEidx

    # - number of drag object source panels (TODO LATER)

    ## -- GET PANEL EDGES, CENTERS AND TYPES -- ##

    # - wall panels
    #initialize from counts above.
    wall_panel_edge_x = [(0.0, 0.0) for i in 1:numwallpan]
    wall_panel_edge_r = [(0.0, 0.0) for i in 1:numwallpan]
    wall_panel_center = [(0.0, 0.0) for i in 1:numwallpan]
    wall_panel_normal = [(0.0, 0.0) for i in 1:numwallpan]

    # need to combine the inner and outer geometry coordinates
    # TODO: the user input decides what order things are in.  Need to make sure that the splines are defined in the correct directions, so probably need to add some checks at the point of the ductgeometry generation.
    wallinnerxcoordinates = wakegrid.wall_xstations
    wallouterxcoordinates = wakegrid.wall_xstations
    wallxcoordinates = [reverse(wallinnerxcoordinates)[1:(end - 1)]; wallouterxcoordinates]

    #similar for r coordinates
    wallinnerrcoordinates = ductsplines.wallinnerspline(wallinnerxcoordinates)
    wallouterrcoordinates = ductsplines.wallouterspline(wallouterxcoordinates)
    wallrcoordinates = [reverse(wallinnerrcoordinates)[1:(end - 1)]; wallouterrcoordinates]
    #TODO: need to make sure this is the correct direction. (what is the correct direction??)

    # loop through the total number of wall coordinates
    for i in 1:numwallpan
        if ductgeometry.wallbluntTE && i == numwallpan

            #create the trailing edge panel if needed
            wall_panel_edge_x[i] = (wallxcoordinates[end], wallxcoordinates[1])

            wall_panel_edge_r[i] = (wallrcoordinates[end], wallrcoordinates[1])

        else
            #get the x coordinates from the geometry and save as the panel edges.
            wall_panel_edge_x[i] = (wallxcoordinates[i], wallxcoordinates[i + 1])

            #similar for r coordinates
            wall_panel_edge_r[i] = (wallrcoordinates[i], wallrcoordinates[i + 1])
        end

        #no matter what case, get the average of the x and r coordinates to find the center point on the panel.
        wall_panel_center[i] = (
            sum(wall_panel_edge_x[i]) / 2.0, sum(wall_panel_edge_r[i]) / 2.0
        )

        # calculate unit normal vectors on panels
        # #TODO: It's important that the panels are defined in the correct direction, otherwise, these normals will be the wrong directions.  need to figure out which way dfdc wants it.
        wall_panel_normal[i] = calc_normal(wall_panel_edge_x[i], wall_panel_edge_r[i])
    end

    # create wall panels object
    wall_panels = Panels(
        wall_panel_edge_x, wall_panel_edge_r, wall_panel_center, wall_panel_normal, "w"
    )

    # - hub panels
    #initialize from counts above
    hub_panel_edge_x = [(0.0, 0.0) for i in 1:numhubpan]
    hub_panel_edge_r = [(0.0, 0.0) for i in 1:numhubpan]
    hub_panel_center = [(0.0, 0.0) for i in 1:numhubpan]
    hub_panel_normal = [(0.0, 0.0) for i in 1:numhubpan]

    hubxcoordinates = wakegrid.hub_xstations

    #similar for r coordinates
    hubrcoordinates = ductsplines.hubspline(hubxcoordinates)

    #loop through hub panel count
    for i in 1:numhubpan
        if ductgeometry.hubbluntTE && i == numhubpan

            #define trailing edge panel if needed
            hub_panel_edge_x[i] = (hubxcoordinates[end], hubxcoordinates[end])
            hub_panel_edge_r[i] = (hubrcoordinates[end], 0.0)

        else
            # get the x geometry coordinates and set as the panel edges
            hub_panel_edge_x[i] = (hubxcoordinates[i], hubxcoordinates[i + 1])

            #similar for r coordinates
            hub_panel_edge_r[i] = (hubrcoordinates[i], hubrcoordinates[i + 1])
        end

        #no mater the case, get the average of the x and r coordinates to set the center points of the panels
        hub_panel_center[i] = (
            sum(hub_panel_edge_x[i]) / 2.0, sum(hub_panel_edge_r[i]) / 2.0
        )

        #calculate unit normals on panels
        hub_panel_normal[i] = calc_normal(hub_panel_edge_x[i], hub_panel_edge_r[i])
    end

    # create the hub panels object
    hub_panels = Panels(
        hub_panel_edge_x, hub_panel_edge_r, hub_panel_center, hub_panel_normal, "w"
    )

    # - rotor source panels
    rotor_source_panels = Array{Panels}(undef, 2)

    #rename pertinent values for convenience
    rotoridxs = wakegrid.rotoridxs
    gridxs = wakegrid.x_grid_points
    gridrs = wakegrid.r_grid_points

    #loop through each rotor
    for i in 1:length(rotors)

        #initialze based on the counts above
        rotor_panel_edge_x = [(0.0, 0.0) for j in 1:numrotorsourcepan]
        rotor_panel_edge_r = [(0.0, 0.0) for j in 1:numrotorsourcepan]
        rotor_panel_center = [(0.0, 0.0) for j in 1:numrotorsourcepan]
        rotor_panel_normal = [(0.0, 0.0) for j in 1:numrotorsourcepan]

        #create a blade object for each rotor (to get dimensional radial data)
        blade = initialize_blade_dimensions(ductgeometry, ductsplines, rotors[i])

        #loop through each of the radial stations
        for j in 1:numrotorsourcepan

            # for now, the xlocation is the x coordinate of the panels
            #TODO: this will change when rake is added.
            rotor_panel_edge_x[j] = (gridxs[rotoridxs[i], j], gridxs[rotoridxs[i], j + 1])

            #use radial stations to define panel edges.
            rotor_panel_edge_r[j] = (gridrs[rotoridxs[i], j], gridrs[rotoridxs[i], j + 1])

            #as before, average the edges to get the centers
            rotor_panel_center[j] = (
                sum(rotor_panel_edge_x[j]) / 2.0, sum(rotor_panel_edge_r[j]) / 2.0
            )

            #calculate unit normals on panels
            rotor_panel_normal[i] = calc_normal(
                rotor_panel_edge_x[i], rotor_panel_edge_r[i]
            )
        end


        # create the rotor source panels object
        rotor_source_panels[i] = Panels(
            rotor_panel_edge_x,
            rotor_panel_edge_r,
            rotor_panel_center,
            rotor_panel_normal,
            "s",
        )
    end

    # - wake panels

    #add together all the various wake panel counts from above
    numwakepan = numrotorwakepan + numhubTEwakepan + numwallTEwakepan

    #initialize from the total count
    wake_panel_edge_x = [(0.0, 0.0) for i in 1:numwakepan]
    wake_panel_edge_r = [(0.0, 0.0) for i in 1:numwakepan]
    wake_panel_center = [(0.0, 0.0) for i in 1:numwakepan]
    wake_panel_normal = [(0.0, 0.0) for i in 1:numwakepan]

    #loop through the hub trailing edge vortex sheet panels first.
    for i in 1:(numhubTEwakepan)

        #rename for convenience: position starts at the TE index
        hidx = wakegrid.hubTEidx + i - 1

        #use the grid coordinates from the hub TE to the end of the grid at the first radial position of the grid as panel edges.
        wake_panel_edge_x[i] = (
            wakegrid.x_grid_points[hidx, 1], wakegrid.x_grid_points[hidx + 1, 1]
        )

        #same for r positions
        wake_panel_edge_r[i] = (
            wakegrid.r_grid_points[hidx, 1], wakegrid.r_grid_points[hidx + 1, 1]
        )

        #again, average for centers.
        wake_panel_center[i] = (
            sum(wake_panel_edge_x[i]) / 2.0, sum(wake_panel_edge_r[i]) / 2.0
        )

        #calculate unit normals on panels
        wake_panel_normal[i] = calc_normal(wake_panel_edge_x[i], wake_panel_edge_r[i])
    end

    #create a custom index for going forward since there's not an easier way to count.
    wpanidx = numhubTEwakepan + 1

    # loop through the grid radial stations
    for i in 2:(wakegrid.nr - 1)

        #loop through the grid x stations (from first rotor to end of wake)
        for j in 1:(wakegrid.nx - 1)

            # get the panel edges from the grid points
            wake_panel_edge_x[wpanidx] = (
                wakegrid.x_grid_points[j, i], wakegrid.x_grid_points[j + 1, i]
            ) #TODO: check that the indexing here is correct.  is it [r,x] or [x,r]?

            # same for r points
            wake_panel_edge_r[wpanidx] = (
                wakegrid.r_grid_points[j, i], wakegrid.r_grid_points[j + 1, i]
            ) #TODO: check that the indexing here is correct.  is it [r,x] or [x,r]?

            #again, average for centers
            wake_panel_center[wpanidx] = (
                sum(wake_panel_edge_x[wpanidx]) / 2.0, sum(wake_panel_edge_r[wpanidx]) / 2.0
            )

            #calculate unit normals on panels
            wake_panel_normal[wpanidx] = calc_normal(
                wake_panel_edge_x[wpanidx], wake_panel_edge_r[wpanidx]
            )

            #update custom index
            wpanidx += 1
        end
    end

    #loop through wall trailing edge vortex sheet panels
    for i in 1:(numwallTEwakepan)

        #rename index for convenience
        widx = wakegrid.wallTEidx + i - 1

        #get grid points at top radial position starting at duct wall trailing edge index.
        wake_panel_edge_x[wpanidx] = (
            wakegrid.x_grid_points[widx, end], wakegrid.x_grid_points[widx + 1, end]
        )

        #similar for r coordinates
        wake_panel_edge_r[wpanidx] = (
            wakegrid.r_grid_points[widx, end], wakegrid.r_grid_points[widx + 1, end]
        )

        #average for centers
        wake_panel_center[wpanidx] = (
            sum(wake_panel_edge_x[wpanidx]) / 2.0, sum(wake_panel_edge_r[wpanidx]) / 2.0
        )

        #calculate unit normals on panels
        wake_panel_normal[wpanidx] = calc_normal(
            wake_panel_edge_x[wpanidx], wake_panel_edge_r[wpanidx]
        )

        # increment custom index
        wpanidx += 1
    end

    # create wake (vortex sheet) panels object
    wake_panels = Panels(
        wake_panel_edge_x, wake_panel_edge_r, wake_panel_center, wake_panel_normal, "v"
    )

    #return the 3 types of panel objects.
    return wall_panels, hub_panels, wake_panels, rotor_source_panels
end

"""
    generate_panel_system(ductgeometry, ductsplines, rotors, wakegrid)

Put all the various panel objects together for convenience.

**Arguments:**
 - `ductgeometry::DuctTAPE.DuctGeometry` : Duct Geometry object
 - `ductsplines::DuctTAPE.DuctSplines` : Duct Splines object
 - `rotors::Array{DuctTAPE.Rotor}` : Array of rotor objects
 - `wakegrid::DuctTAPE.WakeGridGeometry` : Wake Grid object

**Returns:**
 - `panelsystem::DuctTAPE.PanelSystem` : All System Panels
"""
function generate_panel_system(ductgeometry, ductsplines, rotors, wakegrid)

    #generate individual panel objects
    wall_panels, hub_panels, wake_panels, rotor_source_panels = generate_paneling(
        ductgeometry, ductsplines, rotors, wakegrid
    )

    #return panel system containing all panel objects
    return PanelSystem(wall_panels, hub_panels, wake_panels, rotor_source_panels)
end
