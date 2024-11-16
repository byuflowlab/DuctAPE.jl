"""
    ProblemDimensions{Int}

Struct containing dimensions of the problem used throughout the analysis.

- `nrotor`    : number of rotors
-  nwn`       : number of wake nodes
-  nwp`       : number of wake panels
-  ncp`       : number of casing panels
-  ndn`       : number of duct nodes
-  ncbn`      : number of centerbody nodes
-  nbn`       : number of body nodes
-  nbp`       : number of body panels
-  nws`       : number of wake sheets (also rotor nodes)
-  nbe`       : number of blade elements (also rotor panels)
-  nwsn`      : number of nodes in each wake sheet
-  nwsp`      : number of panels in each wake sheet
-  ndwin`     : number of duct-wake interfacing nodes
-  ncbwin`    : number of centerbody-wake interfacing nodes
-  nbodies=2` : number of bodies (currently hardcoded to 2)
"""
@kwdef struct ProblemDimensions{TI}
    nrotor::TI      # number of rotors
    nwn::TI         # number of wake nodes
    nwp::TI         # number of wake panels
    ncp::TI         # number of casing panels
    ndn::TI         # number of duct nodes
    ncbn::TI        # number of centerbody nodes
    nbn::TI         # number of body nodes
    nbp::TI         # number of body panels
    nws::TI         # number of wake sheets (also rotor nodes)
    nbe::TI         # number of blade elements (also rotor panels)
    nwsn::TI        # number of nodes in each wake sheet
    nwsp::TI        # number of panels in each wake sheet
    ndwin::TI       # number of duct-wake interfacing nodes
    ncbwin::TI      # number of centerbody-wake interfacing nodes
    nbodies::TI = 2 # hard code this for now.
end

"""
    get_problem_dimensions(paneling_constants::PanelingConstants)
    get_problem_dimensions(body_vortex_panels, rotor_source_panels, wake_vortex_panels)

Determine all relevant dimensions to the problem based either on the paneling_constants or the panels themselves.

# Arguments
- `paneling_constants::PanelingConstants` : Rotor (and possibly stator) geometric paramters.

# Returns
- `problem_dimensions::ProblemDimensions` : ProblemDimensions object.
"""
function get_problem_dimensions(paneling_constants::PanelingConstants)

    # - Extract Paneling Constants - #
    (; npanels, ncenterbody_inlet, nduct_inlet, wake_length, nwake_sheets, dte_minus_cbte) =
        paneling_constants

    # number of rotors is one less than the length of npanels if the duct and hub trailing edges line up, and is two less if they don't
    nrotor = iszero(dte_minus_cbte) ? length(npanels) - 1 : length(npanels) - 2

    # number of wake sheets (blade nodes)
    nws = nwake_sheets

    # number of blade elements (panels)
    nbe = nws - 1

    # number of body panels
    ncp = nduct_inlet
    ncbp = ncenterbody_inlet
    # add rest of panels mutual between centerbody and duct
    if iszero(dte_minus_cbte)
        ncp += sum(npanels[1:(end - 1)])
        ncbp += sum(npanels[1:(end - 1)])
    else
        ncp += sum(npanels[1:(end - 2)])
        ncbp += sum(npanels[1:(end - 2)])
    end

    # add additional duct or centerbody panels if one extends further back
    if dte_minus_cbte > 0
        ncp += npanels[end - 1]
    elseif dte_minus_cbte < 0
        ncbp += npanels[end - 1]
    end

    # duct panels are 2x the number of nacelle panels
    ndp = 2 * ncp

    # duct and center body nodes are 1 more than number of panels
    ndn = ndp + 1
    ncbn = ncbp + 1

    # number of body panels is sum of duct and centerbody panels
    nbp = ndp + ncbp
    # number of body nodes is number of panels + number of bodies
    nbn = ndn + ncbn

    # number of panels in each wake sheet
    nwsp = sum(npanels)
    # number of nodes in each wake sheet
    nwsn = nwsp + 1

    # number of wake panels is the total number of npanels times the number of wake sheets
    nwp = sum(npanels) * nwake_sheets

    # number of wake nodes is one more than the number of panels for each wake sheet
    nwn = nwp + nwake_sheets

    # number of duct-wake and centerbody-wake interface nodes
    if iszero(dte_minus_cbte)
        ndwin = sum(npanels[1:(end - 1)]) + 1
        ncbwin = sum(npanels[end - 1]) + 1
    elseif dte_minus_cbte < 0
        ndwin = sum(npanels[1:(end - 2)]) + 1
        ncbwin = sum(npanels[end - 1]) + 1
    else
        ndwin = sum(npanels[1:(end - 1)]) + 1
        ncbwin = sum(npanels[end - 2]) + 1
    end

    return ProblemDimensions(;
        nrotor,     # number of rotors
        nwn,    # number of wake nodes
        nwp,    # number of wake panels
        ncp,    # number of casing panels
        ndn,    # number of duct nodes
        ncbn,   # number of centerbody nodes
        nbn,    # number of body nodes
        nbp,    # number of body panels
        nws,    # number of wake sheets (also rotor nodes)
        nbe,    # number of blade elements (also rotor panels)
        nwsn,   # number of nodes in each wake sheet
        nwsp,   # number of panels in each wake sheet
        ndwin,  # number of duct-wake interfacing nodes
        ncbwin, # number of centerbody-wake interfacing nodes
        nbodies=2, #hard code this for now.
    )
end

"""
"""
function get_problem_dimensions(body_vortex_panels, rotor_source_panels, wake_vortex_panels)
    # number of rotors
    nrotor = rotor_source_panels.nbodies[]

    # number of wake sheets (blade nodes)
    nws = wake_vortex_panels.nbodies[]
    # number of blade elements (panels)
    nbe = nws - 1

    # number of body panels
    ncp =
        findmin(@view(body_vortex_panels.node[1, 1:Int(body_vortex_panels.nnode[1])]))[2] -
        1 #TODO check this is correct
    ndp = body_vortex_panels.npanel[1]
    ncbp = body_vortex_panels.npanel[2]

    # duct and center body nodes
    ndn = body_vortex_panels.nnode[1]
    ncbn = body_vortex_panels.nnode[2]

    # number of body panels
    nbp = body_vortex_panels.totpanel[]
    # number of body nodes
    nbn = body_vortex_panels.totnode[]

    # number of panels in each wake sheet
    nwsp = wake_vortex_panels.npanel[1]
    # number of nodes in each wake sheet
    nwsn = wake_vortex_panels.nnode[1]

    # number of wake panels total
    nwp = wake_vortex_panels.totpanel[]
    # number of wake nodes total
    nwn = wake_vortex_panels.totnode[]

    # number of duct-wake and centerbody-wake interface nodes
    ndwin = length(
        intersect(
            findall(
                in(
                    body_vortex_panels.node[
                        1,
                        Int(body_vortex_panels.endnodeidxs[1, 1]):Int(
                            body_vortex_panels.endnodeidxs[2, 1]
                        ),
                    ],
                ),
                wake_vortex_panels.node[
                    1,
                    Int(wake_vortex_panels.endnodeidxs[1, end]):Int(
                        wake_vortex_panels.endnodeidxs[2, end]
                    ),
                ],
            ),
            findall(
                in(
                    body_vortex_panels.node[
                        2,
                        Int(body_vortex_panels.endnodeidxs[1, 1]):Int(
                            body_vortex_panels.endnodeidxs[2, 1]
                        ),
                    ],
                ),
                wake_vortex_panels.node[
                    2,
                    Int(wake_vortex_panels.endnodeidxs[1, end]):Int(
                        wake_vortex_panels.endnodeidxs[2, end]
                    ),
                ],
            ),
        ),
    )

    ncbwin = length(
        intersect(
            findall(
                in(
                    body_vortex_panels.node[
                        1,
                        Int(body_vortex_panels.endnodeidxs[1, 2]):Int(
                            body_vortex_panels.endnodeidxs[2, 2]
                        ),
                    ],
                ),
                wake_vortex_panels.node[
                    1,
                    Int(wake_vortex_panels.endnodeidxs[1, 1]):Int(
                        wake_vortex_panels.endnodeidxs[2, 1]
                    ),
                ],
            ),
            findall(
                in(
                    body_vortex_panels.node[
                        2,
                        Int(body_vortex_panels.endnodeidxs[1, 2]):Int(
                            body_vortex_panels.endnodeidxs[2, 2]
                        ),
                    ],
                ),
                wake_vortex_panels.node[
                    2,
                    Int(wake_vortex_panels.endnodeidxs[1, 1]):Int(
                        wake_vortex_panels.endnodeidxs[2, 1]
                    ),
                ],
            ),
        ),
    )

    return ProblemDimensions(;
        nrotor, # number of rotors
        nwn,    # number of wake nodes
        nwp,    # number of wake panels
        ncp,    # number of casing panels
        ndn,    # number of duct nodes
        ncbn,   # number of centerbody nodes
        nbn,    # number of body nodes
        nbp,    # number of body panels
        nws,    # number of wake sheets (also rotor nodes)
        nbe,    # number of blade elements (also rotor panels)
        nwsn,   # number of nodes in each wake sheet
        nwsp,   # number of panels in each wake sheet
        ndwin,  # number of duct-wake interfacing nodes
        ncbwin, # number of centerbody-wake interfacing nodes
        nbodies=2, #hardcode for now.
    )
end
