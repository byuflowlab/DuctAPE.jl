function run_body_only(panels, Vinf, prescribedpanels)
    #---------------------------------#
    #       Operating Conditions      #
    #---------------------------------#
    # Define freestream on panels
    # Vinf = 1.0 #magnitude doesn't matter yet.
    Vs = Vinf * [1.0 0.0] # axisymmetric, so no radial component
    Vsmat = repeat(Vs, panels.npanels) # need velocity on each panel

    # prescribe a panel for the least squares solve:
    # choose the first panel to be prescirbed to zero (assumes first panel is not hub leading/traling edge).
    # prescribedpanels = [(1, 0.0)]

    #---------------------------------#
    #        Induced Velocities       #
    #---------------------------------#

    # - Initial System Matrices - #
    LHS = dt.doublet_panel_influence_matrix(panels.nodes, panels)
    RHS = dt.freestream_influence_vector(panels.normal, Vsmat)
    LHSnokutta = deepcopy(LHS)
    RHSnokutta = deepcopy(RHS)

    # - Adding Kutta Condition - #
    dt.body_lhs_kutta!(LHS, panels; tol=1e1 * eps(), verbose=true)

    # - Prepping for Least Sqaures Solve - #
    # without kutta condition
    LHSlsq_nokutta, RHSlsq_nokutta = dt.prep_leastsquares(
        LHSnokutta, RHSnokutta, prescribedpanels
    )
    # with kutta condition
    LHSlsq, RHSlsq = dt.prep_leastsquares(LHS, RHS, prescribedpanels)

    #---------------------------------#
    #             Solving             #
    #---------------------------------#

    mured = LHSlsq \ RHSlsq
    mured_nokutta = LHSlsq_nokutta \ RHSlsq_nokutta

    # - Solving Without Kutta - #
    # mu_nokutta = LHSlsq_nokutta\RHSlsq_nokutta
    mu_nokutta = dt.mured2mu(mured_nokutta, prescribedpanels)

    # - Solving With Kutta - #
    mu = dt.mured2mu(mured, prescribedpanels)

    #---------------------------------#
    #         Post-Processing         #
    #---------------------------------#

    ### --- Velocity Contributions --- ###
    # - Body-induced Surface Velocity - #
    Vb_nokutta = dt.vfromdoubletpanels(panels.controlpoint, panels.nodes, mu_nokutta)
    Vb = dt.vfromdoubletpanels(panels.controlpoint, panels.nodes, mu)

    # - "Wake"-induced Surface Velocity - #
    Vb_nokutta_wake = dt.vfromTE(
        panels.controlpoint, panels.endpoints, panels.endpointidxs, mu_nokutta
    )
    Vb_wake = dt.vfromTE(panels.controlpoint, panels.endpoints, panels.endpointidxs, mu)

    # - ∇μ/2 surface velocity - #
    Vb_gradmu = dt.vfromgradmu(panels, mu)

    # - Total Velocity - #
    V_nokutta = Vb_nokutta .+ Vb_nokutta_wake
    V_nogradmu = Vb .+ Vb_wake
    Vtot = Vsmat .+ Vb .+ Vb_wake .+ Vb_gradmu

    # ### --- Velocity Tangent to Surface --- ###
    # Vtan_nokutta = [dt.dot(v,t) for (v,t) in zip(eachrow(V_nokutta), eachrow(panels.tangent))]
    # Vtan_nogradmu = [dt.dot(v,t) for (v,t) in zip(eachrow(V_nogradmu), eachrow(panels.tangent))]
    # Vtan = [dt.dot(v,t) for (v,t) in zip(eachrow(Vtot), eachrow(panels.tangent))]

    ### --- Steady Surface Pressure --- ###
    cp_nokutta = 1.0 .- (dt.norm.(eachrow(V_nokutta)) / Vinf) .^ 2
    cp_nogradmu = 1.0 .- (dt.norm.(eachrow(V_nogradmu)) / Vinf) .^ 2
    cp = 1.0 .- (dt.norm.(eachrow(Vtot)) / Vinf) .^ 2

    return mu
end
