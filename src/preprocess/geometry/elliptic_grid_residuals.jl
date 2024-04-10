function elliptic_grid_residual!(r, y, x, p)

    # - extract parameters - #
    (; xi, eta, gridshape, itshape, proposed_grid_cache) = p
    proposed_grid = get_tmp(proposed_grid_cache, y)
    proposed_grid .= reshape(x, gridshape)

    # dimensions
    nxi = length(xi)
    neta = length(eta)

    # - Initialize Residual - #
    # note: the boundary values are kept constant, so want to set those to zero before getting started
    r .= y .- reshape(proposed_grid[:, 2:end, 2:(end - 1)], :)

    # - reshape vector into 3D array - #
    resid = reshape(r, itshape)

    #overwrite proposed_grid internals with the current state variables
    proposed_grid[:, 2:end, 2:(end - 1)] .= reshape(y, itshape)

    # separate x and r components of proposed wake_grid
    xr = view(proposed_grid, 1, :, :)
    rr = view(proposed_grid, 2, :, :)

    # for cache testing
    # TF = eltype(y)
    # x_xi = TF[1.0]

    # loop over interior streamlines
    for j in 2:(neta - 1)

        # --- Outlet Residuals --- #

        # Neumann Boundary Condition (u^2 + v^2 = constant)

        # eta unit vector components in x and r directions
        dx = xr[end, j + 1] - xr[end, j - 1]
        dr = rr[end, j + 1] - rr[end, j - 1]

        enorm = sqrt(dx^2 + dr^2)
        xhat = iszero(enorm) ? 0.0 : dx / enorm
        rhat = iszero(enorm) ? 0.0 : dr / enorm

        # calculate 1st derivatives using finite differencing on irregular wake_grid
        dxi1 = xi[end] - xi[end - 1]
        dxi2 = xi[end - 1] - xi[end - 2]

        dxdx1 = (xr[end, j] - xr[end - 1, j]) / dxi1
        dxdx2 = (xr[end - 1, j] - xr[end - 2, j]) / dxi2

        drdx1 = (rr[end, j] - rr[end - 1, j]) / dxi1
        drdx2 = (rr[end - 1, j] - rr[end - 2, j]) / dxi2

        x_xi = dxdx1 + dxi1 * (dxdx1 - dxdx2) / (dxi1 + dxi2)
        r_xi = drdx1 + dxi1 * (drdx1 - drdx2) / (dxi1 + dxi2)

        # populate residual vector
        tmp = xhat * x_xi + rhat * r_xi
        resid[1, end, j - 1] = xhat * tmp
        resid[2, end, j - 1] = rhat * tmp

        # -- Interior Residuals --- #

        for i in 2:(nxi - 1)

            # calculate 1st derivatives using finite differencing on an irregular wake_grid
            dximinus = xi[i] - xi[i - 1]
            dxiplus = xi[i + 1] - xi[i]
            dxiavg = 0.5 * (dximinus + dxiplus)

            detaminus = eta[j] - eta[j - 1]
            detaplus = eta[j + 1] - eta[j]
            detaavg = 0.5 * (detaminus + detaplus)

            x_eta = (xr[i, j + 1] - xr[i, j - 1]) / (2.0 * detaavg)
            r_eta = (rr[i, j + 1] - rr[i, j - 1]) / (2.0 * detaavg)
            x_xi = (xr[i + 1, j] - xr[i - 1, j]) / (2.0 * dxiavg)
            r_xi = (rr[i + 1, j] - rr[i - 1, j]) / (2.0 * dxiavg)

            # calculate 2nd derivatives using finite differencing on an irregular wake_grid
            ravgplus = 0.5 * (rr[i, j] + rr[i, j + 1])
            ravgminus = 0.5 * (rr[i, j] + rr[i, j - 1])
            ravg = 0.5 * (ravgplus + ravgminus)

            ximinuscoeff = detaminus * detaplus / (dximinus * dxiavg)
            xipluscoeff = detaminus * detaplus / (dxiplus * dxiavg)
            etaminuscoeff = detaplus / detaavg * ravgminus / ravg
            etapluscoeff = detaminus / detaavg * ravgplus / ravg

            x_xixi =
                (xr[i + 1, j] - xr[i, j]) * xipluscoeff -
                (xr[i, j] - xr[i - 1, j]) * ximinuscoeff
            x_etaeta =
                (xr[i, j + 1] - xr[i, j]) * etapluscoeff -
                (xr[i, j] - xr[i, j - 1]) * etaminuscoeff
            x_xieta =
                detaminus *
                detaplus *
                (xr[i + 1, j + 1] - xr[i - 1, j + 1] - xr[i + 1, j - 1] + xr[i - 1, j - 1]) /
                (4.0 * dxiavg * detaavg)
            r_xixi =
                (rr[i + 1, j] - rr[i, j]) * xipluscoeff -
                (rr[i, j] - rr[i - 1, j]) * ximinuscoeff
            r_etaeta =
                (rr[i, j + 1] - rr[i, j]) * etapluscoeff -
                (rr[i, j] - rr[i, j - 1]) * etaminuscoeff
            r_xieta =
                detaminus *
                detaplus *
                (rr[i + 1, j + 1] - rr[i - 1, j + 1] - rr[i + 1, j - 1] + rr[i - 1, j - 1]) /
                (4.0 * dxiavg * detaavg)

            # calculate alpha, beta, and gamma (as defined in theory document)
            alpha = x_eta^2 + r_eta^2
            beta = x_eta * x_xi + r_eta * r_xi
            gamma = x_xi^2 + r_xi^2

            # populate residual vector
            resid[1, i - 1, j - 1] =
                alpha * x_xixi - 2.0 * beta * x_xieta + gamma * x_etaeta -
                beta * r_xi * x_eta * detaminus * detaplus / ravg
            resid[2, i - 1, j - 1] =
                alpha * r_xixi - 2.0 * beta * r_xieta + gamma * r_etaeta -
                beta * r_xi * r_eta * detaminus * detaplus / ravg
        end
    end

    return r
end

function solve_elliptic_grid_iad(x, p)

    # - Wrap residual - #
    function rwrap!(r, y)
        return elliptic_grid_residual!(r, y, x, p)
    end

    wake_grid = reshape(x, p.gridshape)

    # - Call NLsolve - #
    # df = OnceDifferentiable(rwrap!, x, similar(x))
    result = NLsolve.nlsolve(
        # df,
        rwrap!,
        reshape(@view(wake_grid[:, 2:end, 2:(end - 1)]), :);
        # method=:newton,
        # autodiff=:forward,
        method=p.algorithm,
        autodiff=p.autodiff,
        linsolve=(x, A, b) -> x .= ImplicitAD.implicit_linear(A, b),
        ftol=p.atol,
        iterations=p.iteration_limit,
        show_trace=p.verbose,
    )

    p.converged[1] = NLsolve.converged(result)

    # - overwrite output in place - #
    wake_grid[:, 2:end, 2:(end - 1)] .= reshape(result.zero, p.itshape)

    return x
end

"""
TODO: will want to pass a convergence flag that get's updated so we can break out if the wake geometry did not converge and pass a fail flag to the optimizer
"""
function solve_elliptic_grid!(
    wake_grid;
    algorithm=:newton,
    autodiff=:forward,
    atol=1e-14,
    iteration_limit=10,
    converged=[false],
    verbose=false,
)

    # - dimensions - #
    gridshape = size(wake_grid)
    nx = gridshape[2]
    nr = gridshape[3]

    # - precomputation - #
    eta = zeros(eltype(wake_grid), nr)
    for j in 2:nr
        @views eta[j] = wake_grid[2, 1, j] .^ 2 - wake_grid[2, 1, j - 1] .^ 2 + eta[j - 1]
    end
    xi = @view(wake_grid[1, :, 1]) .- @view(wake_grid[1, 1, 1])

    # - set up solve - #
    p = (;
        eta,
        xi,
        gridshape,
        proposed_grid_cache=DiffCache(wake_grid),
        itshape=(gridshape[1], gridshape[2] - 1, gridshape[3] - 2),
        algorithm,
        autodiff,
        atol,
        iteration_limit,
        converged,
        verbose,
    )

    # - solve - #
    y = implicit(solve_elliptic_grid_iad, elliptic_grid_residual!, reshape(wake_grid, :), p)

    for g in eachrow(view(y, :))
        if g[1] < eps()
            g[1] = 0.0
        end
    end

    # - format outputs - #
    wake_grid .= reshape(y, gridshape)

    return wake_grid
end

#---------------------------------#
#    DFDC-like wake relaxation    #
#---------------------------------#
"""
    relax_grid!(xg, rg, nxi, neta; relaxation_iteration_limit, relaxation_atol)

Relax wake_grid using elliptic wake_grid solver.

# Arguments:
 - `xg::Matrix{Float64}` : Initial x wake_grid points guess
 - `rg::Matrix{Float64}` : Initial r wake_grid points guess
 - `nxi::Int` : number of xi (x) stations in the wake_grid
 - `neta::Int` : number of eta (r) stations in the wake_grid

# Keyword Arguments:
 - `relaxation_iteration_limit::Int` : maximum number of iterations to run, default=100
 - `relaxation_atol::Float` : convergence tolerance, default = 1e-9

# Returns:
 - `x_relax_points::Matrix{Float64}` : Relaxed x wake_grid points
 - `r_relax_points::Matrix{Float64}` : Relaxed r wake_grid points
"""
function relax_grid!(
    wake_grid;
    relaxation_iteration_limit=100,
    relaxation_atol=1e-9,
    converged,
    verbose=false,
    silence_warnings=true,
    ntab=1,
    tabchar="    ",
)
    xr = view(wake_grid, 1, :, :)
    rr = view(wake_grid, 2, :, :)

    TF = eltype(wake_grid)
    _, nxi, neta = size(wake_grid)

    # initialize output and computation arrays
    # xr = [xg[i, j] for i in 1:nxi, j in 1:neta]
    # rr = [rg[i, j] for i in 1:nxi, j in 1:neta]
    C = zeros(TF, nxi)
    D = zeros(TF, 2, nxi)

    #set up relaxation factors
    #TODO: how are these decided?
    if relaxation_iteration_limit > 0
        relaxfactor1 = 1.0
        relaxfactor2 = 1.1
        relaxfactor3 = 1.4
    else
        relaxfactor1 = 1.0
        relaxfactor2 = 1.0
        relaxfactor3 = 1.0
        relaxation_iteration_limit = 1
    end

    relaxfactor = relaxfactor1

    # convergence tolerances for each phase
    # #TODO: where do these come from?
    dset1 = 1e-1
    dset2 = 5e-3
    dset3 = 5e-7

    #initialize streamline values (used to calculate derivatives later)
    eta = zeros(TF, neta)
    for j in 2:neta
        # DFDC comment for this is roughly: streamline values (axisymmetric streamfunction). eta(j) = (j-1)**2 gives uniform spacing in r
        # not sure why it doesn't match the actual expression though.
        eta[j] = view(rr, 1, j) .^ 2 - view(rr, 1, j - 1) .^ 2 + eta[j - 1]
    end

    # Get the x positions along the x axis (used later for various ratios compared to internal points as they move)
    xi = [xr[i, 1] for i in 1:nxi]
    # Get the r positions along the inner radial boundary
    # centerbodyrstations = [rr[i, 1] for i in 1:nxi]
    centerbodyrstations = view(rr, 1:nxi, 1)

    # used as scaling for tolerance
    dxy = max(
        xr[1, 1],
        rr[1, 1],
        abs(xi[end] - xi[1]),
        abs(centerbodyrstations[end] - centerbodyrstations[1]),
    )

    #next dfdc goes to AXELL function
    #skip over most of the stuff since the intlet and walls don't get relaxed

    for iterate in 1:relaxation_iteration_limit
        if verbose
            println(tabchar^(ntab) * "iteration $iterate")
        end

        #initialize max step
        dmax = 0.0

        #loop over interior streamlines:
        for j in 2:(neta - 1)
            #rename indices for convenience
            jm1 = j - 1
            jp1 = j + 1

            ## -- Relax Outlet -- ##

            #rename outlet points for convenience
            xej = xr[end, j]
            xem1j = xr[end - 1, j]
            xem2j = xr[end - 2, j]

            rej = rr[end, j]
            rem1j = rr[end - 1, j]
            rem2j = rr[end - 2, j]

            #x arc distance between j+1 and j-1 radial positions of outlet
            xarc = xr[end, jp1] - xr[end, jm1]
            #r arc distances between j+1 and j-1 radial positions of outlet
            rarc = rr[end, jp1] - rr[end, jm1]

            # Check if sqrt of zero will take place
            if isapprox(xarc, rarc)
                xhat = TF(0.0)
                rhat = TF(0.0)
            else
                # normalized x component of arc length (x direction) of current position on outlet
                xhat = xarc / sqrt(xarc^2 + rarc^2)

                # normalized r component of arc length (r direction) of current position on outlet
                rhat = rarc / sqrt(xarc^2 + rarc^2)
            end

            #xi distance between outlet and second to last xi station before outlet
            dxem1 = xi[end] - xi[end - 1]
            #xi distance along inner bound between second and third to last xi stations before outlet.
            dxem2 = xi[end - 1] - xi[end - 2]

            #ratios of x and r spacing comparing current interior r position vs inner bound x position
            dxem1dxinb = (xej - xem1j) / dxem1
            dxem2dxinb = (xem1j - xem2j) / dxem2
            drem1dxinb = (rej - rem1j) / dxem1
            drem2dxinb = (rem1j - rem2j) / dxem2

            #2nd-order 3-point difference to get tangential velocity
            #TODO: figure out what is happening here.
            rez =
                xhat * (dxem1dxinb + dxem1 * (dxem1dxinb - dxem2dxinb) / (dxem1 + dxem2)) +
                rhat * (drem1dxinb + dxem1 * (drem1dxinb - drem2dxinb) / (dxem1 + dxem2))

            x_xij = xhat * (1.0 / dxem1 + 1.0 / (dxem1 + dxem2))
            x_rij = rhat * (1.0 / dxem1 + 1.0 / (dxem1 + dxem2))
            x_xi = x_xij * xhat + x_rij * rhat

            # outletrelax = 1.0 * relaxfactor
            doutlet = -relaxfactor * rez / x_xi

            xr[end, j] += doutlet * xhat
            rr[end, j] += doutlet * rhat

            ## -- END Outlet Relaxation -- ##

            ## -- Relax points on the jth streamline -- ##
            #TODO:NEED TO RECONCILE WRITTEN THEORY WITH CODE...
            #TODO: can all the derivatives be replaced with FLOWMath or similar?

            #march down xi direction
            for i in 2:(nxi - 1)

                #rename indices for convenience
                im1 = i - 1
                ip1 = i + 1

                #part of the first derivative calculations
                dximinus = xi[i] - xi[i - 1]
                dxiplus = xi[i + 1] - xi[i]
                dxiavg = 0.5 * (dximinus + dxiplus)

                #part of the first derivative calculations
                detaminus = eta[j] - eta[j - 1]
                detaplus = eta[j + 1] - eta[j]
                detaavg = 0.5 * (detaminus + detaplus)

                # - Calculate 1st Derivatives (non-uniform spacing)
                x_eta = (xr[i, j + 1] - xr[i, j - 1]) / (2.0 * detaavg)
                r_eta = (rr[i, j + 1] - rr[i, j - 1]) / (2.0 * detaavg)
                x_xi = (xr[i + 1, j] - xr[i - 1, j]) / (2.0 * dxiavg)
                r_xi = (rr[i + 1, j] - rr[i - 1, j]) / (2.0 * dxiavg)

                #alpha, beta, and gamma as defined in theory doc.
                alpha = x_eta^2 + r_eta^2
                beta = x_eta * x_xi + r_eta * r_xi
                gamma = x_xi^2 + r_xi^2
                # J = r_eta * x_xi - x_eta * r_xi

                # - Calculate 2nd Derivatives

                #Used in calculating second derivatives with non-constant intervals.
                ravgplus = 0.5 * (rr[i, j] + rr[i, j + 1])
                ravgminus = 0.5 * (rr[i, j] + rr[i, j - 1])
                ravg = 0.5 * (ravgplus + ravgminus)

                ximinuscoeff = detaminus * detaplus / (dximinus * dxiavg)
                xipluscoeff = detaminus * detaplus / (dxiplus * dxiavg)
                etaminuscoeff = detaplus / detaavg * ravgminus / ravg
                etapluscoeff = detaminus / detaavg * ravgplus / ravg

                #x_ξξ
                x_xixi =
                    (xr[i + 1, j] - xr[i, j]) * xipluscoeff -
                    (xr[i, j] - xr[i - 1, j]) * ximinuscoeff

                #x_ηη
                x_etaeta =
                    (xr[i, j + 1] - xr[i, j]) * etapluscoeff -
                    (xr[i, j] - xr[i, j - 1]) * etaminuscoeff

                #x_ξη
                x_xieta =
                    detaminus *
                    detaplus *
                    (
                        xr[i + 1, j + 1] - xr[i - 1, j + 1] - xr[i + 1, j - 1] +
                        xr[i - 1, j - 1]
                    ) / (4.0 * dxiavg * detaavg)

                #r_ξξ
                r_xixi =
                    (rr[i + 1, j] - rr[i, j]) * xipluscoeff -
                    (rr[i, j] - rr[i - 1, j]) * ximinuscoeff

                #r_ηη
                r_etaeta =
                    (rr[i, j + 1] - rr[i, j]) * etapluscoeff -
                    (rr[i, j] - rr[i, j - 1]) * etaminuscoeff

                #r_ξη
                r_xieta =
                    detaminus *
                    detaplus *
                    (
                        rr[i + 1, j + 1] - rr[i - 1, j + 1] - rr[i + 1, j - 1] +
                        rr[i - 1, j - 1]
                    ) / (4.0 * dxiavg * detaavg)

                #D's look to be from eqns 92 and 93 in theory doc, but don't actually match up completely.  They don't even match up with the notes in the code...
                D[1, i] =
                    alpha * x_xixi - 2.0 * beta * x_xieta + gamma * x_etaeta -
                    beta * r_xi * x_eta * detaminus * detaplus / ravg

                D[2, i] =
                    alpha * r_xixi - 2.0 * beta * r_xieta + gamma * r_etaeta -
                    beta * r_xi * r_eta * detaminus * detaplus / ravg

                #SLOR Stuff
                #TODO: replace with linearsolve package?
                #TODO: what are A, B, and C?
                A =
                    alpha * (ximinuscoeff + xipluscoeff) +
                    gamma * (etaminuscoeff + etapluscoeff)

                if i == 2
                    B = 0 #TODO: Why?
                else
                    B = -alpha * ximinuscoeff
                end

                C[i] = -alpha * xipluscoeff

                ainv = 1.0 / (A - B * C[i - 1])
                C[i] *= ainv
                D[1, i] = (D[1, i] - B * D[1, i - 1]) * ainv
                D[2, i] = (D[2, i] - B * D[2, i - 1]) * ainv
            end #for i (streamwise stations)

            # Rest of SLOR Stuff
            # run through streamline backwards?
            for ib in 2:(nxi - 1)
                i = nxi - ib + 1

                D[1, i] -= C[i] * D[1, i + 1]
                D[2, i] -= C[i] * D[2, i + 1]

                #over relaxation step?
                xr[i, j] += relaxfactor * D[1, i]
                rr[i, j] += relaxfactor * D[2, i]

                # stuff for convergence checking
                ad1 = abs(D[1, i])
                ad2 = abs(D[2, i])

                dmax = max(dmax, ad1, ad2)
            end #for ib
        end #for j (radial stations)

        # -- Update relaxation factors
        #TODO: need to figure out how the dset numbers and relaxation factors are chosen.  Why are they the values that they are? Does it have something to do with the SLOR setup?
        if dmax < relaxation_atol * dxy
            if verbose
                println(tabchar^(ntab) * "Total iterations: $iterate")
            end
            converged[1] = true
            return wake_grid
        end

        relaxfactor = relaxfactor1

        if dmax < dset1 * dxy
            relaxfactor = relaxfactor2
        end

        if dmax < dset2 * dxy
            relaxfactor = relaxfactor3
        end

        if dmax < dset3 * dxy
            if verbose
                println(tabchar^(ntab) * "Total iterations = $iterate")
            end
            converged[1] = true
            return wake_grid
        end
    end

    if verbose
        println(tabchar^(ntab) * "Total iterations = ", relaxation_iteration_limit, "\n")
    end

    if !silence_warnings
        @warn "Wake grid relaxation did not converge, iteration limit of $(relaxation_iteration_limit) met."
    end

    return wake_grid
end
