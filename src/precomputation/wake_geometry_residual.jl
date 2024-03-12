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
        method=:newton,
        autodiff=:forward,
        linsolve=(x, A, b) -> x .= ImplicitAD.implicit_linear(A, b), #used in newton method, unused otherwise so fine to keep it here.
        ftol=p.wake_nlsolve_ftol,
        iterations=p.wake_max_iter,
        show_trace=p.verbose,
    )

    # - overwrite output in place - #
    wake_grid[:, 2:end, 2:(end - 1)] .= reshape(result.zero, p.itshape)

    return x
end

"""
TODO: will want to pass a convergence flag that get's updated so we can break out if the wake geometry did not converge and pass a fail flag to the optimizer
"""
function solve_elliptic_grid!(
    wake_grid; wake_nlsolve_ftol=1e-14, wake_max_iter=10, verbose=false
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
        wake_nlsolve_ftol,
        wake_max_iter,
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
