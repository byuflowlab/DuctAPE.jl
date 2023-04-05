function wake_grid_residuals!(resid, proposed_grid, original_grid, xi, eta)

    # wake grid dimensions
    _, nxi, neta = size(grid)

    # separate x and r components of proposed grid
    xr = view(proposed_grid, 1, :, :)
    rr = view(proposed_grid, 2, :, :)

    # initialize residuals
    resid .= proposed_grid .- original_grid

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

        # calculate 1st derivatives using finite differencing on irregular grid
        dxi1 = xi[end] - xi[end - 1]
        dxi2 = xi[end - 1] - xi[end - 2]

        dxdx1 = (xr[end, j] - xr[end - 1, j]) / dxi1
        dxdx2 = (xr[end - 1, j] - xr[end - 2, j]) / dxi2

        drdx1 = (rr[end, j] - rr[end - 1, j]) / dxi1
        drdx2 = (rr[end - 1, j] - rr[end - 2, j]) / dxi2

        x_xi = dxdx1 + dxi1 * (dxdx1 - dxdx2) / (dxi1 + dxi2)
        r_xi = drdx1 + dxi1 * (drdx1 - drdx2) / (dxi1 + dxi2)

        # # populate residual vector
        tmp = xhat * x_xi + rhat * r_xi
        resid[1, end, j] = xhat * tmp
        resid[2, end, j] = rhat * tmp

        # -- Interior Residuals --- #

        for i in 2:(nxi - 1)

            # calculate 1st derivatives using finite differencing on an irregular grid
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

            # calculate 2nd derivatives using finite differencing on an irregular grid
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
            resid[1, i, j] =
                alpha * x_xixi - 2.0 * beta * x_xieta + gamma * x_etaeta -
                beta * r_xi * x_eta * detaminus * detaplus / ravg
            resid[2, i, j] =
                alpha * r_xixi - 2.0 * beta * r_xieta + gamma * r_etaeta -
                beta * r_xi * r_eta * detaminus * detaplus / ravg
        end
    end

    return resid
end

nx = 10
nr = 10

x = range(0, 1; length=nx)
grid = zeros(2, nx, nr)
for j in 1:nr
    grid[1, :, j] .= x
    for i in 1:nx
        rt = 0.3 + 0.1 * (x[i])^2
        rb = 0.2 - 0.1 * (x[i])^2
        grid[2, i, j] = (j - 1) / (nr - 1) * rt + (nr - j) / (nr - 1) * rb
    end
end

vtk_save(vtk_grid("original", grid))

fresid = (x) -> reshape(wake_grid_residuals(reshape(x, size(grid)), grid), :)

x0 = reshape(grid, :)

result = nlsolve(fresid, x0)

x = result.zero

new_grid = reshape(x, size(grid))

vtk_save(vtk_grid("revised", new_grid))
