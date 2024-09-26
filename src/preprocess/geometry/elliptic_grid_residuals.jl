"""
    backward_stencil_1(f, x)

Three-point, backward difference scheme on a non-uniform grid.

# Arguments:
- `f::AbstractMatrix{Float}` : f(x) for x_{i-2}:x_{i}
- `x::AbstractVector{Float}` : [x_{i-2}:x_{i}]

# Returns:
- `f'(x)::Float` : ∂f/∂x at x_i
"""
function backward_stencil_1(f, x)
    #1 = n-2
    #2 = n-1
    #3 = n

    # First Term
    n1 = (x[3] - x[2])
    d1 = (x[3] - x[1]) * (x[2] - x[1])
    t1 = f[1] * n1 / d1

    # Second Term
    n2 = (x[1] - x[3])
    d2 = (x[2] - x[1]) * (x[3] - x[2])
    t2 = f[2] * n2 / d2

    # Third Term
    n3 = (2.0 * x[3] - x[2] - x[1])
    d3 = (x[3] - x[2]) * (x[3] - x[1])
    t3 = f[3] * n3 / d3

    return t1 + t2 + t3
end

"""
    center_stencil_1(f, x)

Three-point, central difference scheme on a non-uniform grid, for first derivatives.

# Arguments:
- `f::AbstractMatrix{Float}` : f(x) for x_{i-1}:x_{i+1}
- `x::AbstractVector{Float}` : [x_{i-1}:x_{i+1}]

# Returns:
- `f'(x)::Float` : ∂f/∂x at x_i
"""
function center_stencil_1(f, x)

    #1 = i-1
    #2 = i
    #3 = i+1

    return (f[3] - f[1]) / (x[3] - x[1])
end

"""
    center_stencil_2(f, x)

Three-point, central difference scheme on a non-uniform grid, for second derivatives.

# Arguments:
- `f::AbstractMatrix{Float}` : f(x) for x_{i-1}:x_{i+1}
- `x::AbstractVector{Float}` : [x_{i-1}:x_{i+1}]

# Returns:
- `f''(x)::Float` : ∂^2f/∂x^2 at x_i
"""
function center_stencil_2(f, x)

    #1 = i-1
    #2 = i
    #3 = i+1

    h2 = (x[2] - x[1])
    h3 = (x[3] - x[2])

    n1 = h2 * (f[3] - f[2])
    n2 = h3 * (f[2] - f[1])
    d1 = h2 * h3 * (h3 + h2)

    return 2.0 * (n1 - n2) / d1
end

"""
    center_stencil_2_mixed(f, x, y)

Three-point, central difference scheme on a non-uniform grid, for mixed derivatives.
Uses a consecutive first derivative central difference (`center_stencil_1`) in the x-dimension then another centeral difference of the f'(x) values in the y-dimension.

# Arguments:
- `f::AbstractMatrix{Float}` : f(x,y) for x_{i-1}:x_{i+1} and y_{i-1}:y_{i+1}
- `x::AbstractVector{Float}` : [x_{i-1}:x_{i+1}]
- `y::AbstractVector{Float}` : [y_{i-1}:y_{i+1}]

# Returns:
- `f'(x,y)::Float` : ∂^2f/∂x∂y at x_ij
"""
function center_stencil_2_mixed(f, x, y)
    TF = promote_type(eltype(f), eltype(x), eltype(y))
    dfdx = zeros(TF, 3)

    for j in 1:3
        dfdx[j] = center_stencil_1(f[:, j], x)
    end

    return center_stencil_1(dfdx, y)
end

"""
assume grid has radial lines of constant axial position such that z'(ξ)=1.0 and z''(ξ)=z`(η)=z''(η)=0.0
"""
function elliptic_grid_residual!(r, y, x, p)

    # # - extract parameters - #
    # allocating version
    # (; x_caching, itshape) = p
    # (; x_dims) = x_caching
    # TF = promote_type(eltype(x), eltype(y))
    # const_grid, ξ, η = withdraw_grid_parameter_cache(x, x_dims)
    # proposed_grid = zeros(TF, size(const_grid))
    # proposed_grid .= const_grid

    # - extract parameters - #
    # use preallocated proposed grid to save a few allocations
    (; x_caching, itshape) = p
    (; x_cache, x_dims) = x_caching
    x_vec = PreallocationTools.get_tmp(x_cache, promote_type(eltype(y), eltype(x))(1.0))
    x_vec .= x
    proposed_grid, ξ, η = withdraw_grid_parameter_cache(x_vec, x_dims)

    # dimensions
    nxi = x_dims.xi.shape[1]
    neta = x_dims.eta.shape[1]

    # - Initialize Residual - #
    r .= 0

    # - reshape vector into 3D array - #
    resid = reshape(r, itshape)

    #overwrite proposed_grid internals with the current state variables
    proposed_grid[2, 2:end, 2:(end - 1)] .= reshape(y, itshape)

    # separate x and r components of proposed grid
    zs = view(proposed_grid, 1, :, :)
    rs = view(proposed_grid, 2, :, :)

    # - Populate Residual Vector - #

    # loop over interior streamlines (center body and duct wall points are set by Dirichlet Boundary conditions)
    for j in 2:(neta - 1)

        # --- Outlet Residuals --- #
        #=
        Note: Neumann Boundary Condition:
        Velocity is tangent to streamline
        Streamlines are along lines of constant η
        ∴ Velocity components in η̂ direction (perpendicular to streamlines) must go to zero
        =#

        # - η unit vector components in z and r directions - #

        # change in z between points radially above and below current point
        dz = zs[end, j + 1] - zs[end, j - 1]
        # change in r between points radially above and below current point
        dr = rs[end, j + 1] - rs[end, j - 1]
        # distance between grid points radially above and below current point
        enorm = sqrt(dz^2 + dr^2)

        # unit vector components
        # note: set to zero if you'd divide by zero
        η̂z = iszero(enorm) ? 0.0 : dz / enorm
        η̂r = iszero(enorm) ? 0.0 : dr / enorm

        # Calculate derivatives using spline interpolation
        r_ξ_outlet = backward_stencil_1(rs[(end - 2):end, j], ξ[(end - 2):end])

        # Boundary condition is that velocity is aligned with streamlines
        #=
        In other words, we want the velocity perpendicular to the streamline to be zero
        Cz = dz_dξ_outlet/(ρrJ)
        Cr = dr_dξ_outlet/(ρrJ)
        V⟂ = η̂z*Cz + η̂r*Cr (= 0)
        V⟂ = η̂z*dz_dξ_outlet/(ρrJ) + η̂r*dr_dξ_outlet/(ρrJ) (= 0)
        V⟂*ρ*r*J = η̂z*dz_dξ_outlet + η̂r*dr_dξ_outlet (= 0)
        V⟂ = η̂z*dz_dξ_outlet + η̂r*dr_dξ_outlet (= 0)
        =#
        Vperp = η̂z + η̂r * r_ξ_outlet

        # for residual, need to return a component relative to the associated state
        # (z for index 1, r for index 2)
         # resid[1, end, j - 1] = η̂z * Vperp
         # resid[2, end, j - 1] = η̂r * Vperp
         resid[end, j - 1] = η̂r * Vperp

        # - Interior Residuals - #
        # Loop over interior points between rotor plane and outlet plane
        # (rotor plane is set by Dirichlet boundary condition
        # and outlet plane was just handled above with Neumann boundary condition)
        for i in 2:(nxi - 1)

            ### --- Rename For Readability --- ###

            # - First Derivative Terms - #
            r_ξ = center_stencil_1(rs[(i - 1):(i + 1), j], ξ[(i - 1):(i + 1)])
            r_η = center_stencil_1(rs[i, (j - 1):(j + 1)], η[(j - 1):(j + 1)])

            # assemble alpha, beta, and gamma (convenience variables)
            α = r_η^2
            β = r_η * r_ξ
            γ = 1.0 + r_ξ^2

            # - Second and Multivariate Derivative Terms - #
            # second
            r_ξξ = center_stencil_2(rs[(i - 1):(i + 1), j], ξ[(i - 1):(i + 1)])
            r_ηη = center_stencil_2(rs[i, (j - 1):(j + 1)], η[(j - 1):(j + 1)])

            # multivariate
            r_ξη = center_stencil_2_mixed(
                rs[(i - 1):(i + 1), (j - 1):(j + 1)], ξ[(i - 1):(i + 1)], η[(j - 1):(j + 1)]
            )

            resid[i - 1, j - 1] =
                α * r_ξξ - 2.0 * β * r_ξη + γ * r_ηη + γ * r_η^2 / rs[i, j] -
                β * r_ξ * r_η / rs[i, j]
        end
    end

    return r
end

"""
"""
function elliptic_grid_residual_dfdc!(r, y, x, p)

    # - extract parameters - #
    (; x_caching, itshape) = p
    (; x_dims) = x_caching

    const_grid, xi, eta = withdraw_grid_parameter_cache(x, x_dims)
    proposed_grid = zeros(promote_type(eltype(x), eltype(y)), size(const_grid))
    proposed_grid .= const_grid

    # dimensions
    nxi = x_dims.xi.shape[1]
    neta = x_dims.eta.shape[1]

    # - Initialize Residual - #
    r .= 0

    # - reshape vector into 3D array - #
    resid = reshape(r, itshape)

    #overwrite proposed_grid internals with the current state variables
    proposed_grid[:, 2:end, 2:(end - 1)] .= reshape(y, itshape)

    # separate x and r components of proposed grid
    zs = view(proposed_grid, 1, :, :)
    rs = view(proposed_grid, 2, :, :)

    # loop over interior streamlines
    for j in 2:(neta - 1)

        # --- Outlet Residuals --- #

        # Neumann Boundary Condition (u^2 + v^2 = constant)

        # - η unit vector components in z and r directions - #
        # change in z between points radially above and below current point
        # note that for vertical lines, this will be zero
        dz = zs[end, j + 1] - zs[end, j - 1]
        # change in r between points radially above and below current point
        dr = rs[end, j + 1] - rs[end, j - 1]

        # distance between grid points radially above and below current point
        enorm = sqrt(dz^2 + dr^2)

        # unit vector components
        # note: set to zero if you'd divide by zero
        zhat = iszero(enorm) ? 0.0 : dz / enorm
        rhat = iszero(enorm) ? 0.0 : dr / enorm

        # - calculate 1st derivatives using finite differencing on irregular grid - #
        # Note: we're at the end of the grid, so we have to take derivatives based on points ahead axially

        # change in xi coordinates from second to last to last grid point on the streamline (axial direction)
        dxi1 = xi[end] - xi[end - 1]
        # change in xi coordinates from third to last to second to last grid point on the streamline (axial direction)
        dxi2 = xi[end - 1] - xi[end - 2]

        # change in z coordinates from second to last to last grid point on the streamline (axial direction)
        dz1 = zs[end, j] - zs[end - 1, j]
        # change in z coordinates from third to last to second to last grid point on the streamline (axial direction)
        dz2 = zs[end - 1, j] - zs[end - 2, j]

        # change in r coordinates from second to last to last grid point on the streamline (axial direction)
        dr1 = rs[end, j] - rs[end - 1, j]
        # change in r coordinates from third to last to second to last grid point on the streamline (axial direction)
        dr2 = rs[end - 1, j] - rs[end - 2, j]

        # rename for convenience
        dzdxi_1 = dz1 / dxi1
        dzdxi_2 = dz2 / dxi2
        drdxi_1 = dr1 / dxi1
        drdxi_2 = dr2 / dxi2

        # - 3-point, second order backward finite diff to get tangential velocity components - #
        # cz = dz_dξ/(ρrJ)
        dxdxi = dzdxi_1 + dxi1 * (dzdxi_1 - dzdxi_2) / (dxi1 + dxi2)
        drdxi = drdxi_1 + dxi1 * (drdxi_1 - drdxi_2) / (dxi1 + dxi2)


        # - populate residual vector - #
        tmp = zhat * dxdxi + rhat * drdxi

        resid[1, end, j - 1] = zhat * tmp
        resid[2, end, j - 1] = rhat * tmp

        # -- Interior Residuals --- #

        for i in 2:(nxi - 1)

            # calculate 1st derivatives using finite differencing on an irregular grid
            dximinus = xi[i] - xi[i - 1]
            dxiplus = xi[i + 1] - xi[i]
            dxiavg = 0.5 * (dximinus + dxiplus)

            detaminus = eta[j] - eta[j - 1]
            detaplus = eta[j + 1] - eta[j]
            detaavg = 0.5 * (detaminus + detaplus)

            z_eta = (zs[i, j + 1] - zs[i, j - 1]) / (2.0 * detaavg)
            r_eta = (rs[i, j + 1] - rs[i, j - 1]) / (2.0 * detaavg)
            z_xi = (zs[i + 1, j] - zs[i - 1, j]) / (2.0 * dxiavg)
            r_xi = (rs[i + 1, j] - rs[i - 1, j]) / (2.0 * dxiavg)

            # calculate 2nd derivatives using finite differencing on an irregular grid
            ravgplus = 0.5 * (rs[i, j] + rs[i, j + 1])
            ravgminus = 0.5 * (rs[i, j] + rs[i, j - 1])
            ravg = 0.5 * (ravgplus + ravgminus)

            ximinuscoeff = detaminus * detaplus / (dximinus * dxiavg)
            xipluscoeff = detaminus * detaplus / (dxiplus * dxiavg)
            etaminuscoeff = detaplus / detaavg * ravgminus / ravg
            etapluscoeff = detaminus / detaavg * ravgplus / ravg

            z_xixi =
                (zs[i + 1, j] - zs[i, j]) * xipluscoeff -
                (zs[i, j] - zs[i - 1, j]) * ximinuscoeff
            z_etaeta =
                (zs[i, j + 1] - zs[i, j]) * etapluscoeff -
                (zs[i, j] - zs[i, j - 1]) * etaminuscoeff
            z_xieta =
                detaminus *
                detaplus *
                (zs[i + 1, j + 1] - zs[i - 1, j + 1] - zs[i + 1, j - 1] + zs[i - 1, j - 1]) /
                (4.0 * dxiavg * detaavg)
            r_xixi =
                (rs[i + 1, j] - rs[i, j]) * xipluscoeff -
                (rs[i, j] - rs[i - 1, j]) * ximinuscoeff
            r_etaeta =
                (rs[i, j + 1] - rs[i, j]) * etapluscoeff -
                (rs[i, j] - rs[i, j - 1]) * etaminuscoeff
            r_xieta =
                detaminus *
                detaplus *
                (rs[i + 1, j + 1] - rs[i - 1, j + 1] - rs[i + 1, j - 1] + rs[i - 1, j - 1]) /
                (4.0 * dxiavg * detaavg)

            # calculate alpha, beta, and gamma
            alpha = z_eta^2 + r_eta^2
            beta = z_eta * z_xi + r_eta * r_xi
            gamma = z_xi^2 + r_xi^2

            # populate residual vector

            resid[1, i - 1, j - 1] =
                alpha * z_xixi - 2.0 * beta * z_xieta + gamma * z_etaeta -
                beta * r_xi * z_eta * detaminus * detaplus / ravg
            resid[2, i - 1, j - 1] =
                alpha * r_xixi - 2.0 * beta * r_xieta + gamma * r_etaeta -
                beta * r_xi * r_eta * detaminus * detaplus / ravg

        end
    end

    return r
end

"""
"""
function solve_elliptic_grid(x, p)

    # - Wrap residual - #
    function rwrap!(r, y)
        return elliptic_grid_residual!(r, y, x, p)
    end

    wake_grid, _, _ = withdraw_grid_parameter_cache(x, p.x_caching.x_dims)

    # - Call NLsolve - #
    result = NLsolve.nlsolve(
        rwrap!,
        reshape(@view(wake_grid[2, 2:end, 2:(end - 1)]), :);
        method=p.algorithm,
        autodiff=p.autodiff,
        linsolve=(x, A, b) -> x .= ImplicitAD.implicit_linear(A, b),
        ftol=p.atol,
        iterations=p.iteration_limit,
        show_trace=p.verbose,
    )

    p.converged[1] = NLsolve.converged(result)

    return result.zero
end

"""
"""
function solve_elliptic_grid!(
    wake_grid;
    algorithm=:trust_region,
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

    x_caching = allocate_grid_parameter_cache(wake_grid, xi, eta)

    # - set up solve - #
    itshape = (gridshape[2] - 1, gridshape[3] - 2)
    constants = (;
        x_caching, itshape, algorithm, autodiff, atol, iteration_limit, converged, verbose
    )

    # non ImplicitAD version
    # grid_internals = solve_elliptic_grid([reshape(wake_grid, :); xi; eta], constants)

    grid_internals = ImplicitAD.implicit(
        solve_elliptic_grid,
        elliptic_grid_residual!,
        [reshape(wake_grid, :); xi; eta],
        constants,
    )

    # Place solved grid internals into wake grid
    wake_grid[2, 2:end, 2:(end - 1)] .= reshape(grid_internals, itshape)

    # Check that no weird things happened with values being placed beneath the axis of rotation
    for g in eachindex(wake_grid[2, :, :])
        if wake_grid[g] < eps()
            wake_grid[g] = 0.0
        end
    end

    return wake_grid
end

#---------------------------------#
#    DFDC-like wake relaxation    #
#---------------------------------#
"""
    relax_grid!(xg, rg, nxi, neta; iteration_limit, atol)

Relax wake_grid using elliptic wake_grid solver.

# Arguments:

# Keyword Arguments:
 - `iteration_limit::Int` : maximum number of iterations to run, default=100
 - `atol::Float` : convergence tolerance, default = 1e-9

# Returns:
"""
function relax_grid!(
    wake_grid;
    iteration_limit=1000,
    atol=eps(),
    converged=[false],
    verbose=false,
    silence_warnings=true,
    ntab=1,
    tabchar="    ",
)
    zs = view(wake_grid, 1, :, :)
    rs = view(wake_grid, 2, :, :)

    TF = eltype(wake_grid)
    _, nxi, neta = size(wake_grid)

    # initialize output and computation arrays
    # zs = [xg[i, j] for i in 1:nxi, j in 1:neta]
    # rs = [rg[i, j] for i in 1:nxi, j in 1:neta]
    C = zeros(TF, nxi)
    D = zeros(TF, 2, nxi)

    #set up relaxation factors
    if iteration_limit > 0
        relaxfactor1 = 1.0
        relaxfactor2 = 1.1
        relaxfactor3 = 1.4
    else
        relaxfactor1 = 1.0
        relaxfactor2 = 1.0
        relaxfactor3 = 1.0
        iteration_limit = 1
    end

    relaxfactor = relaxfactor1

    # convergence tolerances for each phase
    dset1 = 1e-1
    dset2 = 5e-3
    dset3 = 5e-7

    #initialize streamline values (used to calculate derivatives later)
    eta = zeros(TF, neta)
    for j in 2:neta
        # DFDC comment for this is roughly: streamline values (axisymmetric streamfunction). eta(j) = (j-1)**2 gives uniform spacing in r
        # not sure why it doesn't match the actual expression though.
        # NOTE: does this require uniform blade element spacing?
        eta[j] = view(rs, 1, j) .^ 2 - view(rs, 1, j - 1) .^ 2 + eta[j - 1]
    end

    # Get the x positions along the x axis (used later for various ratios compared to internal points as they move)
    xi = [zs[i, 1] for i in 1:nxi]
    # Get the r positions along the inner radial boundary
    # centerbodyrstations = [rs[i, 1] for i in 1:nxi]
    centerbodyrstations = view(rs, 1:nxi, 1)

    # used as scaling for tolerance
    dxy = max(
        zs[1, 1],
        rs[1, 1],
        abs(xi[end] - xi[1]),
        abs(centerbodyrstations[end] - centerbodyrstations[1]),
    )

    #next dfdc goes to AXELL function
    #skip over most of the stuff since the intlet and walls don't get relaxed

    for iterate in 1:iteration_limit
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
            xej = zs[end, j]
            xem1j = zs[end - 1, j]
            xem2j = zs[end - 2, j]

            rej = rs[end, j]
            rem1j = rs[end - 1, j]
            rem2j = rs[end - 2, j]

            #x arc distance between j+1 and j-1 radial positions of outlet
            xarc = zs[end, jp1] - zs[end, jm1]
            #r arc distances between j+1 and j-1 radial positions of outlet
            rarc = rs[end, jp1] - rs[end, jm1]

            # Check if sqrt of zero will take place
            if isapprox(xarc, rarc)
                zhat = TF(0.0)
                rhat = TF(0.0)
            else
                # normalized x component of arc length (x direction) of current position on outlet
                zhat = xarc / sqrt(xarc^2 + rarc^2)

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
            rez =
                zhat * (dxem1dxinb + dxem1 * (dxem1dxinb - dxem2dxinb) / (dxem1 + dxem2)) +
                rhat * (drem1dxinb + dxem1 * (drem1dxinb - drem2dxinb) / (dxem1 + dxem2))

            z_xj = zhat * (1.0 / dxem1 + 1.0 / (dxem1 + dxem2))
            z_rj = rhat * (1.0 / dxem1 + 1.0 / (dxem1 + dxem2))
            z_s = z_xj * zhat + z_rj * rhat

            # outletrelax = 1.0 * relaxfactor
            doutlet = -relaxfactor * rez / z_s

            zs[end, j] += doutlet * zhat
            rs[end, j] += doutlet * rhat

            ## -- END Outlet Relaxation -- ##

            ## -- Relax points on the jth streamline -- ##

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
                z_eta = (zs[i, j + 1] - zs[i, j - 1]) / (2.0 * detaavg)
                r_eta = (rs[i, j + 1] - rs[i, j - 1]) / (2.0 * detaavg)
                z_xi = (zs[i + 1, j] - zs[i - 1, j]) / (2.0 * dxiavg)
                r_xi = (rs[i + 1, j] - rs[i - 1, j]) / (2.0 * dxiavg)

                #alpha, beta, and gamma as defined in theory doc.
                alpha = z_eta^2 + r_eta^2
                beta = z_eta * z_xi + r_eta * r_xi
                gamma = z_xi^2 + r_xi^2
                # J = r_eta * z_xi - z_eta * r_xi

                # - Calculate 2nd Derivatives

                #Used in calculating second derivatives with non-constant intervals.
                ravgplus = 0.5 * (rs[i, j] + rs[i, j + 1])
                ravgminus = 0.5 * (rs[i, j] + rs[i, j - 1])
                ravg = 0.5 * (ravgplus + ravgminus)

                ximinuscoeff = detaminus * detaplus / (dximinus * dxiavg)
                xipluscoeff = detaminus * detaplus / (dxiplus * dxiavg)
                etaminuscoeff = detaplus / detaavg * ravgminus / ravg
                etapluscoeff = detaminus / detaavg * ravgplus / ravg

                #z_ξξ
                z_xixi =
                    (zs[i + 1, j] - zs[i, j]) * xipluscoeff -
                    (zs[i, j] - zs[i - 1, j]) * ximinuscoeff

                #z_ηη
                z_etaeta =
                    (zs[i, j + 1] - zs[i, j]) * etapluscoeff -
                    (zs[i, j] - zs[i, j - 1]) * etaminuscoeff

                #z_ξη
                z_xieta =
                    detaminus *
                    detaplus *
                    (
                        zs[i + 1, j + 1] - zs[i - 1, j + 1] - zs[i + 1, j - 1] +
                        zs[i - 1, j - 1]
                    ) / (4.0 * dxiavg * detaavg)

                #r_ξξ
                r_xixi =
                    (rs[i + 1, j] - rs[i, j]) * xipluscoeff -
                    (rs[i, j] - rs[i - 1, j]) * ximinuscoeff

                #r_ηη
                r_etaeta =
                    (rs[i, j + 1] - rs[i, j]) * etapluscoeff -
                    (rs[i, j] - rs[i, j - 1]) * etaminuscoeff

                #r_ξη
                r_xieta =
                    detaminus *
                    detaplus *
                    (
                        rs[i + 1, j + 1] - rs[i - 1, j + 1] - rs[i + 1, j - 1] +
                        rs[i - 1, j - 1]
                    ) / (4.0 * dxiavg * detaavg)

                #D's look to be from eqns 92 and 93 in theory doc, but don't actually match up completely.  They don't even match up with the notes in the code...
                D[1, i] =
                    alpha * z_xixi - 2.0 * beta * z_xieta + gamma * z_etaeta -
                    beta * r_xi * z_eta * detaminus * detaplus / ravg

                D[2, i] =
                    alpha * r_xixi - 2.0 * beta * r_xieta + gamma * r_etaeta -
                    beta * r_xi * r_eta * detaminus * detaplus / ravg

                #SLOR Stuff
                A =
                    alpha * (ximinuscoeff + xipluscoeff) +
                    gamma * (etaminuscoeff + etapluscoeff)

                if i == 2
                    B = 0
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
                zs[i, j] += relaxfactor * D[1, i]
                rs[i, j] += relaxfactor * D[2, i]

                # stuff for convergence checking
                ad1 = abs(D[1, i])
                ad2 = abs(D[2, i])

                dmax = max(dmax, ad1, ad2)
            end #for ib
        end #for j (radial stations)

        # -- Update relaxation factors
        if dmax < atol #* dxy
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

        if dmax < atol #* dxy
            if verbose
                println(tabchar^(ntab) * "Total iterations = $iterate")
            end
            converged[1] = true
            return wake_grid
        end
    end

    if verbose
        println(tabchar^(ntab) * "Total iterations = ", iteration_limit, "\n")
    end

    if !silence_warnings
        @warn "Wake grid relaxation did not converge, iteration limit of $(iteration_limit) met."
    end

    return wake_grid
end
