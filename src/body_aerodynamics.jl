"""
    calculate_body_strengths_residual(gamb, A_bb, b_bf, vx_bw, vr_bw, gamw)

Calculate body vortex strengths
kid is kutta indices, where kid[1] is the row/column to keep and kid[2] is the row/column to subtract from kid[1] and then delete (1 is the first duct panel, and 2 is the Nth duct panel)
"""
function calculate_body_vortex_strengths!(
    gamb, A_bb, b_bf, kidx, A_bw, gamw, A_br, sigr, ductwakeidx; debug=false
)

    # problem dimensions
    nwake = length(A_bw) # number of wake sheets
    nrotor = length(A_br) # number of rotors

    # add freestream contributions to right hand side
    # note: the negative was already included in the precomputation for the freestream.
    RHS = similar(gamb) .= b_bf

    if debug
        #initialize extra outputs
        bfree = copy(b_bf)
        bwake = similar(bfree) .= 0.0
        brotor = similar(bfree) .= 0.0
    end

    # add wake vortex sheet contributions to right hand side
    # note: the subtraction was not included in the other coefficients in the precomputation, so we need to subract
    for jwake in 1:nwake

        # get induced velocity in the x-direction
        RHS .-= A_bw[jwake] * view(gamw, jwake, :)
        if debug
            bwake .-= A_bw[jwake] * view(gamw, jwake, :)
        end
    end

    # add rotor source sheet contributions to right hand side
    # note: the subtraction is not included in the coefficients here, so we need to subract
    for jrotor in 1:nrotor

        # get induced velocity in the x-direction
        RHS .-= A_br[jrotor] * view(sigr, :, jrotor)
        if debug
            brotor .-= A_br[jrotor] * view(sigr, :, jrotor)
        end
    end

    #return the residual, rather than solving the linear system.
    #TODO: in the solve functions, need to make sure that we don't subtract things again.
    # return A * gamb - RHS

    # TODO: do these later after things are working if you need to get things faster.
    # TODO: decide if linear system needs to be solved or not
    # TODO: precompute factorization, test if implicit_linear is faster
    # return ldiv!(gamb, factorize(A_bb), RHS)

    # bk = [RHS[1:kidx[end]]; -gamw[ductwakeidx[1], ductwakeidx[2]]; RHS[(kidx[end] + 1):end]]
    bk = [RHS[1:kidx[end]]; 0.0; RHS[(kidx[end] + 1):end]]

    view(gamb, :) .= (A_bb \ bk)[1:end .∉ kidx[end] + 1]
    # view(gamb, :) .= solve_body_system(A_bb, RHS, kidx)

    if debug
        return bfree, bwake, brotor
    else
        return nothing
    end
end

"""
Apply Kutta Condition and Solve Linear System
"""
function solve_body_system(A_bb, RHS, kidx)

    # - Apply the Subtractive Kutta Condition - #
    for i in 1:length(kidx[:, 1])
        # subtract "Nth" row from "1st" row
        A_bb[kidx[i, 1], :] .-= A_bb[kidx[i, 2], :]
        # subtract "Nth" column from "1st" column
        A_bb[:, kidx[i, 1]] .-= A_bb[:, kidx[i, 2]]
        # subtract "Nth" row from "1st" row
        RHS[kidx[i, 1]] -= RHS[kidx[i, 2]]
    end

    # - Solve with subtractive Kutta Condition applied - #
    # x = ImplicitAD.implicit_linear(
    #     A_bb[1:end .∉ [kidx[:, 2]], 1:end .∉ [kidx[:, 2]]], RHS[1:end .∉ [kidx[:, 2]]]
    # )

    x = A_bb[1:end .∉ [kidx[:, 2]], 1:end .∉ [kidx[:, 2]]] \ RHS[1:end .∉ [kidx[:, 2]]]

    # - put things back after solving - #
    for i in 1:length(kidx[:, 1])
        # add "Nth" row back to "1st" row
        A_bb[kidx[i, 1], :] .+= A_bb[kidx[i, 2], :]
        # add "Nth" column back to "1st" column
        A_bb[:, kidx[i, 1]] .+= A_bb[:, kidx[i, 2]]
        # add "Nth" row back to "1st" row
        RHS[kidx[i, 1]] += RHS[kidx[i, 2]]
    end

    # - recover the final body panel strength value - #
    for i in 1:length(kidx[:, 1])
        insert!(x, kidx[i, 2], -x[kidx[i, 1]])
    end

    return x
end

######################################################################
#                                                                    #
#                          New Kutta setup                           #
#                                                                    #
######################################################################

"""
adds wake panel influence (from trailing edge panels) to the LHS matrix for the Kutta condition
"""
function body_lhs_kutta!(LHS, panels; tol=1e1 * eps(), verbose=false)
    (; endpoints, controlpoint, normal) = panels

    for (i, (te1, te2)) in
        enumerate(zip(eachrow(endpoints[:, 1, :]), eachrow(endpoints[:, 2, :])))

        # check that trailing edge points are coincident
        if norm(te2 - te1) < tol
            idxl = panels.endpointidxs[i, 1]
            idxu = panels.endpointidxs[i, 2]

            # Loop through control points being influenced
            for (m, (cp, nhat)) in enumerate(zip(eachrow(controlpoint), eachrow(normal)))

                # influence due to lower TE
                xi, rho, k2, rj = calculate_xrm(te1, cp)
                vx = -vortex_ring_vx(xi, rho, k2, rj, 19.5733 * rho)#lengths shouldn't be needed here, set such that self-induced case returns zero.
                vr = -vortex_ring_vr(xi, rho, k2, rj)
                LHS[m, idxl] -= dot([vx; vr], nhat)

                # influence due to upper TE
                xi, rho, k2, rj = calculate_xrm(te2, cp)
                vx = vortex_ring_vx(xi, rho, k2, rj, 19.5733 * rho)#lengths shouldn't be needed here, set such that self-induced case returns zero.
                vr = vortex_ring_vr(xi, rho, k2, rj)
                LHS[m, idxu] -= dot([vx; vr], nhat)
            end

        elseif verbose
            @warn "Trailing edge points for body $i are not coincident, skipping Kutta condtition."
        end
    end

    return nothing
end

######################################################################
#                                                                    #
#                            NEW RHS SETUP                           #
#                                                                    #
######################################################################
function update_RHS(
    b_bf::AbstractVector{T1}, A_bw, gamw::AbstractVector{T2}, A_br, sigr::AbstractMatrix{T3}
) where {T1,T2,T3}
    T = promote_type(T1, T2, T3)

    RHS = zeros(T, size(b_bf, 1))

    update_RHS!(RHS, b_bf, A_bw, gamw, A_br, sigr)

    return RHS
end

function update_RHS!(RHS, b_bf, A_bw, gamw, A_br, sigr)

    # start with freestream contributions to right hand side
    # note: the negative was already included in the precomputation for the freestream.
    RHS .+= b_bf

    # add wake vortex sheet contributions to right hand side
    # note: the subtraction was not included in the other coefficients in the precomputation, so we need to subract here
    RHS .-= A_bw * gamw

    # add rotor source sheet contributions to right hand side
    # note: the subtraction was not included in the other coefficients in the precomputation, so we need to subract here
    for jrotor in 1:length(A_br)
        # get induced velocity in the x-direction
        RHS .-= A_br[jrotor] * view(sigr, :, jrotor)
    end

    return nothing
end

######################################################################
#                                                                    #
#                      NEW LINEAR SOLVE (LSQ)                        #
#                                                                    #
######################################################################

#TODO; need to check that RHS isn't overwritten here, or if it is, need to make sure that doesn't mess things up elsewhere
"""
Given the original system of equations LHS*mu = RHS, it converts it into its
equivalent least-squares problem by prescribing the strengths of the
panels under `prescribedpanels` and returning Glsq and blsq.

NOTE: RHS is modified in place to become RHS-bp.
NOTES 2: Gred is an auxiliary matrix used to build Glsq and blsq.
"""
function prep_leastsquares!(Gred, Glsq, blsq, LHS, RHS, prescribedpanels)

    #=
        prescribedpanels = [
                                [prescribed_panel_index1, prescribed_strength1],
                                [prescribed_panel_index2, prescribed_strength2],
                                ...
                            ]
    =#

    # Total number of panels
    n = length(RHS)

    # Number of prescribed panels
    npres = length(prescribedpanels)

    # Error cases
    @assert size(LHS, 1) == n && size(LHS, 2) == n "" *
        "Invalid $(size(LHS, 1))x$(size(LHS, 2)) matrix LHS; expected $(n)x$(n)"
    @assert size(Gred, 1) == n && size(Gred, 2) == n - npres "" *
        "Invalid $(size(Gred, 1))x$(size(Gred, 2)) matrix Gred; expected $(n)x$(n-npres)"
    @assert size(Glsq, 1) == n - npres && size(Glsq, 2) == n - npres "" *
        "Invalid $(size(Glsq, 1))x$(size(Glsq, 2)) matrix Glsq; expected $(n-npres)x$(n-npres)"

    @assert length(RHS) == n "Invalid RHS length $(length(RHS)); expected $(n)"
    @assert length(blsq) == n - npres "Invalid blsq length $(length(blsq)); expected $(n-npres)"

    # Sort prescribed elements by index
    sort!(prescribedpanels; by=x -> x[1])

    # Move influence of prescribed panels to right-hand side
    for (paneli, strength) in prescribedpanels
        for i in 1:length(RHS)
            RHS[i] -= strength * LHS[i, paneli]
        end
    end

    # Reduce LHS: copy LHS into Gred without the prescribed panels
    prev_paneli = 0
    for (i, (paneli, strength)) in enumerate(prescribedpanels)
        Gred[:, (prev_paneli + 2 - i):(paneli - i)] .= view(
            LHS, :, (prev_paneli + 1):(paneli - 1)
        )

        if i == length(prescribedpanels) && paneli != size(LHS, 2)
            Gred[:, (paneli - i + 1):end] .= view(LHS, :, (paneli + 1):size(LHS, 2))
        end

        prev_paneli = paneli
    end

    tGred = transpose(Gred)

    # Store Gred'*(RHS - bp) under blsq
    mul!(blsq, tGred, RHS)

    # Store Gred'*Gred under Glsq
    mul!(Glsq, tGred, Gred)

    return Glsq, blsq
end

function prep_leastsquares(
    LHS::AbstractMatrix{T1},
    RHS::AbstractVector{T2},
    prescribedpanels::AbstractArray{Tuple{Int,T3}},
) where {T1,T2,T3}
    T = promote_type(T1, T2, T3)

    n = length(RHS)
    npres = length(prescribedpanels)

    LHSred = zeros(T, n, n - npres)
    LHSlsq = zeros(T, n - npres, n - npres)
    RHSlsq = zeros(T, n - npres)

    prep_leastsquares!(LHSred, LHSlsq, RHSlsq, LHS, RHS, prescribedpanels)

    return LHSlsq, RHSlsq
end

"""
Converts the reduced vector of strengths to the full vector of strengths
"""
function mured2mu!(mu, mured, prescribedpanels)

    # Total number of panels
    n = length(mu)

    # Number of prescribed panels
    npres = length(prescribedpanels)

    # Case of no prescrbied panels
    if npres == 0
        mu .= mured
        return mu
    end

    prev_paneli = 0

    # Iterate over prescribed panels building the full vector
    for (i, (paneli, strength)) in enumerate(prescribedpanels)
        mu[(prev_paneli + 1):(paneli - 1)] .= view(
            mured, (prev_paneli + 2 - i):(paneli - i)
        )
        mu[paneli] = strength

        if i == npres && paneli != n
            mu[(paneli + 1):end] .= view(mured, (paneli - i + 1):length(mured))
        end

        prev_paneli = paneli
    end

    return mu
end

function mured2mu(
    mured::AbstractVector{T1}, prescribedpanels::AbstractArray{Tuple{Int,T2}}
) where {T1,T2}
    T = promote_type(T1, T2)

    n = length(mured) + length(prescribedpanels)
    mu = zeros(T, n)

    mured2mu!(mu, mured, prescribedpanels)

    return mu
end

#---------------------------------#
#             Solvers             #
#---------------------------------#
function solve_body_strengths(
    LHS::AbstractMatrix{T1}, RHS::AbstractVector{T2}, prescribedpanels
) where {T1,T2}
    T = promote_type(T1, T2)

    mu = zeros(T, size(RHS, 1))
    solve_body_strengths!(mu, LHS, RHS, prescribedpanels)

    return mu
end

function solve_body_strengths!(mu, LHS, RHS, prescribedpanels)
    LHSlsq, RHSlsq = prep_leastsquares(LHS, RHS, prescribedpanels)

    mured = LHSlsq \ RHSlsq

    mured2mu!(mu, mured, prescribedpanels)

    return nothing
end

#TODO: for when cache is implemented
# function solve_body_strengths!(mu, body_system_matrices)

#     # - Extract Matrices - #
#     (; LHS, LHSred, LHSlsq, RHS, RHSlsq, prescribedpanels) = body_system_matrices

#     # - Prepping for Least Sqaures Solve - #
#     prep_leastsquares!(LHSred, LHSlsq, RHSlsq, LHS, RHS, prescribedpanels)

#     mured = LHSlsq \ RHSlsq

#     mured2mu!(mu, mured, prescribedpanels)

#     return nothing
# end

#---------------------------------#
#            Post Solve           #
#---------------------------------#
function update_stagnation_point_index!(mu, prescribedpanels; prescribedstrength=0.0)

    ## -- find the index at which the sign of mu changes -- ##

    # check that the sign actually changes
    if !all(sign.(mu) .== 1.0) && !all(sign.(mu) .== -1.0)
        # find the change in sign

        # set that index to be the prescribed panel index
        prescribedpanels[1] = (stagidx, prescribedstrength)
    end

    return nothing
end
