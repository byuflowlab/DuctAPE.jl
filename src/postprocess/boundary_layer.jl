using FLOWMath
using Plots
using DifferentialEquations

#---------------------------------#
#       Pre-solve Functions       #
#---------------------------------#

"""
"""
function sa1(alt; hardness=5)
    # Sigmoid blend of fits for <11000 and <25000
    T1(alt) = 15.04 - 0.00649 * alt + 273.1 #units: K
    P1(alt) = 101.29 * (T1(alt) / 288.08)^5.256 #units: kPa
    T2(alt) = -56.46 + 273.1 #units: K
    P2(alt) = 22.65 * exp(1.73 - 0.000157 * alt) #units: kPa

    return FLOWMath.sigmoid_blend(T1(alt), T2(alt), alt, 11000, hardness),
    FLOWMath.sigmoid_blend(P1(alt), P2(alt), alt, 11000, hardness)
end

"""
"""
function sa2(alt; hardness=5)
    # Sigmoid blend of <25000 and >25000
    T2(alt) = -56.46 + 273.1 #units: K
    P2(alt) = 22.65 * exp(1.73 - 0.000157 * alt) #units: kPa
    T3(alt) = -131.21 + 0.00299 * alt + 273.1 #units: K
    P3(alt) = 2.488 / ((T3(alt) / 216.6)^11.388) #units: kPa

    return FLOWMath.sigmoid_blend(T2(alt), T3(alt), alt, 25000, hardness),
    FLOWMath.sigmoid_blend(P2(alt), P3(alt), alt, 25000, hardness)
end

"""
"""
function ideal_gas_rho(P, T)
    return @. P / (0.2869 * T)
end

"""
    standard_atmosphere(alt; hardness=25)

Smoothed fits to the Standard Atmosphere model.

Assumes calorically imperfect gas.

# Arguments
- `alt::Float` : Altitude

# Keyword Arguments:
- `hardness::float` : hardness factor for sigmoid blend

# Returns
- `static_temperature::Float` : Static temperature
- `static_pressure::Float` : Static pressure
- `static_density::Float` : Static density
"""
function standard_atmosphere(alt; hardness=25)
    #Get temperature (T) and pressure (P) from table fits
    if alt < (11000 + 25000) / 2.0
        T, P = sa1(alt; hardness=hardness)
    else
        T, P = sa2(alt; hardness=hardness)
    end

    # return T, P, rho
    return T, P * 1000, ideal_gas_rho(P, T)
end

"""
"""
function speed_of_sound(static_pressure, static_density; gamma=1.4)
    return sqrt(gamma * static_pressure / static_density)
end

"""
"""
function calc_edge_mach(edge_velocity, speed_of_sound)
    return edge_velocity / speed_of_sound
end

"""
"""
function total_pressure(static_pressure, M; gamma=1.4)
    return static_pressure * (1.0 + (gamma - 1.0) * M^2 / 2.0)^(gamma / (gamma - 1.0))
end

"""
"""
function total_temperature(static_temperature, M; gamma=1.4)
    return static_temperature * (1.0 + (gamma - 1.0) * M^2 / 2.0)
end

# """
# """
# function t0_ts(minf; gamma=1.4)
#     return 1.0 + minf^2 * (gamma - 1.0) / 2.0
# end

# """
# assumes calorically perfect gas
# """
# function static_temperature(M, minf; gamma=1.4)
#     return t0_ts(minf; gamma=gamma) / t0_ts(M; gamma=gamma)
# end

"""
"""
function static_pressure(total_pressure, M; gamma=1.4)
    return total_pressure / ((1.0 + (gamma - 1.0) * M^2 / 2.0)^(gamma / (gamma - 1.0)))
end

"""
"""
function static_temperature(total_temperature, M; gamma=1.4)
    return total_temperature / (1.0 + (gamma - 1.0) * M^2 / 2.0)
end

"""
"""
function static_density(static_pressure, speed_of_sound; gamma=1.4)
    return gamma * static_pressure / speed_of_sound^2
end

"""
"""
function calc_edge_viscosity(
    static_temperure, mu_sea_level=0.0000181206, T_sea_level=288.15, S=110.4
)
    return mu_sea_level * (static_temperure / T_sea_level)^(3.0 / 2.0) * (T_sea_level + S) /
           (static_temperure + S)
end

#---------------------------------#
#         solve functions         #
#---------------------------------#
#=
note: function variables are in the following order:
- state variables
- variables dependent on state variables
- precomputed parameters
=#

"""
"""
function Fc(M)
    return sqrt(1.0 + 0.2 * M^2)
end

"""
"""
function FR(M)
    return 1.0 + 0.056 * M^2
end

"""
"""
function calc_Re(rhoe, Ue, d2, mue)
    return rhoe * Ue * d2 / mue
end

"""
"""
function calc_Cf0(Red2, M; hardness=5)
    logblend = FLOWMath.sigmoid_blend(
        1.05, log(10, FR(M) * Red2), FR(M) * Red2, 10^1.02, hardness
    )

    return (0.01013 / (logblend - 1.02) - 0.00075) / Fc(M)
end

"""
TODO: is there ever any way this could end up with a zero in the denominator?
"""
function calc_H12bar0(Cf0, M)
    return 1.0 / (1.0 - 6.55 * sqrt(Cf0 / 2 * (1.0 + 0.04 * M^2)))
end

"""
TODO: need to make sure that H12bar0 is never initialzied to zero
"""
function calc_Cf(H12bar, H12bar0, Cf0)
    return Cf0 * (0.9 / (H12bar / H12bar0 - 0.4) - 0.5)
end

"""
"""
function calc_H12(H12bar, M, Pr=1.0)
    return (H12bar + 1.0) * (1.0 + Pr^(1.0 / 3.0) * M^2 / 5) - 1.0
end

"""
"""
function calc_Ctau(CE, Cf0, M)
    return (0.024 * CE + 1.2 * CE^2 + 0.32Cf0) * (1.0 + 0.1 * M^2)
end

"""
TODO: CE cannot be -0.01 or we'll have a divide by zero.
Green talks about handling this.
"""
function calc_F(CE, Cf0)
    return (0.02 * CE + CE^2 + 0.8 * Cf0 / 3) / (0.01 + CE)
end

"""
TODO: H12bar can't be 1, or we divide by zero
can H12bar ever get that low?
"""
function calc_H1(H12bar)
    return 3.15 + 1.72 / (H12bar - 1.0) - 0.01 * (H12bar - 1.0)^2
end

"""
TODO: at what value(s) of H12bar does the denominator go to zero?
"""
function calc_dH12bardH1(H12bar)
    return -(H12bar - 1.0)^2 / (1.72 + 0.02 * (H12bar - 1.0)^3)
end

"""
TODO: can't evaluate this when r=0 or we'll divide by zero
"""
function calc_richardson_number(H12bar, d2, H12, H1, R)
    return 2.0 * d2 / (3.0 * R) * (H12 + H1) * (H1 / H12bar + 0.3)
end

"""
    beta(Ri; hardness=100.0)

Sigmoind blended version of piecewise function:
    | 7.0 if Ri > 0
β = |
    | 4.5 if Ri < 0

# Arguments:
- `Ri::float` : Richardson number

# Keyword Arguments:
- `hardness::float` : hardness factor for sigmoid blend

# Returns:
- `beta::float` : factor used in secondary influence from longitudinal curvature.
"""
function beta(Ri, hardness=100.0)
    return flowmath.sigmoid_blend(4.5, 7.0, Ri, 0.0, hardness)
end

"""
"""
function longitudinal_curvature_influence(M, Ri)
    return 1 + beta(Ri) * (1.0 + M^2 / 5) * Ri
end

"""
"""
function lateral_strain_influence(H12bar, d2, H12, H1, r, drds)
    return 1.0 - 7.0 / 3.0 * (H1 / H12bar + 1.0) * (H12 + H1) * d2 * drds / r
end

"""
TODO: Ue can't be zero or we divide by zero
"""
function dilation_influence(H12bar, d2, H12, H1, M, Ue, dUedx)
    return 1.0 + 7.0 / 3.0 * M^2 * (H1 / H12bar + 1.0) * (H12 + H1) * d2 * dUedx / Ue
end

"""
"""
function calc_lambda(args...)
    return *(args...)
end

"""
TODO; make sure H12 can't be zero
"""
function calc_d2dUedxUeeq0(H12bar, H12, Cf, M)
    return 1.25 * (Cf / 2.0 - ((H12bar - 1.0) / (6.432 * H12bar))^2 / (1.0 + 0.04 * M^2)) /
           H12
end

"""
"""
function calc_CEeq0(H1, H12, Cf, d2dUedxUeeq0)
    return H1 * (Cf / 2.0 - (H12 + 1.0) * d2dUedxUeeq0)
end

"""
"""
function calc_Ctaueq0(CEeq0, Cf0, M)
    return (0.24 * CEeq0 + 1.2 * CEeq0^2 + 0.32 * Cf0) * (1.0 + 0.1 * M^2)
end

"""
TODO: c can be negative here.
"""
function calc_CEeq(Ctaueq0, M, lambda, Cf0)
    c = Ctaueq0 / ((1.0 + 0.1 * M^2) * lambda^2) - 0.32 * Cf0
    return sqrt(c / 1.2 + 0.0001) - 0.01
end

"""
TODO: make sure H12 can't be -1
"""
function calc_d2dUedxUeeq(H1, H12, Cf, CEeq)
    return (Cf / 2.0 - CEeq / H1) / (H12 + 1.0)
end


#---------------------------------#
#          Viscous Drag           #
#---------------------------------#

"""
"""
function squire_young(d2,Ue, Uinf, H12, c)
    return 2.0*d2/c *(Ue/Uinf)^((5.0+H12)/2.0)
end

#---------------------------------#
#  state initilization functions  #
#---------------------------------#

"""
initialize momentum thickness state with flat plate schlichting solution given a dx from the stagnation point to the next panel edge.
"""
function d2_init(dx, Rex)
    return 0.036 * dx / Rex^0.2
end

"""
todo: initialize with H12bar0
"""
function H12bar_init(Cf0, M)
    return calc_H12bar0(Cf0, M)
end

"""
todo: initialize with CEeq
"""
function CE_init(Ctaueq0, M, Cf0)
    return calc_CEeq(Ctaueq0, M, 1, Cf0)
end

"""
"""
function f(u, params, s)
    du = zeros(length(u))
    return f!(du, u, params, s)
end

"""
"""
function f!(du, u, params, s)

    # - unpack variables - #
    rd2, H12bar, CE = u
    (;
        r_coords,
        edge_velocity,
        edge_density,
        edge_viscosity,
        edge_mach,
        edge_acceleration,
        surface_curvature,
        surface_derivative,
    ) = params

    # - rename for convenience - #
    r = r_coords(s)
    Ue = edge_velocity(s)
    rhoe = edge_density(s)
    mue = edge_viscosity(s)
    M = edge_mach(s)
    dUedx = edge_acceleration(s)
    R = surface_curvature(s)
    drds = surface_derivative(s)

    # - calculate intermediate variables - #
    d2 = rd2 / r
    Red2 = calc_Re(d2, Ue, rhoe, mue)
    Cf0 = calc_Cf0(Red2, M)
    H12bar0 = calc_H12bar0(Cf0, M)
    Cf = calc_Cf(H12bar, H12bar0, Cf0)
    H12 = calc_H12(H12bar, M)

    dH12bardH1 = calc_dH12bardH1(H12bar) #yes, this should be negative according to example
    H1 = calc_H1(H12bar)

    F = calc_F(CE, Cf0)

    d2dUedxUeeq0 = calc_d2dUedxUeeq0(H12bar, H12, Cf, M)
    CEeq0 = calc_CEeq0(H1, H12, Cf, d2dUedxUeeq0)
    Ctaueq0 = calc_Ctaueq0(CEeq0, Cf0, M)

    Ctau = calc_Ctau(CE, Cf0, M)

    if params.lambda
        if params.longitudinal_curvature
            Ri = calc_richardson_number(H12bar, d2, H12, H1, r)
            l1 = longitudinal_curvature_influence(M, Ri)
        else
            l1 = 1
        end

        if params.lateral_strain
            l2 = lateral_strain_influence(H12bar, d2, H12, H1, r, drds)
        else
            l2 = 1
        end

        if params.dilation
            l3 = dilation_influence(H12bar, d2, H12, H1, M, Ue, dUedx)
        else
            l3 = 1
        end

        lambda = calc_lambda(l1, l2, l3)
    else
        lambda = 1
    end
    CEeq = calc_CEeq(Ctaueq0, M, lambda, Cf0)
    d2dUedxUeeq = calc_d2dUedxUeeq(H1, H12, Cf, CEeq)

    # - system of equations - #

    # momentum
    du[1] = r * Cf / 2.0 - (H12 + 2.0 - M^2) * rd2 * dUedx / Ue

    # entrainment
    du[2] = dH12bardH1 * (CE - H1 * (Cf / 2.0 - (H12 + 1.0) * d2 * dUedx / Ue)) / d2

    # lag
    du[3] =
        F / d2 * (
            2.8 / (H12 + H1) * (sqrt(Ctaueq0) - lambda * sqrt(Ctau)) + d2dUedxUeeq -
            d2 * dUedx / Ue * (1.0 + 0.075 * M^2 * ((1.0 + 0.2 * M^2) / (1.0 + 0.1 * M^2)))
        )

    return du, Cf, H12
end

"""
"""
function initialize_turbulent_boundary_layer_states(r_init, rhoe, Ue, dx, mue, M)

    # - Initialize States - #
    # initialize momentum thickness (d2) using flat plate schlichting model
    Rex = calc_Re(rhoe, Ue, dx, mue)
    d2 = d2_init(dx, Rex)

    # initialize H12bar (compressible shape factor) using _0 equations
    Red2 = calc_Re(rhoe, Ue, d2, mue)
    Cf0 = calc_Cf0(Red2, M)
    H12bar0 = H12bar_init(Cf0, M)

    # initialize the entrainment coefficient (CE) using equilibrium equations
    H1 = calc_H1(H12bar0)
    H12 = calc_H12(H12bar0, M)
    Cf = calc_Cf(H12bar0, H12bar0, Cf0)
    d2dUedxUeeq0 = calc_d2dUedxUeeq0(H12bar0, H12, Cf, M)
    CEeq0 = calc_CEeq0(H1, H12, Cf, d2dUedxUeeq0)
    Ctaueq0 = calc_Ctaueq0(CEeq0, Cf0, M)
    CE = CE_init(Ctaueq0, M, Cf0)

    # initial states
    return [r_init * d2, H12bar0, CE], Cf, H12
end

"""
"""
function RK2(f, u, params, s, ds)
    k1, _, _ = f(u, params, s)
    k2, Cf, H12 = f(u .+ (ds / 2) .* k1, params, s + (ds / 2))
    unext =@. u + k2 * ds
    return unext, Cf, H12
end

"""
"""
function RK4(f, u, params, s, ds)
    k1, Cf1, H12_1 = f(u, params, s)
    k2, Cf2, H12_2 = f(u .+ (ds / 2) .* k1, params, s + (ds / 2))
    k3, Cf3, H12_3 = f(u .+ (ds / 2) .* k2, params, s + (ds / 2))
    k4, Cf4, H12_4 = f(u .+ ds .* k3, params, s + ds)
    unext = @. u + (k1 + k2 * 2 + k3 * 2 + k4) * ds / 6
    Cfnext = (Cf1 + Cf2 * 2 + Cf3 * 2 + Cf4) / 6
    H12next = (H12_1 + H12_2 * 2 + H12_3 * 2 + H12_4) / 6
    return unext, Cfnext, H12next
end

"""
"""
function solve_turbulent_boundary_layer_rk4!(f, init, params, svec; s_init_default=0.01)
    u0, Cf0, H12_0 = init
    s = svec[1] > s_init_default ? [svec[1]] : [s_init_default]
    sep = [false]
    sepid = [1]
    us = zeros(eltype(u0), length(u0), length(svec))
    Cfs = zeros(eltype(u0), length(svec))
    H12s = zeros(eltype(u0), length(svec))
    us[:, 1] = u0
    Cfs[1] = Cf0
    H12s[1] = H12_0

    for i in 1:(length(svec) - 1)

        # take step
        us[:, i + 1], Cfs[i + 1], H12s[i+1] = RK4(f, us[:, i], params, svec[i], svec[i + 1] - svec[i])

        if Cfs[i + 1] <= 0.0
            sep[1] = true
            sepid[1] = i
            break
        end
    end

    if sep[1] == true
        u1sep = FLOWMath.akima(Cfs[(sepid - 1):sepid], us[1, (sepid - 1):sepid], 0.0)
        u2sep = FLOWMath.akima(Cfs[(sepid - 1):sepid], us[2, (sepid - 1):sepid], 0.0)
        u3sep = FLOWMath.akima(Cfs[(sepid - 1):sepid], us[3, (sepid - 1):sepid], 0.0)
        usep = [u1sep; u2sep; u3sep]
        H12sep = FLOWMath.akima(Cfs[(sepid - 1):sepid], H12s[(sepid - 1):sepid], 0.0)
        ssep = FLOWMath.akima(Cfs[(sepid - 1):sepid], svec[(sepid - 1):sepid], 0.0)
    else
        usep = us[:, end]
        H12sep = H12s[end]
        ssep = svec[end]
    end

    # return states at separate, and separation state, and separation postition
    return us, usep, Cfs, H12sep, ssep, sepid[1]
end

"""
"""
function solve_turbulent_boundary_layer(f!, u0, xspan, params)
    prob = DifferentialEquations.ODEProblem(f!, u0, xspan, params)
    # alg = DifferentialEquations.radau()
    alg = DifferentialEquations.Tsit5()
    sol = DifferentialEquations.solve(prob, alg)
    return sol
end

"""
"""
function test_flat_plate(s0, L; ds=1e-2)

    # - Parameters - #
    lambda = false
    longitudinal_curvature = false
    lateral_strain = false
    dilation = false
    s = range(s0 + ds, s0 + L; step=ds)
    N = length(s)
    r_coords = ones(N)
    drds = zeros(N)
    Uinf = 10.0
    edge_velocity = Uinf * ones(N)
    edge_acceleration = zeros(N)
    surface_curvature = zeros(N)
    surface_derivative = zeros(N)
    alt = 0.0

    Tinf, Pinf, rhoinf = standard_atmosphere(alt; hardness=25)

    a = speed_of_sound(Pinf, rhoinf)

    edge_mach = calc_edge_mach(edge_velocity, a)
    Minf = calc_edge_mach(Uinf, a)

    P0 = total_pressure(Pinf, Minf)

    T0 = total_temperature(Tinf, Minf)

    Pe = static_pressure.(Ref(P0), edge_mach)

    Te = static_temperature.(Ref(T0), edge_mach)

    edge_density = static_density(Pe, a)

    edge_viscosity = calc_edge_viscosity.(Te)

    # - Spline things in terms of curve length - #
    params = (;
        lambda,
        longitudinal_curvature,
        lateral_strain,
        dilation,
        edge_velocity=FLOWMath.Akima(s, edge_velocity),
        edge_mach=FLOWMath.Akima(s, edge_mach),
        edge_density=FLOWMath.Akima(s, edge_density),
        edge_viscosity=FLOWMath.Akima(s, edge_viscosity),
        edge_acceleration=FLOWMath.Akima(s, edge_acceleration),
        surface_curvature=FLOWMath.Akima(s, surface_curvature),
        surface_derivative=FLOWMath.Akima(s, surface_derivative),
        r_coords=FLOWMath.Akima(s, r_coords),
        Uinf,
    )

    u0, Cf0, H12_0 = initialize_turbulent_boundary_layer_states(
        r_coords[2],
        edge_density[1],
        edge_velocity[1],
        s[1],
        edge_viscosity[1],
        edge_mach[1],
    )

    # return solve_turbulent_boundary_layer(f!, u0, (s[1], s[end]), params)
    return u0, Cf0, H12_0, params, s
end

ds = 1e-4
s0 = 0.01
L = 1.0
u0, Cf0, H12_0, params, svec = test_flat_plate(s0, L; ds=ds)

"""
"""
function schlichting(x, Ue, rhoe, mue)
    Rex = rhoe * Ue * x / mue
    d1 = 0.046 * x / Rex^0.2
    d2 = 0.036 * x / Rex^0.2
    cf = 0.0592 / Rex^0.2
    return d1, d2, cf
end

schs = svec
sch =
    schlichting.(
        schs,
        params.edge_velocity.(schs),
        params.edge_density.(schs),
        params.edge_viscosity.(schs),
    )

sd1 = getindex.(sch, 1)
sd2 = getindex.(sch, 2)
scf = getindex.(sch, 3)

plotstep = Int(1 / ds * L * 5 * 1e-3)
p1 = plot(; xlabel="s", ylabel="momentum thickness")
plot!(
    p1,
    schs[1:plotstep:end],
    sd2[1:plotstep:end];
    label="Schlicting",
    linewidth=3,
    linestyle=:dot,
)

us, usep, Cfs, H12sep, ssep, separation_id = solve_turbulent_boundary_layer_rk4!(
    f, [u0, Cf0, H12_0], params, svec; s_init_default=0.01
)

cd = squire_young(usep[1],params.edge_velocity(ssep), params.Uinf, H12sep, L)

plot!(p1, svec[1:plotstep:end], us[1, 1:plotstep:end]; label="Green RK4")
savefig("mom_thick_schlichting_green.png")

# pcf = plot(xlabel="s", ylabel="Cf")
# plot!(pcf, svec[1:plotstep:end], Cfs[1:plotstep:end]; label="Green RK4")
