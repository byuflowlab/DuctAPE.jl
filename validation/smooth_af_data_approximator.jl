using Roots
using FLOWMath

function xrotor_lift(
    clmax=2.0,
    amax=10.0,
    clmin=-1.5,
    amin=-10.0,
    alpha0=0.0,
    dclda=2.0 * pi,
    dclda_stall=-0.1,
    blend_hardness=5,
)
    function fmax(cl)
        return maximum(
            smooth_cl(cl, amax, clmin, amin, alpha0, dclda, dclda_stall, blend_hardness)[2]
        ) - clmax
    end

    function fmin(cl)
        return minimum(
            smooth_cl(clmax, amax, cl, amin, alpha0, dclda, dclda_stall, blend_hardness)[2]
        ) - clmin
    end

    trueclmax = Roots.find_zero(fmax, clmax)
    trueclmin = Roots.find_zero(fmin, clmin)

    return smooth_cl(
        trueclmax, amax, trueclmin, amin, alpha0, dclda, dclda_stall, blend_hardness
    )
end

"""
"""
function smooth_cl(
    clmax=2.0,
    amax=10.0,
    clmin=-1.5,
    amin=-10.0,
    alpha0=0.0,
    dclda=2.0 * pi,
    dclda_stall=-0.1,
    blend_hardness=5,
)

    #initialize range of alphas
    alpha = range(-pi, pi, 361)

    #get function for lift curve slope linear region
    cl_alpha = dclda * (alpha .- alpha0) .+ clmin

    #get function for negative stall region
    cl_ns = dclda_stall * (alpha .- alpha0) .+ clmin

    #get function for positive stall region.
    cl_ps = dclda_stall * (alpha .- alpha0) .+ clmax

    # use sigmoid flowmath function to put them together
    cl1 = FLOWMath.sigmoid_blend.(cl_ns, cl_alpha, alpha, amin, blend_hardness)
    cl = FLOWMath.sigmoid_blend.(cl1, cl_ps, alpha, amax, blend_hardness)

    return alpha, cl
end

"""

from xrotor docs:
                |                   2|           f
         CD  =  |CD  +  b (CL  - CL) | (Re/Re   )
                |  o         o       |       ref

where

   CD    =  minimum drag coefficient        ; Fortran name:  ( CDMIN )
     o

   CL    =  CL at which CD = CD                              ( CLDMIN )
     o                         o

   b     =  quadratic CD(CL) coefficient  d(CD)/d(CL**2)     ( DCDCL2 )


   Re    =  Reynolds Number at which CD formula applies      ( REREF )
     ref

   f     =  Reynolds Number scaling exponent.                ( REXP )
            Typically:
              f = -0.1 to -0.2  for high-Re turbulent flow Re > 2M
              f = -0.5 to -1.5  for "low-Re" regime Re ~ 200K..800K
              f = -0.3 to -0.5  for mostly-laminar airfoils at Re < 100K

"""
function xrotor_drag(cl, Re, cdo=5e-3, clo=0.15, b=4e-3, Re_ref=5e5, f=nothing)
    if isnothing(f)
        #=
        Typically:
              f = -0.1 to -0.2  for high-Re turbulent flow Re > 2M
              f = -0.5 to -1.5  for "low-Re" regime Re ~ 200K..800K
              f = -0.3 to -0.5  for mostly-laminar airfoils at Re < 100K
        =#
        #TODO: create smooth function for f(re)
        f = -1.0
    end

    cd = (cdo .+ b * (cl .- clo) .^ 2) * (Re / Re_ref)^f

    return cd
end

"""
"""
function generate_airfoil_data(lift_params, drag_params, Re; filename=nothing)
    # extract parameters
    (; clmax, clmin, alpha0, dclda, dclda_stall, blend_hardness) = lift_params
    (; cdo, clo, b, Re_ref, f) = drag_params

    alpha, cl = xrotor_lift(clmax, clmin, alpha0, dclda, dclda_stall, blend_hardness)

    cd = xrotor_drag(cl, Re, cdo, clo, b, Re_ref, f)

    # write file
    if !isnothing(filename)
        af = open(filename, "w")
        write(
            af,
            "Xrotor-like Airfoil, clmax=$clmax, clmin=$clmin, alpha0=$alpha0, dcdla=$dclda, dclda_stall=$dclda_stall, blend_hardness=$blend_hardness, cdo=$cdo, clo=$clo, b=$b, Re_ref=$Re_ref, f=$f\n",
        )
        write(af, "$Re\n")
        write(af, "0\n") #mach number

        for i in 1:length(alpha)
            write(af, "$(alpha[i]) $(cl[i]) $(cd[i])\n")
        end
        close(af)
    end

    return alpha, cl, cd
end
