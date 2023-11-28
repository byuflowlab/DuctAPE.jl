#=

compare intermediate ins/outs with dfdc's values

=#

project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

datapath = project_dir * "/dev_debug_archive/dfdc_comp/"
savepath = project_dir * "/dev_debug_archive/dfdc_comp/"

include(project_dir * "/plots_default.jl")

using DuctAPE
const dt = DuctAPE

# CCBlade used for it's airfoils function objects here.
using CCBlade
const ccb = CCBlade

using FLOWMath
const fm = FLOWMath

#---------------------------------#
# - Include DFDC DATA for Comparison
#---------------------------------#
# a bunch of the things to compare from dfdc,
include(project_dir * "/dev_debug_archive/dfdc_comp/ALL_THE_STUFF.jl")

#---------------------------------#
#  check velocity in, angles out  #
#---------------------------------#

Wtheta = dfdcwtheta
Wm = sqrt.(dfdcvx .^ 2 .+ dfdcvr .^ 2)

phi, alpha = dt.calculate_inflow_angles(Wm, Wtheta, dfdcbeta * pi / 180.0)

plot(; xlabel="Angles", ylabel="r")
plot!(phi * 180.0 / pi, dfdcr; label="inflow: DuctAPE functions")
plot!(alpha * 180.0 / pi, dfdcr; label="aoa: DuctAPE functions")
plot!(dfdcphi, dfdcr; linewidth=2, linestyle=:dash, label="inflow: DFDC")
plot!(dfdcalpha, dfdcr; linewidth=2, linestyle=:dash, label="aoa: DFDC")
savefig(project_dir * "/dev_debug_archive/dfdc_comp/inflow-io-comp.pdf")

#------------------------------------------#
# check geom&vel in, rey, solid, stagr out #
#------------------------------------------#

B = 5

alpha0 = 0.0
clmax = 1.5
clmin = -1.0
dclda = 2.0 * pi
dclda_stall = 0.5
dcl_stall = 0.2
cdmin = 0.012
cldmin = 0.1
dcdcl2 = 0.005
cmcon = 0.0
Re_ref = 2e5
Re_exp = 0.35
mcrit = 0.7

asound = 340.0
rho = 1.226
mu = 1.78e-5

afparams = (;
    alpha0,
    clmax,
    clmin,
    dclda,
    dclda_stall,
    dcl_stall,
    cdmin,
    cldmin,
    dcdcl2,
    cmcon,
    Re_ref,
    Re_exp,
    mcrit,
)

#CHORD,YRC,ALF,WWB,REY,SECSIG,SECSTAGR,CL,CD
include(datapath * "CLCD_CHECK.jl")

# local_reynolds = clcdcheck[:, 1] .* abs.(clcdcheck[:, 4]) * rho / mu
local_reynolds = dt.calc_reynolds.(clcdcheck[:, 1], clcdcheck[:, 4], rho, mu)
# local_solidity = B * clcdcheck[:, 1] ./ (2.0 * pi * clcdcheck[:, 2])
local_solidity = dt.get_local_solidity(B, clcdcheck[:, 1], clcdcheck[:, 2])
# local_stagger = 0.5 * pi .- dfdcbeta * pi / 180.0
local_stagger = dt.get_stagger(dfdcbeta * pi / 180.0)

plot(; xlabel="Reynolds", ylabel="r")
plot!(local_reynolds, clcdcheck[:, 2]; label="DuctAPE Functions")
plot!(clcdcheck[:, 5], clcdcheck[:, 2]; linewidth=2, linestyle=:dash, label="DFDC")
savefig(savepath * "re-io-comp.pdf")

plot(; xlabel="Local Solidity", ylabel="r")
plot!(local_solidity, clcdcheck[:, 2]; label="DuctAPE Functions")
plot!(clcdcheck[:, 6], clcdcheck[:, 2]; linewidth=2, linestyle=:dash, label="DFDC")
savefig(savepath * "solidity-io-comp.pdf")

plot(; xlabel="Local Stagger", ylabel="r")
plot!(local_stagger, clcdcheck[:, 2]; label="DuctAPE Functions")
plot!(clcdcheck[:, 7], clcdcheck[:, 2]; linewidth=2, linestyle=:dash, label="DFDC")
savefig(savepath * "stagger-io-comp.pdf")

#---------------------------------#
#  check angles in, cls, cds out  #
#---------------------------------#

cl = zeros(length(dfdcr))
cd = zeros(length(dfdcr))
cm = zeros(length(dfdcr))

for i in 1:length(dfdcr)
    cl[i], cd[i], cm[i] = dt.dfdc_clcdcm(
        clcdcheck[i, 4],
        clcdcheck[i, 5],
        clcdcheck[i, 6],
        clcdcheck[i, 7],
        clcdcheck[i, 3],
        afparams,
        asound,
    )
end

plot(; xlabel=L"c_\ell", ylabel="r")
plot!(cl, dfdcr; label="DuctAPE functions")
plot!(clcdcheck[:, 8], clcdcheck[:, 2]; linewidth=2, linestyle=:dash, label="DFDC")
savefig(savepath * "cl-io-comp.pdf")

plot(; xlabel=L"c_d", ylabel="r")
plot!(cd, clcdcheck[:, 2]; label="DuctAPE functions")
plot!(clcdcheck[:, 9], clcdcheck[:, 2]; linewidth=2, linestyle=:dash, label="DFDC")
savefig(savepath * "cd-io-comp.pdf")

#-----------------------------------#
# check cls/cds in, gamma/sigma out #
#-----------------------------------#

#sanity check from prints: Wm, chord, cd, sigr
dblchk = [
    58.5465393 8.91420543E-02 1.62045527E-02 3.36496867E-02
    64.9546051 7.97843859E-02 1.64908767E-02 3.39254625E-02
    72.1485596 7.12995827E-02 1.68367121E-02 3.43352258E-02
    79.9760818 6.39787391E-02 1.69718377E-02 3.46038528E-02
    88.2188873 5.77768683E-02 1.69386417E-02 3.45387086E-02
    96.7162323 5.25387935E-02 1.69188082E-02 3.43542881E-02
    105.385757 4.81022969E-02 1.69089083E-02 3.42206843E-02
    114.170273 4.43163998E-02 1.69053841E-02 3.41252387E-02
    123.028114 4.10622135E-02 1.69049483E-02 3.40549685E-02
    131.925644 3.82424332E-02 1.69075467E-02 3.40022333E-02
]

Gamr = similar(dfdccl) .= 0.0
sigr = similar(dfdccd) .= 0.0

for ir in 1:length(dfdcr)
    dt.gamma_sigma_from_coeffs!(
        view(Gamr, ir),
        view(sigr, ir),
        # dfdcWmag[ir],
        dblchk[ir, 1],
        B,
        # dfdcchord[ir],
        dblchk[ir, 2],
        dfdcr[ir],
        dfdccl[ir],
        # dfdccd[ir],
        dblchk[ir, 3],
    )
end

dcsigrcalc = similar(sigr) .= 0.0
for i in 1:length(dfdcr)
    if i == 1
        dcsigrcalc[i] = B / (4.0 * pi) .* dblchk[i, 1] .* dblchk[i, 2] .* dblchk[i, 3]
    else
        wave = (dblchk[i, 1] + dblchk[i - 1, 1]) / 2.0
        cave = (dblchk[i, 2] + dblchk[i - 1, 2]) / 2.0
        cdave = (dblchk[i, 3] + dblchk[i - 1, 3]) / 2.0
        dcsigrcalc[i] = B / (4.0 * pi) * wave * cave * cdave
    end
end
dcsigr = dblchk[:, 4]

plot(; xlabel="rotor circulation", ylabel="r")
plot!(Gamr, dfdcr; label="DuctAPE functions")
plot!(dfdcGamr, dfdcr; label="DFDC")
savefig(savepath * "rotorcirc-io-comp.pdf")

plot(; xlabel="rotor source strengths", ylabel="r")
plot!(sigr, dfdcr; label="DuctAPE functions")
# plot!(dfdcsigma, dfdcr; linewidth=2, linestyle=:dash, label="DFDC")
plot!(dcsigr, dfdcr; label="DFDC from prints")
plot!(dcsigrcalc, dfdcr; label="DuctAPE double check from prints")
savefig(savepath * "rotorsigma-io-comp.pdf")

#---------------------------------#
# check Gamr,sigr in, H and S out #
#---------------------------------#

Omega = 8000.0 * pi / 30.0
#get gamr in correct shape
dfdccirc = zeros(length(dfdcr), 1)
dfdccirc[:, 1] .= dfdcGamr

dS = dt.calculate_entropy_jumps(dfdcsigma, Wm)
dH = dt.calculate_enthalpy_jumps(dfdccirc, Omega, B)

plot(; xlabel="entropy jump across disk", ylabel="r")
plot!(dS, dfdcr; label="DuctAPE functions")
plot!(deltaS, dfdcr; linewidth=2, linestyle=:dash, label="DFDC")
savefig(savepath * "entropyjump-io-comp.pdf")

plot(; xlabel="enthalpy jump across disk", ylabel="r")
plot!(dH, dfdcr; label="DuctAPE functions")
plot!(deltaH, dfdcr; linewidth=2, linestyle=:dash, label="DFDC")
savefig(savepath * "enthalpyjump-io-comp.pdf")

#---------------------------------#
#  check everything in, gamw out  #
#---------------------------------#
