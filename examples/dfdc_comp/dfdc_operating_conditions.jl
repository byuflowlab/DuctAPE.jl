#=
OPER
!        Vinf         Vref          RPM
   0.0000       50.000       8000.0
!         Rho          Vso          Rmu           Alt
   1.2260       340.00      0.17800E-04   0.0000
!       XDwake        Nwake
  0.80000               20
!       Lwkrlx
            F
ENDOPER
=#

Vinf = 0.0 #TODO: probably want to change that...
Vref = 50.0 #TODO: look up how this is used.  is this why the Cp values are low in the dfdc example???
RPM = 8000
rho = 1.226
asound = 340.0
mu = 1.78e-5
wake_length = 0.8 #times duct chord?
nwake = 20 #20 wake sheets? TODO: look this up.
