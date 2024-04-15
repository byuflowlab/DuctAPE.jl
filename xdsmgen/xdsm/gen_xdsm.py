from pyxdsm.XDSM import XDSM, OPT, SOLVER, FUNC, LEFT

# Change `use_sfmath` to False to use computer modern
x = XDSM(use_sfmath=True)


# - Blocks - #
x.add_system("pre", FUNC, r"\text{Pre-computation}")
x.add_system("solve", SOLVER, r"\text{Non-Linear Solve}")
x.add_system("body", FUNC, r"\text{Body Aero}")
x.add_system("rotor", FUNC, r"\text{Rotor(s) Aero}")
x.add_system("wake", FUNC, r"\text{Wake Aero}")
x.add_system("post", FUNC, r"\text{Post-process}")


# - Solver I/O - #
x.add_input("pre", r"\text{Geometry, Operating Conditions}")
x.add_output("post", r"\text{Aero Performance, Blade Loads}", side=LEFT)


# - Connections - #
# Body
x.connect("solve", "body", "\gamma_b, \gamma_w, \sigma_r")
# x.connect("body", "rotor", "\mathbf{V}_b")
# x.connect("body", "wake", "\mathbf{V}_b")
x.connect("body", "post", "\gamma_b^*")
x.connect("body", "solve", "\gamma_b")
x.connect("body", "solve", "\gamma_b^*")

# Rotors
x.connect("solve", "rotor", "\gamma_b, \Gamma_r, \sigma_r, \gamma_w")
# x.connect("rotor", "body", "\mathbf{V}_r")
# x.connect("rotor", "wake", "\mathbf{V}_r, \Gamma_r^*")
x.connect("rotor", "post", "\Gamma_r^*, \sigma_r^*")
x.connect("rotor", "solve", "\Gamma_r^*, \sigma_r^*")

# Wake
x.connect("solve", "wake", "\gamma_b, \Gamma_r, \sigma_r, \gamma_w")
# x.connect("wake", "body", "\mathbf{V}_w")
# x.connect("wake", "rotor", "\mathbf{V}_w")
x.connect("wake", "post", "\gamma_w^*")
x.connect("wake", "solve", "\gamma_w^*")

# Post-process
x.connect("pre", "post", r"\text{Geometry, Operating Condtions}")

# - Parameters - #
x.connect("pre","body", r"L, U, \text{AIC's}")
x.connect("pre","rotor", r"\text{Blade Elements}")
# x.add_input("rotor", "\Omega")
x.connect("pre","wake", r"\text{``Grid''}")
x.connect("pre","solve", r"\gamma_b^{(o)},\Gamma_r^{(o)},\sigma_r^{(o)},\gamma_w^{(o)}")


## -- Write it -- ##
x.write("dtsolve")
