from pyxdsm.XDSM import XDSM, OPT, SOLVER, FUNC, LEFT

# Change `use_sfmath` to False to use computer modern
x = XDSM(use_sfmath=True)


# - Blocks - #
x.add_system("solve", SOLVER, r"\text{Non-Linear Solve}")
x.add_system("body", FUNC, r"\text{Body Aero}")
x.add_system("rotor", FUNC, r"\text{Rotor(s) Aero}")
x.add_system("wake", FUNC, r"\text{Wake Aero}")
x.add_system("post", FUNC, r"\text{Post-process}")


# - Solver I/O - #
x.add_input("solve", r"Geometry, Operating Conditions")
x.add_output("solve", r"\text{Aero Performance, Blade Loads}", side=LEFT)


# - Connections - #
# Body
x.connect("solve", "body", r"AIC's")
x.connect("body", "rotor", "\mathbf{V}_b")
x.connect("body", "wake", "\mathbf{V}_b")

# Rotors
x.connect("solve", "rotor", r"\Omega")
x.connect("rotor", "body", "\mathbf{V}_r")
x.connect("rotor", "wake", "\mathbf{V}_r, \Gamma")

# Wake
x.connect("solve", "wake", r"x, r")
x.connect("wake", "body", "\mathbf{V}_w")
x.connect("wake", "rotor", "\mathbf{V}_w")


# - Parameters - #
# x.add_input("body", r"")
# x.add_input("rotor", r"")
# x.add_input("wake", r"")


## -- Write it -- ##
x.write("dtmdo")
