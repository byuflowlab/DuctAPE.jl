from pyxdsm.XDSM import XDSM, OPT, SOLVER, FUNC, LEFT

# Change `use_sfmath` to False to use computer modern
x = XDSM(use_sfmath=True)


# - Blocks - #
x.add_system("precomp", FUNC, r"\text{Precomputations}")
x.add_system("solver", SOLVER, r"\text{Solver}")
x.add_system("body", FUNC, r"\text{Duct / Center Body}")
x.add_system("rotor", FUNC, r"\text{Rotor / Stator}")
x.add_system("wake", FUNC, r"\text{Wake}")


# - I/O - #
# x.add_input("solver", r"x^{(0)}")
x.add_output("solver", "x^*")
x.add_output("solver", "\mu^*,\Gamma^*,\sigma^*,\gamma^*")


# - Connections - #
# solver ->
x.connect("solver", "body", r"\mu, \Gamma, \sigma, \gamma")
x.connect("solver", "rotor", r"\mu, \Gamma, \sigma, \gamma")
x.connect("solver", "wake", r"\mu, \Gamma, \sigma, \gamma")

# body ->
x.connect("body", "solver", "\Delta \mu")
# x.connect("body", "rotor", "\mu^*")
# x.connect("body", "wake", "\mu^*")

# Rotor ->
# x.connect("rotor", "wake", "\Gamma^*, \sigma^*")
x.connect("rotor", "solver", "\Delta \Gamma, \Delta \sigma")

# Wake ->
x.connect("wake", "solver", "\Delta \gamma")


# Precomputations ->
x.connect("precomp", "body", [r"\text{LU Decomposed LHS}", r"\text{Freestream Boundary Condition}"])
x.connect("precomp", "rotor", [r"\text{Unit Induced Velocities}", r"\text{Interpolated Geometry}", r"\text{Operating Conditions}"])
x.connect("precomp", "wake", [r"\text{Unit Induced Velocities}",r"\text{Freestream Velocity}"])
x.connect("precomp", "solver", "\mu^{(0)},\Gamma^{(0)},\sigma^{(0)},\gamma^{(0)}")


# - Parameters - #
x.add_input("precomp", r"\text{Inputs}")





## -- Write it -- ##
x.write("dtxdsm")
