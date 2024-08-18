using Pkg
Pkg.add("PyCall")
Pkg.add("Conda")
pyimport_conda("scipy.integrate", "scipy")

using PyCall
# Import scipy's solve_ivp
solve_ivp = pyimport("scipy.integrate").solve_ivp