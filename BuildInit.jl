using Pkg
Pkg.add("PyCall")
Pkg.add("Conda")
Pkg.add("Plots")
Pkg.add("PlotlyJS")

Pkg.build("PyCall")
Pkg.build("Conda")

using PyCall
using Conda

Conda.add("scipy")

# Ensure scipy and its modules are imported correctly
scipy_integrate = pyimport("scipy.integrate")
scipy_constants = pyimport("scipy.constants")

# Access solve_ivp directly
solve_ivp = scipy_integrate.solve_ivp
