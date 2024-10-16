using Pkg
Pkg.add("PyCall")
Pkg.add("Conda")
Pkg.add("Plots")
Pkg.add("PlotlyJS")
Pkg.add("Evolutionary")
Pkg.add("LinearAlgebra")

Pkg.build("PyCall")
Pkg.build("LinearAlgebra")
Pkg.build("Conda")
Pkg.build("Plots")
Pkg.build("PlotlyJS")
Pkg.build("Evolutionary")

using PyCall
using Conda

Conda.add("scipy")

# Ensure scipy and its modules are imported correctly
scipy_integrate = pyimport("scipy.integrate")
scipy_constants = pyimport("scipy.constants")

# Access solve_ivp directly
solve_ivp = scipy_integrate.solve_ivp
