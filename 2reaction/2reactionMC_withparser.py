import matplotlib.pyplot as plt

from solver import solve_withN
from input_parser import parse_input_file

filename = "2reaction_2N.input"
all_species, parameters, reactions, tables = parse_input_file(filename)
parameters['e'] = parameters['Ar^+'] = 1e5
parameters['Ar'] = 1e10

for N in [100, 500, 1000, 5000]:
    parameters_copy = parameters.copy()
    parameters_copy['N'] = N  # sets N manually, so there is no need to parse many input files, differing only in N
    times, values = solve_withN(all_species, parameters_copy, reactions)
    plt.plot(times, values['e'], label=f"N = {N}")

plt.xlabel("time [s]")
plt.ylabel("electron concentration")
plt.title(f"initial electron concentration = {parameters['e']}")
plt.legend()
plt.show()
