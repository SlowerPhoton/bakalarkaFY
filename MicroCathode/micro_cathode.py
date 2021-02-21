from solver import solve
from input_parser import parse_file
from plot import plot

eV_to_K = 1.16045052e4
q_elem = 1.60217662e-19
k_B = 1.38064852e-23

filename = "micro_cathode.input"
all_species, parameters, reactions = parse_file(filename)

parameters['gap_area'] = 3.1415 * parameters['radius'] ** 2
parameters['gas_density'] = 1.0e-6 * (101325.0 * parameters['gas_pressure']/760.0) \
                            / (1.38065e-23 * parameters['gas_temperature'])
parameters['Ar'] = parameters['gas_density']
parameters['Ar^+'] = parameters['e']
parameters['EN'] = parameters['voltage'] / parameters['gap_length'] / parameters['gas_density'] * 1.0e17

#############################
#DYNAMIC
#############################
def update(parameters_old):
    return parameters_old

J  = 1.6d-19 * gap_area * density(species_electrons) * Vdr
V  = voltage - resistance * J

times, values = solve(all_species, parameters, reactions, update=update)
plot(times, values, all_species)
