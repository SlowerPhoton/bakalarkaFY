from solver import solve
from input_parser import parse_input_file, parse_table
from plot import plot
import numpy as np


#eV_to_K = 1.16045052e4
q_elem = 1.60217662e-19
k_B = 1.38064852e-23
#Td_to_Vcm2 = 1e-17
Td_to_Vm2 = 1e-21
Vm2_to_Td = 1 / Td_to_Vm2

bulk = 1e10

filename = "micro_cathode_test.input"
all_species, parameters, reactions, tables = parse_input_file(filename)

# constants
parameters['gas_density'] = parameters['gas_pressure'] / (k_B * parameters['gas_temperature'])
parameters['gap_area'] = np.pi * parameters['radius'] ** 2

########################## initialization ##########################
parameters['Ar'] = parameters['gas_density']  # Ar: 3.2188255028634424e+24
parameters['e'] = bulk
parameters['Ar^+'] = parameters['e']
parameters['EN'] = parameters['voltage'] / parameters['gap_length'] / parameters['gas_density'] * Vm2_to_Td

parameters["Ar*"] = 0

# equilibrium values
#parameters['EN'] = 50
#parameters['Ar^+'] = 1e9
#parameters['Ar2^+'] = 1e11
#parameters['Ar*'] = 1e11

E = parse_table("mean energy", tables)  # mean energy in eV

mobility = parse_table("mobility", tables)  # in 1/m


def Te(EN):
    # return 2 / 3 * E(EN) * q_elem / k_B  # scale from eV
    return E(EN) * 11_600

parameters['Te'] = Te(parameters['EN'])


def recombination_rate(prmtrs):
    return 8.75e-27 * (prmtrs["Te"] / 11600.0) ** (-4.5)
    #return 1e-37

reactions[1].rate_fun = recombination_rate


def update(prmtrs, time):
    neutral_particles = prmtrs["Ar"] + prmtrs["Ar*"]
    J = q_elem * prmtrs['gap_area'] * prmtrs['e'] * mobility(prmtrs['EN']) * prmtrs['EN'] * Td_to_Vm2
    prmtrs['EN'] = Vm2_to_Td * prmtrs['voltage'] / neutral_particles / \
                   (prmtrs['gap_length'] + prmtrs['resistance'] * J / (prmtrs['EN'] * Td_to_Vm2 * neutral_particles))
    prmtrs['EN'] = 0.5 * (prmtrs['EN'] + abs(prmtrs['EN']))
    parameters['Te'] = Te(parameters['EN'])
    return prmtrs

times, values = solve(all_species, parameters, reactions, update=update, recompute_N=False)
plot(times, values, all_species)
