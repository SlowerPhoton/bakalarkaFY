from solver import solve
from input_parser import parse_input_file, parse_table
from plot import plot
import numpy as np


eV_to_K = 1.16045052e4
q_elem = 1.60217662e-19
k_B = 1.38064852e-23
Td_to_Vcm2 = 1e-17
Td_to_Vm2 = 1e-21
Vm2_to_Td = 1 / Td_to_Vm2

bulk = 1e18

filename = "micro_cathode.input"
all_species, parameters, reactions, tables = parse_input_file(filename)

# constants
parameters['gas_density'] = parameters['gas_pressure'] / (k_B * parameters['gas_temperature'])
# in ZDPlasKin: 3.2188309472845553e+18
parameters['gap_area'] = np.pi * parameters['radius'] ** 2

########################## initialization ##########################
parameters['Ar'] = parameters['gas_density']  # Ar: 3.2188255028634424e+24
parameters['e'] = bulk
parameters['Ar^+'] = parameters['e']
parameters['EN'] = parameters['voltage'] / parameters['gap_length'] / parameters['gas_density'] * Vm2_to_Td

E = parse_table("mean energy", tables)  # mean energy in eV


def Te(EN):
    return 2 / 3 * E(EN) * q_elem / k_B  # scale from eV


parameters['Te'] = Te(parameters['EN'])

mobility = parse_table("mobility", tables)  # in 1/m


########################## rate function initialization ##########################

def reaction4_rate(prmtrs):
    # in m^3/s
    return 8.5e-13 * (prmtrs["Te"] / 300.0) ** (-0.67)


reactions[4].rate_fun = reaction4_rate


def reaction5_rate(prmtrs):
    return 6.06e-12 / prmtrs['gas_temperature'] * np.exp(-15130.0 / prmtrs['gas_temperature'])


reactions[5].rate_fun = reaction5_rate


def reaction7_rate(prmtrs):
    # TODO: remove when ready
    # return 8.75e-27 * (prmtrs["Te"] / 11600.0) ** (-4.5)
    return 1e-37


reactions[7].rate_fun = reaction7_rate


def reaction9_rate(prmtrs):
    return 2.25e-43 * (prmtrs['gas_temperature'] / 300.0) ** (-0.4)


reactions[9].rate_fun = reaction9_rate


def reaction_surface_rate(prmtrs):
    # TODO: definitely needs polishing up the units
    return 1.52 * (760 / prmtrs['gas_pressure']) * (prmtrs['gas_temperature'] / 273.16) * \
           (prmtrs['Te'] / 11600) * ((2.405 / prmtrs['radius']) ** 2 + (3.141 / prmtrs['gap_length']) ** 2)


reactions[10].rate_fun = reaction_surface_rate
reactions[11].rate_fun = reaction_surface_rate
reactions[12].rate_fun = reaction_surface_rate

# TODO: remove when ready
reactions = [reactions[0], reactions[7]]

# scale the rate functions
for reaction in reactions:
    # correctly parsed reaction rates are already set: we overwrite them in order to not scale them twice
    reaction.scale(parameters, overwrite_ratio=True)


#############################
# DYNAMIC
#############################
def update(prmtrs):
    neutral_particles = prmtrs["Ar"] + prmtrs["Ar*"]
    J = q_elem * prmtrs['gap_area'] * prmtrs['e'] * mobility(prmtrs['EN']) * prmtrs['EN'] * Td_to_Vm2
    # TODO: remove when done
    J = 0
    prmtrs['EN'] = prmtrs['voltage'] / neutral_particles * Vm2_to_Td / \
                   (prmtrs['gap_length'] +
                    prmtrs['resistance'] * J / (prmtrs['EN'] * Td_to_Vm2 * neutral_particles + 1.0e-99))
    prmtrs['EN'] = 0.5 * (prmtrs['EN'] + abs(prmtrs['EN']))
    parameters['Te'] = Te(parameters['EN'])
    return prmtrs


times, values = solve(all_species, parameters, reactions, update=update, bulk=bulk)
plot(times, values, all_species)
