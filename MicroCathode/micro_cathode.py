from solver import solve
from input_parser import parse_input_file, parse_table
from plot import plot
import numpy as np

eV_to_K = 1.16045052e4
q_elem = 1.60217662e-19
k_B = 1.38064852e-23

filename = "micro_cathode.input"
all_species, parameters, reactions, tables = parse_input_file(filename)

cm_to_micro = 1e4
m_to_micro = 1e6

########################## rescaling parameters ##########################
parameters['radius'] *= cm_to_micro
parameters['e'] *= cm_to_micro ** -3
parameters['gap_length'] *= cm_to_micro

# constants (originally in cm)
parameters['gas_density'] = 1.0e-6 * (101325.0 * parameters['gas_pressure'] / 760.0) \
                            / (1.38065e-23 * parameters['gas_temperature']) * cm_to_micro ** -3
parameters['gap_area'] = 3.1415 * parameters['radius'] ** 2

########################## initialization ##########################
parameters['Ar'] = parameters['gas_density']
parameters['Ar^+'] = parameters['e']
parameters['EN'] = parameters['voltage'] / parameters['gap_length'] / parameters['gas_density'] \
                   * 1.0e17 * cm_to_micro ** -2  # this scales the townsend

E = parse_table("mean energy", tables)


def Te(EN):
    return 2 / 3 * E(EN) * q_elem / k_B


parameters['Te'] = Te(parameters['EN'])

mobility = parse_table("mobility", tables)  # in 1/m

########################## rate function initialization ##########################

def reaction5_rate(prmtrs):
    return 8.5e-7 * (prmtrs["Te"] / 300.0) ** (-0.67)


reactions[4].rate_fun = reaction5_rate


def reaction6_rate(prmtrs):
    return 6.06e-6 / prmtrs['gas_temperature'] * np.exp(-15130.0 / prmtrs['gas_temperature'])


reactions[5].rate_fun = reaction6_rate


def reaction8_rate(prmtrs):
    return 8.75e-27 * (prmtrs["Te"] / 11600.0) ** (-4.5)


reactions[7].rate_fun = reaction8_rate


def reaction10_rate(prmtrs):
    return 2.25e-31 * (prmtrs['gas_temperature'] / 300.0) ** (-0.4)


reactions[9].rate_fun = reaction10_rate


def reaction_surface_rate(prmtrs):
    return 1.52 * (760 / prmtrs['gas_pressure']) * (prmtrs['gas_temperature'] / 273.16) * (
            prmtrs['Te'] / 11600) * \
           ((2.405 / prmtrs['radius']) ** 2 + (3.141 / prmtrs['gap_length']) ** 2)


reactions[10].rate_fun = reaction_surface_rate
reactions[11].rate_fun = reaction_surface_rate
reactions[12].rate_fun = reaction_surface_rate

# scale the rate functions
for reaction in reactions:
    reaction.scale(parameters)


#############################
# DYNAMIC
#############################
def update(prmtrs):
    neutral_particles = prmtrs["Ar"] + prmtrs["Ar*"]
    J = 1.6e-19 * prmtrs['gap_area'] * prmtrs['e'] * mobility(prmtrs['EN']) * 1 / m_to_micro * prmtrs['EN'] * 1e-17
    prmtrs['EN'] = prmtrs['voltage'] / (prmtrs['gap_length'] + prmtrs['resistance'] * J
                                                / (prmtrs['EN'] * neutral_particles / 1.0e17 + 1.0e-99)) \
                       / neutral_particles * 1.0e17
    prmtrs['EN'] = 0.5 * (prmtrs['EN'] + abs(prmtrs['EN']))
    return prmtrs


# V  = voltage - resistance * J

times, values = solve(all_species, parameters, reactions, update=update)
plot(times, values, all_species)
