from solver import solve
from input_parser import parse_input_file, parse_tables, parse_table, parse_table_custom, get_table_scales
from plot import plot
import numpy as np

eV_to_K = 1.16045052e4
q_elem = 1.60217662e-19
k_B = 1.38064852e-23

filename = "micro_cathode.input"
all_species, parameters, reactions = parse_input_file(filename)

cm_to_micro = 1e4

# rescale
parameters['radius'] *= cm_to_micro
parameters['e'] *= cm_to_micro ** -3
parameters['gap_length'] *= cm_to_micro

# constants (originally in cm)
parameters['gas_density'] = 1.0e-6 * (101325.0 * parameters['gas_pressure'] / 760.0) \
                            / (1.38065e-23 * parameters['gas_temperature']) * cm_to_micro ** -3
parameters['gap_area'] = 3.1415 * parameters['radius'] ** 2

# initialization
parameters['Ar'] = parameters['gas_density']
parameters['Ar^+'] = parameters['e']
parameters['EN'] = parameters['voltage'] / parameters['gap_length'] / parameters['gas_density'] \
                   * 1.0e17 * cm_to_micro ** -2  # this scales the townsend

tables = parse_tables(parameters)
E = parse_table_custom("mean energy")
table_scales = get_table_scales(parameters, reactions, tables)

def Te(EN):
    return 2 / 3 * E(EN) * q_elem / k_B


parameters['Te'] = Te(parameters['EN'])


def reaction5_rate(parameters):
    return 8.5e-7 * (parameters["Te"] / 300.0) ** (-0.67)


reactions[4].rate_fun = reaction5_rate


def reaction6_rate(parameters):
    return 6.06e-6 / parameters['gas_temperature'] * np.exp(-15130.0 / parameters['gas_temperature'])


reactions[5].rate_fun = reaction6_rate


def reaction8_rate(parameters):
    return 8.75e-27 * (parameters["Te"] / 11600.0) ** (-4.5)


reactions[7].rate_fun = reaction8_rate


def reaction10_rate(parameters):
    return 2.25e-31 * (parameters['gas_temperature'] / 300.0) ** (-0.4)


reactions[9].rate_fun = reaction10_rate


def reaction_surface_rate(parameters):
    return 1.52 * (760 / parameters['gas_pressure']) * (parameters['gas_temperature'] / 273.16) * (
            parameters['Te'] / 11600) * \
           ((2.405 / parameters['radius']) ** 2 + (3.141 / parameters['gap_length']) ** 2)


reactions[10].rate_fun = reaction_surface_rate
reactions[11].rate_fun = reaction_surface_rate
reactions[12].rate_fun = reaction_surface_rate

mobility = parse_table_custom("mobility")  # in 1/m

#############################
# DYNAMIC
#############################
def update(parameters):
    neutral_particles = parameters["Ar"] + parameters["Ar*"]
    J = 1.6e-19 * parameters['gap_area'] * parameters['e'] * mobility(parameters['EN']) * 1 / 100 * parameters['EN'] * 1e-17
    parameters['EN'] = parameters['voltage'] / (parameters['gap_length'] + parameters['resistance'] * J
                                                / (parameters['EN'] * neutral_particles / 1.0e17 + 1.0e-99)) \
                       / neutral_particles * 1.0e17
    parameters['EN'] = 0.5 * (parameters['EN'] + abs(parameters['EN']))
    return parameters


# V  = voltage - resistance * J

times, values = solve(all_species, parameters, reactions, update=update)
plot(times, values, all_species)
