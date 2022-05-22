from solver import solve_numerical
from input_parser import parse_input_file, parse_table
from solver import solve_generic
import numpy as np

filename = "micro_cathode.input"
all_species, parameters, reactions, tables = parse_input_file(filename)

# -------------------------- PREPARE UPDATE METHOD ------------------------------------

# 1eV = 11600 K
E = parse_table("mean energy", tables)
Te = lambda prmtrs: E(prmtrs['EN']) * 11_600 * 2/3
mobility = parse_table("mobility", tables)  # in 1/m

q_elem = 1.60217662e-19
k_B = 1.38064852e-23
Td_to_Vm2 = 1e-21
Vm2_to_Td = 1 / Td_to_Vm2
gas_density = parameters['gas_pressure'] / (k_B * parameters['gas_temperature'])
gap_area = np.pi * parameters['radius'] ** 2
voltage = 1000.0
resistance = 1.0e5

def update(prmtrs, time):
    neutral_particles = prmtrs['Ar'] + prmtrs['Ar*']
    J = q_elem * gap_area * prmtrs['e'] * mobility(prmtrs['EN']) * prmtrs['EN'] * Td_to_Vm2
    prmtrs['EN'] = Vm2_to_Td * voltage / neutral_particles / \
                   (prmtrs['gap_length'] + resistance * J / (prmtrs['EN'] * Td_to_Vm2 * neutral_particles))
    prmtrs['EN'] = 0.5 * (prmtrs['EN'] + abs(prmtrs['EN']))

parameters['EN'] = voltage / parameters['gap_length'] / gas_density * Vm2_to_Td

# -------------------------- PREPARE REACTION METHODS ------------------------------------

# original values (with not SI units)
gas_pressure = 100  # torr
Tgas = 300  # K
radius = 0.4  # cm
gap_length = 0.4  # cm

def diff_rate(prmtrs):
    return 1240.946565968663 * E(prmtrs['EN'])

# convert from cm^3/s to m^3/s
reactant2conv = 1e-6
# convert from cm^6/s to m^6/s
reactant3conv = 1e-12

def reaction_Ar2withe(prmtrs):
    # Ar2^+ + e => Ar* + Ar          !   8.5d-7 * (Te/300.0d0)**(-0.67d0)
    basic = 8.5e-7 * (Te(prmtrs)/300.0)**(-0.67)
    return basic * reactant2conv

def reaction_ArpAre(prmtrs):
    # Ar^+ + Ar + e => Ar + Ar       !   8.75d-27 * (Te/11600.0d0)**(-4.5d0)
    basic = 8.75e-27 * (Te(prmtrs)/11600.0)**(-4.5)
    return basic * reactant3conv

# -------------------------- SET REACTION RATES ------------------------------------

# Ar2^+ + e => Ar* + Ar          !   8.5d-7 * (Te/300.0d0)**(-0.67d0)
reactions[4].rate_fun = reaction_Ar2withe

# Ar^+ + Ar + e => Ar + Ar       !   8.75d-27 * (Te/11600.0d0)**(-4.5d0)
reactions[7].rate_fun = reaction_ArpAre

# diff_rate: diffusion losses
reactions[10].rate_fun = diff_rate
reactions[11].rate_fun = diff_rate
reactions[12].rate_fun = diff_rate