from solver import solve_numerical
from input_parser import parse_input_file, parse_table
from solver import solve_withN, solve_generic
from plot import simple_plot
import numpy as np

filename = "micro_cathode.input"
all_species, parameters, reactions, tables = parse_input_file(filename)

# -------------------------- PREPARE UPDATE METHOD ------------------------------------

# 1eV = 11600 K
E = parse_table("mean energy", tables)
Te = lambda prmtrs: E(prmtrs['EN']) * 11_600 * 2 / 3


def f(t, amp=120, base=20, x0=1e-7, w=2e-7, c=-1):  # used to compute EN in Td for t in s
    return (amp - base) / (1 + np.exp(-c * (t - x0) / w)) + base


def update(prmtrs, time):
    prmtrs['EN'] = f(time, 75.0, 3.0, 7e-8, 5e-8, -1.0) + f(time, 35.0, 0.0, 1e-6, 1e-6, -1.0)


# set EN at the beginning of the simulation
update(parameters, parameters['time_ini'])

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
    basic = 8.5e-7 * (Te(prmtrs) / 300.0) ** (-0.67)
    return basic * reactant2conv


def reaction_ArpAre(prmtrs):
    # Ar^+ + Ar + e => Ar + Ar       !   8.75d-27 * (Te/11600.0d0)**(-4.5d0)
    basic = 8.75e-27 * (Te(prmtrs) / 11600.0) ** (-4.5)
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


# -------------------------- SOLVE ------------------------------------
def exp_bulk_compute(bulk, run, time, parameters):
    return bulk * 1.01


def print_out(run, time, parameters):
    print(f"run: {run}, time: {time}, EN: {parameters['EN']}, e: {parameters['e']}, Ar: {parameters['Ar^+']}")


solve_generic(all_species, parameters, reactions, update=update, bulk_compute=exp_bulk_compute,
              print_out=print_out, outfile='micro_cathode_improved.txt')
