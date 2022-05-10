from MicroCathode.numerical_simple import solve_numerical
from input_parser import parse_input_file, parse_table
from solver import solve, quick_solve
from plot import plot
import numpy as np

filename = "micro_cathode_simple.input"
all_species, parameters, reactions, tables = parse_input_file(filename)

# -------------------------- PREPARE UPDATE METHOD ------------------------------------

# 1eV = 11600 K
E = parse_table("mean energy", tables)
Te = lambda prmtrs: E(prmtrs['EN']) * 11_600

def f(t, amp=120, base=20, x0=1e-7, w=2e-7, c=-1):  # used to compute EN in Td for t in s
    return (amp - base) / (1 + np.exp(-c * (t - x0) / w)) + base

def update(prmtrs, time):
    prmtrs['EN'] = f(time)

# set EN at the beginning of th esimulation
update(parameters, 0)

# -------------------------- PREPARE REACTION METHODS ------------------------------------

# original values (with not SI units)
gas_pressure = 100  # torr
Tgas = 300  # K
radius = 0.4  # cm
gap_length = 0.4  # cm

def diff_rate(prmtrs):
    assert_tolerance = 1e-9
    assert abs(1.52 * (760 / gas_pressure) * (Tgas / 273.16) * (Te(prmtrs) / 11600) * ((2.405/radius)**2 + (3.141/gap_length)**2) \
           - 1240.946565968663 * E(prmtrs['EN'])) < assert_tolerance
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

def reaction_ArsToArp(prmtrs):
    # e + Ar * = > Ar ^ + + e + e      !   1.5 * ionization
    ionization_rate = parse_table("ionization", tables)
    return 1.5 * ionization_rate(prmtrs['EN'])

# -------------------------- SET REACTION RATES ------------------------------------
# e + Ar* => Ar^+ + e + e        !   1.5 * ionization
reactions[3].rate_fun = reaction_ArsToArp

# Ar2^+ + e => Ar* + Ar          !   8.5d-7 * (Te/300.0d0)**(-0.67d0)
reactions[4].rate_fun = reaction_Ar2withe

# Ar^+ + Ar + e => Ar + Ar       !   8.75d-27 * (Te/11600.0d0)**(-4.5d0)
reactions[7].rate_fun = reaction_ArpAre

# diff_rate: diffusion losses
reactions[10].rate_fun = diff_rate
reactions[11].rate_fun = diff_rate
reactions[12].rate_fun = diff_rate


# -------------------------- SOLVE AND PLOT ------------------------------------
times, values = solve(all_species, parameters, reactions, update=update, recompute_N=False)
'''
times, values = solve_numerical(all_species, parameters, reactions, tables)
n_drop = 5
values_drop = {specie: values[specie][n_drop:] for specie in all_species}
values = values_drop
all_species.remove('Ar2(W)^+')
all_species.remove('Ar2^+')
all_species.remove('Ar*')
times = times[n_drop:]
'''
plot(times, values, all_species)