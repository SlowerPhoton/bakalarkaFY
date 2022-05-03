from input_parser import parse_input_file, parse_table
import numpy as np

filename = "micro_cathode_simple.input"
all_species, parameters, reactions, tables = parse_input_file(filename)

# 1eV = 11600 K
E = parse_table("mean energy", tables)
Te = lambda prmtrs: E(prmtrs['EN']) * 11_600

def f(t, amp=120, base=20, x0=1e-7, w=2e-7, c=-1):  # used to compute EN in Td for t in s
    return (amp - base) / (1 + np.exp(-c * (t - x0) / w)) + base

def diff_rate(prmtrs):
    # 1.52 * (760 / gas_pressure) * ( Tgas / 273.16 ) * ( Te / 11600 ) * ( (2.405/radius)**2 + (3.141/gap_length)**2 )
    return 1240.946565968663 * E(prmtrs['EN'])


