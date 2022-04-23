import matplotlib.pyplot as plt
import random as rnd
import numpy as np

gas_temperature = 301.0  # gas temperature, K
reduced_field = 50.0  # reduced electric field, Td
density_ini_ar = 2.5e19  # initial Ar density, cm-3
density_ini_elec = 1.0  # initial electron density, cm-3, also Ar^+ initial density

time_ini = 0.0
time_end = 3.0e-7

k1 = 0.487e-11  # cm^3/s
k2 = 1.0e-25  # cm^6/s

# MICROMETERS RECALCULATION
density_ini_ar = 2.5e7  # micrometers^-3
density_ini_elec = 1     # micrometers^-3 (raised so that it works)

k1 = 0.487e1  # micrometers^3/s
k2 = 1.0e-1  # micrometers^6/s

electrons = 1
argon_plus = 1
argon = int(2.5e7)

plt.yscale("log")  # set the y-axis in the plot to logarithmic scale
time = time_ini
run = 0
while time < time_end:
    a1 = electrons*argon*k1
    a2 = electrons*argon*argon_plus*k2
    a0 = a1 + a2

    r1 = rnd.uniform(0, 1)
    tau = 1/a0*np.log(1/r1)

    r2 = rnd.uniform(0, 1)
    if r2 < a1/a0:
        electrons += 1
        argon -= 1
        argon_plus += 1
    else:
        electrons -= 1
        argon += 1
        argon_plus -= 1

    plt.plot(time, electrons, "b.")
    time += tau
    run += 1

plt.show()


