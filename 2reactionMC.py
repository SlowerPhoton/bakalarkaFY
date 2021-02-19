from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import random as rnd
import numpy as np

gas_temperature = 301.0  # gas temperature, K
reduced_field = 50.0  # reduced electric field, Td
density_ini_ar = 2.5e19  # initial Ar density, cm-3
density_ini_elec = 1.0  # initial electron density, cm-3, also Ar^+ initial density

time_ini = 0.0
time_end = 3.0e-7

k2 = 1.0e-25  # cm^6/s
k1 = 0.487e-11  # cm^3/s

electrons = 1000
argon_plus = 1000
argon = int(2.5e19)

plt.yscale("log")  # set the y-axis in the plot to logarithmic scale
time = time_ini
while time < time_end:
    a1 = electrons*argon*k1
    a2 = electrons*argon*argon_plus*k2
    a0 = a1 + a2

    r1 = rnd.uniform(0, 1)
    tau = 1/a0*np.log(1/r1)

    r2 = rnd.uniform(0, 1)
    if r2 < a1/a0:
        # compute the probability from the probability density f = a1*exp(-a0*tau)
        Pr = a1/a0*(np.exp(-a0*time) - np.exp(-a0*(time+tau)))
        electrons += int(electrons*Pr)
        argon -= int(argon*Pr)
        argon_plus += int(argon_plus*Pr)
    else:
        Pr = a2/a0*(np.exp(-a0*time) - np.exp(-a0*(time+tau)))
        electrons -= int(electrons * Pr)
        argon += int(argon * Pr)
        argon_plus -= int(argon_plus * Pr)

    plt.plot(time, electrons, "b.")
    time += tau

plt.show()


