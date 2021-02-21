from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

gas_temperature = 301.0  # gas temperature, K
reduced_field = 50.0  # reduced electric field, Td
density_ini_ar = 2.5e19  # initial Ar density, cm-3
density_ini_elec = 1.0  # initial electron density, cm-3, also Ar^+ initial density

time_ini = 0.0
time_end = 3.0e-7
dtime = 1.0e-8  # times, s

k2 = 1.0e-25  # cm^6/s
k1 = 0.487e-11  # cm^3/s


def fun(t, y):
    e, Arp, Ar = y
    return [k1*e*Ar + -k2*e*Arp*Ar, k1*e*Ar + -k2*e*Arp*Ar, - k1*e*Ar + k2*e*Arp*Ar]


if __name__ == '__main__':
    sol = solve_ivp(fun, (time_ini, time_end), [density_ini_elec, density_ini_elec, density_ini_ar])
    plt.yscale("log")
    plt.plot(sol.t, sol.y[0])

