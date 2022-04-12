from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from input_parser import parse_input_file

filename = "2reaction_2N.input"
all_species, parameters, reactions, tables = parse_input_file(filename)

density_ini_ar = parameters['Ar']
density_ini_elec = parameters['e']
density_ini_arp = parameters["Ar^+"]


time_ini = parameters['time_ini']
time_end = parameters['time_end']

k1 = reactions[0].rate_fun("placeholder")
k2 = reactions[1].rate_fun("placeholder")


def fun(t, y):
    e, Arp, Ar = y
    return [k1*e*Ar + -k2*e*Arp*Ar, k1*e*Ar + -k2*e*Arp*Ar, - k1*e*Ar + k2*e*Arp*Ar]


if __name__ == '__main__':
    sol = solve_ivp(fun, (time_ini, time_end), [density_ini_elec, density_ini_arp, density_ini_ar])
    plt.yscale("log")
    plt.plot(sol.t, sol.y[0], "--", label=r"$\left[e^{-}\right]_t$")
    plt.plot(sol.t, sol.y[1], ".", label=r"$\left[Ar^{+}\right]_t$")
    plt.plot(sol.t, sol.y[2], "--", label=r"$\left[Ar\right]_t$")
    plt.legend()
    plt.xlabel("time $t$ [s]")
    plt.ylabel(r"concentration $\left[ \rm cm^{-3} \right]$")
    plt.show()
