from scipy.integrate import solve_ivp


def solve_numerical(all_species, parameters, reactions, tables):
    density_ini_ar = parameters['Ar']
    density_ini_elec = parameters['e']
    density_ini_arp = parameters["Ar^+"]

    time_ini = parameters['time_ini']
    time_end = parameters['time_end']

    k1 = reactions[0].rate_fun("placeholder")
    k2 = reactions[1].rate_fun("placeholder")

    def fun(t, y):
        e, Arp, Ar = y
        return [k1 * e * Ar + -k2 * e * Arp * Ar, k1 * e * Ar + -k2 * e * Arp * Ar, - k1 * e * Ar + k2 * e * Arp * Ar]

    sol = solve_ivp(fun, (time_ini, time_end), [density_ini_elec, density_ini_arp, density_ini_ar])
    times = sol.t
    values = {'e': sol.y[0], 'Ar^+': sol.y[1], 'Ar': sol.y[2]}

    return times, values

'''
if __name__ == '__main__':
    
    plt.yscale("log")
    plt.plot(sol.t, sol.y[0], "--", label=r"$\left[e^{-}\right]_t$")
    plt.plot(sol.t, sol.y[1], ".", label=r"$\left[Ar^{+}\right]_t$")
    plt.plot(sol.t, sol.y[2], "--", label=r"$\left[Ar\right]_t$")
    plt.legend()
    plt.xlabel("time $t$ [s]")
    plt.ylabel(r"concentration $\left[ \rm cm^{-3} \right]$")
    plt.show()
'''