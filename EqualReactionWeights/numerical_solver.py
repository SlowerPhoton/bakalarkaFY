from scipy.integrate import solve_ivp


def solve_numerical(all_species, parameters, reactions, tables):
    density_ini_ar = parameters['Ar']
    density_ini_elec = parameters['e']
    density_ini_arp = parameters["Ar^+"]
    density_ini_eW = parameters['e(W)']

    time_ini = parameters['time_ini']
    time_end = parameters['time_end']

    # since the reaction is constant, it doesn't care about the input parameters
    k1 = reactions[0].rate_fun(None)
    k2 = reactions[1].rate_fun(None)
    k3 = reactions[2].rate_fun(None)

    def fun(t, y):
        e, Arp, Ar, eW = y
        return [k1 * e * Ar + -k2 * e * Arp * Ar -k3 * e,
                k1 * e * Ar + -k2 * e * Arp * Ar,
                - k1 * e * Ar + k2 * e * Arp * Ar,
                k3 * e]

    sol = solve_ivp(fun, (time_ini, time_end), [density_ini_elec, density_ini_arp, density_ini_ar, density_ini_eW])
    times = sol.t
    values = {'e': sol.y[0], 'Ar^+': sol.y[1], 'Ar': sol.y[2], 'e(W)': sol.y[3]}

    return times, values