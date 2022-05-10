from scipy.integrate import solve_ivp
import numpy as np

'''
def f(t, amp=120, base=20, x0=1e-7, w=2e-7, c=-1):  # used to compute EN in Td for t in s
    return (amp - base) / (1 + np.exp(-c * (t - x0) / w)) + base


def solve_numerical(all_species, parameters, reactions, tables):
    initial_concentrations = [parameters[specie] for specie in all_species]

    time_ini = parameters['time_ini']
    time_end = parameters['time_end']

    all_species = list(all_species)  # species need to be in an (any) order

    def fun(t, concentrations):
        differentials = []
        for specie in all_species:
            diff = 0
            for reaction in reactions:
                prod_change = reaction.products[specie] if specie in reaction.products else 0
                reac_change = reaction.reactants[specie] if specie in reaction.reactants else 0
                abs_change = prod_change - reac_change
                rate = reaction.rate_fun({'EN': f(t)})
                product = abs_change * rate
                for reactant in reaction.reactants:
                    reactant_index = all_species.index(reactant)
                    product *= concentrations[reactant_index] ** reaction.reactants[reactant]
                diff += product
            differentials.append(diff)

        return differentials

    sol = solve_ivp(fun, (time_ini, time_end), initial_concentrations)
    times = sol.t
    values = {all_species[i]: sol.y[i] for i in range(len(all_species))}

    return times, values
'''


def solve_numerical(all_species, parameters, reactions, update=None):
    initial_concentrations = [parameters[specie] for specie in all_species]

    time_ini = parameters['time_ini']
    time_end = parameters['time_end']

    all_species = list(all_species)  # species need to be in an (any) order

    def fun(t, concentrations):
        if update:
            update(parameters)

        differentials = []
        for specie in all_species:
            diff = 0
            for reaction in reactions:
                prod_change = reaction.products[specie] if specie in reaction.products else 0
                reac_change = reaction.reactants[specie] if specie in reaction.reactants else 0
                abs_change = prod_change - reac_change
                rate = reaction.rate_fun(parameters)
                product = abs_change * rate
                for reactant in reaction.reactants:
                    reactant_index = all_species.index(reactant)
                    product *= concentrations[reactant_index] ** reaction.reactants[reactant]
                diff += product
            differentials.append(diff)

        return differentials

    sol = solve_ivp(fun, (time_ini, time_end), initial_concentrations)
    times = sol.t
    values = {all_species[i]: sol.y[i] for i in range(len(all_species))}

    return times, values