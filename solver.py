import numpy as np
import random as rnd


def solve(all_species, parameters, reactions, update=None, bulk=1):
    times = []
    values = {species: [] for species in all_species}

    eps = 1e-10
    time = parameters["time_ini"]
    run = 0
    while time < parameters["time_end"]:
        if run % parameters['calc_step'] == 0:
            print(f"run: {run}, time: {time}, EN: {parameters['EN']}, Ar: {parameters['Ar']}, e: {parameters['e']}")
            for species in all_species:
                values[species].append(parameters[species])
            times.append(time)

        a = [reaction.compute_a(parameters) for reaction in reactions]  # TODO: should only update when necessary
        a0 = sum(a)
        if abs(a0) < eps:
            raise ZeroDivisionError("There is no possible reaction given the particle concentrations.")
        a_cum = [aa / a0 for aa in np.cumsum(a)]

        r1 = rnd.uniform(0, 1)
        tau = 1 / a0 * np.log(1 / r1) * bulk

        r2 = rnd.uniform(0, 1)
        for partial_sum in a_cum:
            if r2 < partial_sum:
                reaction_index = a_cum.index(partial_sum)
                break

        parameters = reactions[reaction_index].react(parameters, bulk)
        time += tau
        if update is not None:
            update(parameters)
        run += 1

    return times, values
