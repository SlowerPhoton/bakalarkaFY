import matplotlib.pyplot as plt
import numpy as np
import random as rnd


def solve(all_species, parameters, reactions, update=None):
    times = []
    values = {species : [] for species in all_species}

    time = parameters["time_ini"]
    run = 0
    while time < parameters["time_end"]:
        if run % parameters['calc_step'] == 0:
            print(run, time)
            for species in all_species:
                values[species].append(parameters[species])
            times.append(time)

        a = [reaction.compute_a(parameters) for reaction in reactions]  # TODO: should only update when necessary
        a0 = sum(a)
        a_cum = [aa / a0 for aa in np.cumsum(a)]

        r1 = rnd.uniform(0, 1)
        tau = 1 / a0 * np.log(1 / r1)

        r2 = rnd.uniform(0, 1)
        for partial_sum in a_cum:
            if r2 < partial_sum:
                reaction_index = a_cum.index(partial_sum)
                break

        parameters = reactions[reaction_index].react(parameters)
        time += tau
        if update is not None:
            parameters = update(parameters)
        run += 1

    return times, values
