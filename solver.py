import matplotlib.pyplot as plt
import numpy as np
import random as rnd


def solve(parameters, reactions):
    plt.yscale("log")  # set the y-axis in the plot to logarithmic scale
    time = parameters["time_ini"]
    run = 0  # TODO: function needs to work with number of runs
    while time < parameters["time_end"]:
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

        plt.plot(time, parameters['e'], "b.")  # TODO: needs to be generalized
        time += tau
        run += 1
    plt.show()
