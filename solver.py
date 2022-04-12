import numpy as np
import random as rnd


def bulk_from_N(all_species, parameters):
    """
    if the parameter N was specified in the input file, compute the (super)particle concentrations and compute
    their weights
    """
    if 'N' not in parameters:
        return None
    weight = min([parameters[specie] for specie in all_species]) / parameters['N']
    return weight


def solve(all_species, parameters, reactions, update=None, bulk=1):
    """
    all_species ... is a set of specie names: e.g. {'Ar^+', 'e', 'Ar'}
    parameters ... is a dictionary containing info parsed from the input file, such as time_ini, calc_step, N,
        initial concentrations...
    reactions ... is a list of Reaction objects (i.e. reaction[0] corresponds to the first reaction, stored
        as instance of the class Reaction specified in reactions.py)
    bulk ... (now obsolete and possibly will be removed)
        specifies how many particles are processed at once in each iteration (i.e. bulk=2 means the selected
        reaction is executed twice in each iteration, basically corresponds to particle weight, but the name differs
        because bulk needs to be an integer
    """

    # update concentrations and compute bulk if N specified in parameters
    weight = bulk_from_N(all_species, parameters)
    if weight: bulk = weight  # if N not specified, do nothing

    times = []  # this will store the timestamps
    values = {species: [] for species in all_species}  # stores the concentrations per specie for each timestamp

    eps = 1e-10  # for checking a0 is not too small -> no possible reaction
    time = parameters["time_ini"]
    run = 0
    while time < parameters["time_end"]:
        # only save next timestamp (and concentrations) after calc_step steps
        if run % parameters['calc_step'] == 0:
            # print(f"run: {run}, time: {time}, EN: {parameters['EN']}, Ar: {parameters['Ar']}, e: {parameters['e']}")
            for species in all_species:
                values[species].append(parameters[species])
            times.append(time)

        # sample a reaction and let it react
        a = [reaction.compute_a(parameters) for reaction in reactions]  # TODO: should only update when necessary
        #print(f"a={a}")
        a0 = sum(a)
        if abs(a0) < eps:
            raise ZeroDivisionError("There is no possible reaction given the particle concentrations.")
        a_cum = [aa / a0 for aa in np.cumsum(a)]

        r2 = rnd.uniform(0, 1)
        for partial_sum in a_cum:
            if r2 < partial_sum:
                reaction_index = a_cum.index(partial_sum)
                break
        parameters = reactions[reaction_index].react(parameters, bulk)

        # sample a time delta
        r1 = rnd.uniform(0, 1)
        tau = 1 / a0 * np.log(1 / r1) * weight
        #print("tau=", tau, " for bulk=", bulk, " tau/bulk=", tau/bulk)
        time += tau

        # run the update function on the parameters to modify them
        if update is not None:
            update(parameters)
        run += 1

    return times, values
