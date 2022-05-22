"""
This file contains methods for both stochastic and deterministic simulations.
"""
import numpy as np
import random as rnd

from scipy.integrate import solve_ivp


def solve_generic(selected_params, parameters, reactions, update=None, bulk=1, bulk_compute=None,
                  print_out=None, outfile=None, ERW=False):
    """
    General method for Monte Carlo simulations. Best used for rather simple systems, otherwise implementing
    your own domain-specific function is superior in efficiency.

    selected_params ... subset of parameters keys (typically equal to all_species returned by input_parser)
    parameters ... is a dictionary containing info parsed from the input file, such as time_ini, calc_step, N,
        species concentrations...
    reactions ... is a list of Reaction objects (i.e. reaction[0] corresponds to the first reaction, stored
        as instance of the class Reaction specified in reactions.py)
    update ... method called at the end of each iteration to update parameters: update(parameters, time)
    bulk ... specifies how many reactions are processed at once in each iteration (i.e. bulk=2 means the selected
        reaction is executed twice in each iteration)
    bulk_compute ... if specified, it is called periodically after calc_step iterations to update bulk:
                    bulk = bulk_compute(run, time, parameters, bulk)
    print_out ... if specified, called periodically after calc_step iterations to print/log progress,
                    with signature: print_out(run, time, parameters)
    outfile ... [optional] output filename
    ERW ... if True, uses equal reaction weights for the simulation

    The function returns tuple times, values if outfile is None, otherwise it returns None and saves the output
    continually in the output file. The contents of the output file can then be read and parsed into the times, values
    tuple using the method read_outfile from plot.py.

    times ... list of timestamps (time for each calc_step-th iteration of the algorithm)
    values ... dictionary of values of selected_params from parameters for each calc_step-th iteration of the algorithm
    E.g.; times = [0, 0.5, 1], values = {'e': [100, 98, 96], 'Ar': [1000, 1000, 1000]}
    """

    if outfile:
        out = open(outfile, "w")
    else:
        times = []  # this will store the timestamps
        values = {param: [] for param in selected_params}  # stores the parameters values for each timestamp

    eps = 1e-10  # for checking a0 is not too small -> no possible reaction
    time = parameters["time_ini"]
    run = 0
    while time < parameters["time_end"]:

        # after calc_step iterations
        if run % parameters['calc_step'] == 0:
            if print_out: print_out(run, time, parameters)  # print out computation progress
            if bulk_compute:  # update bulk value
                bulk = bulk_compute(run, time, parameters, bulk)
            if not outfile:  # save current time & parameters including concentrations
                for param in selected_params:
                    values[param].append(parameters[param])
                times.append(time)
            else:  # write current time & selected parameters in the output file
                out.write(f"time: {time}")
                for param in selected_params:
                    out.write(f", {param}: {parameters[param]}")
                out.write("\n")
                out.flush()

        # sample a reaction and let it react
        a = [reaction.compute_a(parameters) for reaction in reactions]
        a0 = sum(a)
        if abs(a0) < eps:
            raise ZeroDivisionError("There is no possible reaction given the particle concentrations.")
        a_cum = [aa / a0 for aa in np.cumsum(a)]

        # choose the reaction
        if ERW:
            chosen_reaction_index = rnd.randrange(len(reactions))
            # the weight of the chosen reaction is given by bulk and its transition rate
            weight = bulk * len(reactions) * a[chosen_reaction_index] / a0
        else:
            r2 = rnd.uniform(0, 1)
            for partial_sum in a_cum:
                if r2 < partial_sum:
                    chosen_reaction_index = a_cum.index(partial_sum)
                    break
            weight = bulk  # the weight of the chosen reaction is given by bulk
        parameters = reactions[chosen_reaction_index].react(parameters, weight)

        # sample a time delta
        r1 = rnd.uniform(0, 1)
        tau = 1 / a0 * np.log(1 / r1) * weight
        time += tau

        if update: update(parameters, time=time)  # run the update function on the parameters to modify them
        run += 1

    if outfile:  # close the file and return None
        out.close()
    else:  # return the computed concentrations
        return times, values


def solve_withN(all_species, parameters, reactions, update=None, bulk=1, recompute_N=True, main_specie='e',
                verbose=False):
    """
    A more specific derivative of the method 'solve_generic' working explicitly with the number of superparticles N.
    The number of superparticles is considered for the main_specie.
    N is computed as parameters[main_specie] / parameters['N'].

    all_species ... is a set of species names: e.g. {'Ar^+', 'e', 'Ar'}
    parameters ... is a dictionary containing info parsed from the input file, such as time_ini, calc_step, N,
        species concentrations...
    reactions ... is a list of Reaction objects (i.e. reaction[0] corresponds to the first reaction, stored
        as instance of the class Reaction specified in reactions.py)
    update ... method called at the end of each iteration to update parameters: update(parameters, time)
    bulk ... specifies how many reactions are processed at once in each iteration (i.e. bulk=2 means the selected
        reaction is executed twice in each iteration)
    recompute_N ... if specified, when actual number of superparticles of main_specie is bigger than 2N or lower
        than 0.5N, the superparticle weight is recomputed so that the actual number of superparticles is N again
    main_specie ... a species name from all_species, N then represents superparticles of this species
    verbose ... if True, prints progress periodically after calc_step iterations

    Returns tuple (times, values).
    times ... list of timestamps (time for each calc_step-th iteration of the algorithm)
    values ... dictionary of values of selected_params from parameters for each calc_step-th iteration of the algorithm
    E.g.; times = [0, 0.5, 1], values = {'e': [100, 98, 96], 'Ar': [1000, 1000, 1000]}
    """

    # compute bulk if N specified in parameters
    if 'N' in parameters:  # if N not specified, do nothing
        bulk = parameters[main_specie] / parameters['N']

    times = []  # this will store the timestamps
    values = {species: [] for species in all_species}  # stores the concentrations per specie for each timestamp

    eps = 1e-10  # for checking a0 is not too small -> no possible reaction
    time = parameters["time_ini"]
    run = 0
    while time < parameters["time_end"]:

        # only save next timestamp (and concentrations) after calc_step steps
        if run % parameters['calc_step'] == 0:
            if verbose:
                print(f"run: {run}, time: {time}, EN: {parameters['EN']}, Ar: {parameters['Ar']}, e: {parameters['e']}")
            for species in all_species:
                values[species].append(parameters[species])
            times.append(time)

        # sample a reaction and let it react
        a = [reaction.compute_a(parameters) for reaction in reactions]
        a0 = sum(a)
        if abs(a0) < eps:
            raise ZeroDivisionError("There is no possible reaction given the particle concentrations.")
        a_cum = [aa / a0 for aa in np.cumsum(a)]

        # choose next reaction
        r2 = rnd.uniform(0, 1)
        for partial_sum in a_cum:
            if r2 < partial_sum:
                reaction_index = a_cum.index(partial_sum)
                break
        parameters = reactions[reaction_index].react(parameters, bulk)

        # sample a time delta
        r1 = rnd.uniform(0, 1)
        tau = 1 / a0 * np.log(1 / r1) * bulk
        time += tau

        # run the update function on the parameters to modify them
        if update: update(parameters, time=time - tau)
        run += 1

        # check if particles need rescaling (N and bulk recomputation)
        if recompute_N:
            actual_N = parameters[main_specie] / bulk
            # if the current number of superparticles differs too much from the original N
            if actual_N > parameters['N'] * 2 or actual_N < parameters['N'] * 0.5:
                bulk = parameters[main_specie] / parameters['N']  # rescale to create N superparticles again

    return times, values


def solve_numerical(all_species, parameters, reactions, update=None, method=None):
    """
    Generic numerical deterministic solver method

    For more complicated systems, this method might be too slow as it doesn't use any domain knowledge
    about the modeled system (it treats the method update as a blackbox).

    all_species ... is a set of specie names: e.g., {'Ar^+', 'e', 'Ar'}
    parameters ... is a dictionary containing info parsed from the input file, such as time_ini, calc_step,
        initial concentrations
    reactions ... is a list of Reaction objects (i.e. reaction[0] corresponds to the first reaction, stored
        as instance of the class Reaction specified in reactions.py), the solver will access it to retrieve
        reaction rates
    update [optional] ... method update(parameters, time) updates any parameters based on the simulation time,
        the parameters can be used as an input for reaction rates
    method [optional] ... if specified, it is passed to solve_ivp, for possible options see docs:
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
    """
    all_species = list(all_species)  # species need to be in an (any) order (deals with Set)

    initial_concentrations = [parameters[specie] for specie in all_species]

    time_ini = parameters['time_ini']
    time_end = parameters['time_end']

    # this function will be integrated (takes concentrations and time and returns concentration diffs)
    def fun(t, concentrations):
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

        if update:
            update(parameters, time=t)

        return differentials

    if method:
        sol = solve_ivp(fun, (time_ini, time_end), initial_concentrations, method=method)
    else:
        sol = solve_ivp(fun, (time_ini, time_end), initial_concentrations)
    times = sol.t
    values = {all_species[i]: sol.y[i] for i in range(len(all_species))}

    return times, values


def solve_numerical_EN(all_species, parameters, reactions, update, precision=1e4, verbose=False):
    """
    Numerical deterministic solver method made specifically for systems depending solely on E/N ratio

    It splits the time interval time_ini - time_end logarithmically into precision steps and solves each of
    them with constant E/N.

    all_species ... is a set of specie names: e.g., {'Ar^+', 'e', 'Ar'}
    parameters ... is a dictionary containing info parsed from the input file, such as time_ini, calc_step,
        initial concentrations
    reactions ... is a list of Reaction objects (i.e. reaction[0] corresponds to the first reaction, stored
        as instance of the class Reaction specified in reactions.py), the solver will access it to retrieve
        reaction rates
    update ... method update(parameters, time) updates any parameters based on the simulation time,
        the parameters can be used as an input for reaction rates, it should update EN value
    precision ... number of time steps approximated as having constant E/N ratio
    verbose ... if True, prints out computation progress
    """
    all_species = list(all_species)  # species need to be in an (any) order

    initial_concentrations = [parameters[specie] for specie in all_species]

    time_ini = parameters['time_ini']
    time_end = parameters['time_end']

    times = [time_ini]
    values = {species: [parameters[species]] for species in all_species}
    values['EN'] = [parameters['EN']]

    timestamps = np.logspace(np.log10(time_ini), np.log10(time_end), num=int(precision))

    # this function will be integrated (takes concentrations and time and returns concentration diffs)
    # it doesn't use the parameter t, as it assumes constant E/N set externally
    def fun(t, concentrations):
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

    # loop through timestamps
    for i in range(int(precision) - 1):
        update(parameters, time=timestamps[i])  # set the current E/N as constant for the timestamp
        sol = solve_ivp(fun, (timestamps[i], timestamps[i + 1]), initial_concentrations, method='Radau')

        # store results
        times += list(sol.t)
        for j in range(len(all_species)):
            values[all_species[j]] += list(sol.y[j])

        # prepare initial_concentrations for the next computation
        initial_concentrations = [values[species][-1] for species in all_species]

        if verbose and i % (int(precision) // 100) == 0:  # print out progress for each percentage
            print(f"{int(i / int(precision)  * 100)} %, run {i}/{precision}, time: {timestamps[i]}")
        values['EN'] += [parameters['EN'] for _ in sol.t]  # store E/N ratio used in this timestamp

    return times, values
