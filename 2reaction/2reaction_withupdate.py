from solver import solve_withN
from input_parser import parse_input_file
from plot import simple_plot


def update(parameters, time):
    if time < 1:
        parameters['EN'] = 50
    elif time < 2:
        parameters['EN'] = 20
    else:
        parameters['EN'] = 40


filename = "2reaction_2N_withupdate.input"
all_species, parameters, reactions, tables = parse_input_file(filename)


times, values = solve_withN(all_species, parameters, reactions, update=update)
simple_plot(times, values, all_species)