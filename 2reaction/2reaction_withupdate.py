from solver import solve
from input_parser import parse_input_file
from plot import plot


def update(parameters, time):
    if time < 1:
        parameters['EN'] = 50
    elif time < 2:
        parameters['EN'] = 20
    else:
        parameters['EN'] = 40


filename = "2reaction_2N_withupdate.input"
all_species, parameters, reactions, tables = parse_input_file(filename)


times, values = solve(all_species, parameters, reactions, update=update)
plot(times, values, all_species)