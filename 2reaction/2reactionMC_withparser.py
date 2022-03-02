from solver import solve
from input_parser import parse_input_file
from plot import plot

filename = "2reaction_2.input"
all_species, parameters, reactions, tables = parse_input_file(filename)

# resize attempt
factor = min([parameters[specie] for specie in all_species])
for specie in all_species:
    parameters[specie] /= factor
for reaction in reactions:
    reaction.ratio *= factor ** (len(reaction.reactants) - 1)

times, values = solve(all_species, parameters, reactions, bulk=1)

# undo resizing
for specie in values:
    values[specie] = list(map(lambda x: x * factor, values[specie]))

plot(times, values, all_species)
