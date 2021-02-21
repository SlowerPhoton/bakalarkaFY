from solver import solve
from input_parser import parse_file
from plot import plot

filename = "inputs/2reaction.input"
all_species, parameters, reactions = parse_file(filename)
times, values = solve(all_species, parameters, reactions)
plot(times, values, all_species)
