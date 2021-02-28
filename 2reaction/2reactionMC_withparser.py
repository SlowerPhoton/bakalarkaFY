from solver import solve
from input_parser import parse_input_file
from plot import plot

filename = "2reaction.input"
all_species, parameters, reactions = parse_input_file(filename)
times, values = solve(all_species, parameters, reactions)
plot(times, values, all_species)
