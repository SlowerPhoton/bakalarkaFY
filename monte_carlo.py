from solver import solve
from input_parser import parse_file

filename = "inputs/2reaction.input"
parameters, reactions = parse_file(filename)
solve(parameters, reactions)