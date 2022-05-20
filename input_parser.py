"""
This file contains utilities for input files (and tables) parsing. Useful methods are:
 * parse_input_file
 * parse_table
"""

from reactions import Reaction
from scipy.interpolate import interp1d
import warnings


class InputFileError(Exception):
    """Base class for exceptions in this module."""
    pass


class InvalidLine(InputFileError):
    """Exception raised for a line that is neither a reaction nor a parameter set up.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


class MissingObligatoryParameter(InputFileError):
    """Exception raised when an obligatory parameter is missing.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


class TableError(Exception):
    """Exception raised when there is a table file syntax error.

        Attributes:
            message -- explanation of the error
        """


def parse_input_file(filename):
    """
    Parses an input file (typically with .input extension).

    filename ... path to the input file

    Returns tuple (all_species, parameters, reactions, tables).
    all_species ... set of all species participating in the reactions, e.g., {'Ar', 'e', 'Ar*'}
    parameters ... dictionary storing parameters from the input file
    reactions ... list of Reaction objects (i.e. reaction[0] corresponds to the first reaction, stored
        as instance of the class Reaction specified in reactions.py)
    tables ... dictionary of tables indexed by their names in the table file, each table is then a list of rows
        empty if no 'table_file' specified
    """
    parameter_lines = []
    reaction_lines = []
    with open(filename) as file:
        for line in file.readlines():
            _parse_line(line, parameter_lines, reaction_lines)

    parameters = {}
    for parameter_line in parameter_lines:
        _parse_parameter(parameter_line, parameters)
    _check_obligatory_parameters(parameters)  # check for any missing obligatory parameter

    _set_default_parameters(parameters)  # set other missing parameters to their default values

    tables = _parse_tables(parameters)

    all_species = set()
    reactions = []
    for reaction_line in reaction_lines:
        _parse_reaction(reaction_line, all_species, parameters, reactions, tables)

    _set_missing_species(parameters, all_species)  # unspecified initial concentrations set to 0
    return all_species, parameters, reactions, tables


def _parse_line(line, parameter_lines, reaction_lines):
    line = line.strip()  # get rid of whitespace L and R
    if line == "":  # skip empty lines
        return
    if line[0] == "#":  # skip comments
        return
    if '=>' in line:
        reaction_lines.append(line)
    elif '=' in line:
        parameter_lines.append(line)
    else:
        raise InvalidLine(f"There is a syntax error on line '{line}'. Each line must be empty, start "
                          f"with the hashtag ('#'), set up a parameter or specify a reaction.")


def _parse_parameter(line, parameters):
    line_split = line.split('=')
    if len(line_split) != 2:
        raise InvalidLine(f"Syntax error on line '{line}': Parameters are set using the '=' equals sign. There can "
                          f"only be one parameter set per one line.")
    parameter = line_split[0].strip()
    try:
        value = float(line_split[1]) if ("." in line_split[1] or "e" in line_split[1]) else int(line_split[1])
    except ValueError:
        value = line_split[1].strip()
    parameters[parameter] = value


def _parse_reaction(line, all_species, parameters, reactions, tables):
    line_split_arrow = line.split('=>')
    if len(line_split_arrow) != 2:
        raise InvalidLine(f"Syntax error on line '{line}': There can only be one reaction per line. In a reaction, "
                          f"the reactants and products are separated by the arrow '=>'.")
    left_side = line_split_arrow[0]
    line_split_excl = line_split_arrow[1].split('!')
    if len(line_split_excl) != 2:
        raise InvalidLine(f"Syntax error on line '{line}': Each reaction must have its rate specified. The rate "
                          f"specification goes after the single exclamation mark '!'. The exclamation mark '!' "
                          f"must not be used elsewhere.")
    right_side = line_split_excl[0]
    rate_spec = line_split_excl[1]

    reaction = Reaction(_parse_species(left_side, all_species), _parse_species(right_side, all_species))
    _parse_rate(reaction, rate_spec, tables, parameters)  # fills attributes rate_fun and table_name
    reactions.append(reaction)


def _parse_species(line, all_species):
    species = line.split()
    # + can be used as a delimiter but it is unnecessary
    # species must be separated by whitespace
    # '+' delimited by whitespace are ignored, otherwise they count as part of the species name
    while '+' in species:
        species.remove('+')
    all_species.update(species)
    return species


def _parse_rate(reaction, rate_spec, tables, parameters):
    fun = None
    try:
        num = float(rate_spec)
        fun = lambda prmtrs: num
    except ValueError:
        rate_spec = rate_spec.strip()
        if rate_spec.startswith("table:"):
            table_name = rate_spec[len("table:"):].strip()
            reaction.table_name = table_name
            fun = _parse_table(table_name, tables, parameters)
        else:
            warnings.warn(f"The rate specification '{rate_spec}' cannot be parsed, because we only support "
                          f"constant rates and rates given by a table. You must write the rate function in "
                          f"the code yourself and assign it to the corresponding reaction.", UserWarning)
    reaction.rate_fun = fun


def _parse_tables(parameters):
    tables = {}
    if "table_file" not in parameters:
        return {}
    with open(parameters["table_file"]) as file:
        expecting_new_table = True
        for line in file.readlines():
            line = line.strip()
            if line == "":  # skip empty lines
                continue
            if line[0] == "#":  # skip comments (lines starting with '#')
                continue
            if expecting_new_table:
                new_table_name = line
                table_content = []
                expecting_new_table = False
                ignore_lines = True
            elif ignore_lines and line == '-' * len(line):  # table content separated by '----------'
                ignore_lines = False
            elif ignore_lines:  # ignore lines between table name and table content
                continue
            elif line == '-' * len(line):  # end of table (specified by '---------')
                tables[new_table_name] = table_content
                expecting_new_table = True
            else:
                table_content.append(line)
    return tables


def _parse_table(table_name, tables, parameters):
    if table_name not in tables:
        err_str = f"There is no table {table_name}."
        if "table_name" not in parameters:
            err_str += " Maybe you forgot to specify the 'table_name' parameter?"
        raise TableError(err_str)

    content = tables[table_name]
    first_col = []
    second_col = []
    for line in content:
        cols = line.split()
        if len(cols) != 2:
            raise TableError(f"Syntax error for table '{table_name}' on table row '{line}': "
                             f"Each table row must have exactly two columns.")
        first_col.append(float(cols[0]))
        second_col.append(float(cols[1]))
    f = interp1d(first_col, second_col)  # the default is linear interpolation
    return lambda prmtrs: f(prmtrs["EN"])


def parse_table(table_name, tables):
    """
    Can be used to parse a table from the variable tables produced by the 'parse_input_file' method.

    It returns a function as a linear interpolation of the table. The first column is the input, the second column is
    the output value of the function.
    """
    content = tables[table_name]
    first_col = []
    second_col = []
    for line in content:
        cols = line.split()
        if len(cols) != 2:
            raise TableError(f"Syntax error for table '{table_name}' on table row '{line}': "
                             f"Each table row must have exactly two columns.")
        first_col.append(float(cols[0]))
        second_col.append(float(cols[1]))
    f = interp1d(first_col, second_col)  # the default is linear interpolation
    return f


def _check_obligatory_parameters(parameters):
    obligatory_parameters = ['time_end']
    for obligatory_parameter in obligatory_parameters:
        if obligatory_parameter not in parameters:
            raise MissingObligatoryParameter(f"Parameter '{obligatory_parameter}' is not specified.")


def _set_default_parameters(parameters):
    if 'time_ini' not in parameters:
        parameters['time_ini'] = 0
        warnings.warn(f"Parameter 'time_ini' is missing: it is set to 0.", UserWarning)
    if 'calc_step' not in parameters:
        parameters['calc_step'] = 1
        warnings.warn(f"Parameter 'calc_step' is missing: it is set to 1.", UserWarning)


def _set_missing_species(parameters, all_species):
    for species in all_species:
        if species not in parameters:
            parameters[species] = 0
            warnings.warn(f"Parameter '{species}' is missing: it is set to 0.", UserWarning)
