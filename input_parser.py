from reactions import Reaction
from scipy.interpolate import interp1d

parameters = {}
reactions = []
all_species = set()


class InputFileError(Exception):
    """Base class for exceptions in this module."""
    pass


class InputFileWarning(Warning):
    """Base class for warnings in this module."""
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


class MissingParameter(InputFileWarning):
    """Warning raised when a parameter is missing and set to a default value.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


class UnsupportedFeature(InputFileWarning):
    """Warning raised for a feature that is not yet supported but it is called in the input file.

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


def parse_file(filename):
    with open(filename) as file:
        for line in file.readlines():
            line = line.strip()  # get rid of whitespace L and R
            if line == "":  # skip empty lines
                continue
            if line[0] == "#":  # skip comments
                continue
            parse_line(line)
    check_parameters()  # check for any missing obligatory parameter
    set_default_parameters()  # set other missing parameters to their default values
    return all_species, parameters, reactions


def parse_line(line):
    if '=>' in line:
        parse_reaction(line)
    elif '=' in line:
        parse_parameter(line)
    else:
        raise InvalidLine(f"There is a syntax error on line '{line}'. Each line must be empty, start "
                          f"with the hashtag ('#'), set up a parameter or specify a reaction.")


def parse_parameter(line):
    line_split = line.split('=')
    if len(line_split) != 2:
        raise InvalidLine(f"Syntax error on line '{line}': Parameters are set using the '=' equals sign. There can "
                          f"only be one parameter set per one line.")
    parameter = line_split[0].strip()
    value = float(line_split[1]) if "." in line_split[1] else int(line_split[1])
    parameters[parameter] = value


def parse_reaction(line):
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

    reaction = Reaction(parse_species(left_side), parse_species(right_side), parse_rate(rate_spec))
    reaction.rate_spec = rate_spec  # remember rate specification for debugging purposes
    reactions.append(reaction)


def parse_species(line):
    species = line.split()
    # + can be used as a delimiter but it is unnecessary
    # species must be separated by whitespace
    # '+' delimited by whitespace are ignored, otherwise they count as part of the species name
    while '+' in species:
        species.remove('+')
    all_species.update(species)
    return species


def parse_rate(rate_spec):
    fun = None
    try:
        num = float(rate_spec)
        fun = lambda energy: num
    except ValueError:
        rate_spec = rate_spec.strip()
        if rate_spec.startswith("table:"):
            table_name = rate_spec[len("table:"):].strip()
            fun = parse_table(table_name)
        else:
            raise UnsupportedFeature(f"The rate specification '{rate_spec}' cannot be parsed, because we only support "
                                     f"constant rates and rates given by a table. You must write the rate function in "
                                     f"the code yourself and assign it to the corresponding reaction.")
    return fun


tables = {}


def parse_tables():
    if "table_file" not in parameters:
        raise MissingObligatoryParameter("Missing parameter 'table_file'. It must be specified before listing "
                                         "reactions.")
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


def parse_table(table_name):
    content = tables[table_name]
    first_col = []
    second_col = []
    for line in content:
        cols = line.split()
        if len(cols) != 2:
            raise TableError(f"Syntax error for table '{table_name}' on table row '{line}': "
                             f"Each tables must have exactly two columns.")
        first_col.append(float(cols[0]))
        second_col.append(float(cols[1]))
    return interp1d(first_col, second_col)  # the default is linear interpolation


def check_parameters():
    obligatory_parameters = ['time_end']
    for obligatory_parameter in obligatory_parameters:
        if obligatory_parameter not in parameters:
            raise MissingObligatoryParameter(f"Parameter '{obligatory_parameter}' is not specified.")


def set_default_parameters():
    for species in all_species:
        if species not in parameters:
            parameters[species] = 0
            raise MissingParameter(f"Parameter '{species}' is missing: it is set to 0.")
    if 'time_ini' not in parameters:
        parameters['time_ini'] = 0
        raise MissingParameter(f"Parameter 'time_ini' is missing: it is set to 0.")
    if 'calc_step' not in parameters:
        parameters['calc_step'] = 1
        raise MissingParameter(f"Parameter 'calc_step' is missing: it is set to 1.")

