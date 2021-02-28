from reactions import Reaction
from scipy.interpolate import interp1d
import warnings

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


def parse_input_file(filename):
    all_species = set()
    parameters = {}
    reactions = []
    with open(filename) as file:
        for line in file.readlines():
            line = line.strip()  # get rid of whitespace L and R
            if line == "":  # skip empty lines
                continue
            if line[0] == "#":  # skip comments
                continue
            parse_line(line, all_species, parameters, reactions)
    check_parameters()  # check for any missing obligatory parameter
    set_default_parameters()  # set other missing parameters to their default values
    return all_species, parameters, reactions


def parse_line(line, all_species, parameters, reactions):
    if '=>' in line:
        parse_reaction(line, all_species, reactions)
    elif '=' in line:
        parse_parameter(line, parameters)
    else:
        raise InvalidLine(f"There is a syntax error on line '{line}'. Each line must be empty, start "
                          f"with the hashtag ('#'), set up a parameter or specify a reaction.")


def parse_parameter(line, parameters):
    line_split = line.split('=')
    if len(line_split) != 2:
        raise InvalidLine(f"Syntax error on line '{line}': Parameters are set using the '=' equals sign. There can "
                          f"only be one parameter set per one line.")
    parameter = line_split[0].strip()
    try:
        value = float(line_split[1]) if "." in line_split[1] else int(line_split[1])
    except ValueError:
        value = line_split[1].strip()
    parameters[parameter] = value


def parse_reaction(line, all_species, reactions):
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

    reaction = Reaction(parse_species(left_side, all_species),
                        parse_species(right_side, all_species),
                        parse_rate(rate_spec))
    # TODO: remove code repetition
    reaction.rate_spec = rate_spec  # remember rate specification for debugging purposes
    if rate_spec.startswith("table:"):
        reaction.table = rate_spec[len("table:"):].strip()
    reactions.append(reaction)


def parse_species(line, all_species):
    species = line.split()
    # + can be used as a delimiter but it is unnecessary
    # species must be separated by whitespace
    # '+' delimited by whitespace are ignored, otherwise they count as part of the species name
    while '+' in species:
        species.remove('+')
    all_species.update(species)
    return species


def parse_rate(rate_spec, ratio=1):
    fun = None
    try:
        num = float(rate_spec)
        fun = lambda energy: num*ratio
    except ValueError:
        rate_spec = rate_spec.strip()
        if rate_spec.startswith("table:"):
            table_name = rate_spec[len("table:"):].strip()
            if len(tables) == 0:
                parse_tables()
            fun = parse_table(table_name, ratio=ratio)
        else:
            warnings.warn(f"The rate specification '{rate_spec}' cannot be parsed, because we only support "
                          f"constant rates and rates given by a table. You must write the rate function in "
                          f"the code yourself and assign it to the corresponding reaction.", UserWarning)
    return fun


tables = {}


def parse_tables(parameters=parameters):
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
        return tables


def parse_table(table_name, ratio=1):
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
    return lambda parameters: f(parameters["EN"]) * ratio


def parse_table_custom(table_name):
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


def check_parameters():
    obligatory_parameters = ['time_end']
    for obligatory_parameter in obligatory_parameters:
        if obligatory_parameter not in parameters:
            raise MissingObligatoryParameter(f"Parameter '{obligatory_parameter}' is not specified.")


def set_default_parameters():
    for species in all_species:
        if species not in parameters:
            parameters[species] = 0
            warnings.warn(f"Parameter '{species}' is missing: it is set to 0.", UserWarning)
    if 'time_ini' not in parameters:
        parameters['time_ini'] = 0
        warnings.warn(f"Parameter 'time_ini' is missing: it is set to 0.", UserWarning)
    if 'calc_step' not in parameters:
        parameters['calc_step'] = 1
        warnings.warn(f"Parameter 'calc_step' is missing: it is set to 1.", UserWarning)


def extract_columns(text, cols_to_extract):
    ret = ""
    for line in text.splitlines():
        line = line.strip()
        if line == "":
            continue
        cols = line.split()
        for col in cols_to_extract:
            ret += cols[col]
            ret += "\t"
        ret += "\n"
    return ret


def get_table_scales(parameters, reactions, tables):
    table_scales = {}
    for reaction in reactions:
        if reaction.table is not None:
            table_scales[reaction.table] = (sum(reaction.reactants.values()) - 1) * 3
    for table_name in tables:
        if table_name in parameters:
            table_scales[table_name] = parameters[table_name]
    return table_scales
