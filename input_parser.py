from reactions import Reaction

parameters = {}
reactions = []


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


class UnsupportedFeature(InputFileError):
    """Exception raised for a feature that is not yet supported but it is called in the input file.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


def parse_file(filename):
    with open(filename) as file:
        for line in file.readlines():
            line = line.strip()  # get rid of whitespace L and R
            if line == "":  # skip empty lines
                continue
            if line[0] == "#":  # skip comments
                continue
            parse_line(line)
    check_parameters(parameters)
    return parameters, reactions


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
    return species


def parse_rate(rate_spec):
    try:
        num = float(rate_spec)
    except ValueError:
        raise UnsupportedFeature(f"The rate specification '{rate_spec}' cannot be parsed, because we only support "
                                 f"constant rates.")
    return (lambda energy: num)


def check_parameters(parameters):
    # TODO
    pass
