"""
This file contains the class Reaction and specifies its methods.
"""
from collections import Counter


class Reaction:
    """
    When instantiated from the input_parser, the attribute rate_fun will be None only if it cannot parse the string
    following the character '!' after specifying the reaction in the input file.
    Reactants and products are always known. They are dictionaries with keys being a subset of species and values being
    natural numbers. E.g., reactants = {"Ar" : 1, "e" : 2}.

    Attribute table_name is set by the input_parser if the reaction rate is described by a table.
    """
    reactants = {}
    products = {}
    rate_fun = None

    table_name = None

    def __init__(self, reactants, products, rate_fun=None):
        # Counter: convert list to a dictionary, example:
        # ["Ar", "e", "e"] -> {"Ar" : 1, "e" : 2}
        self.reactants = Counter(reactants)
        self.products = Counter(products)
        self.rate_fun = rate_fun

    def compute_a(self, parameters):
        """
        Compute transition rate a_mu of the reaction R_mu.
        """
        a = self.rate_fun(parameters)
        for reactant in self.reactants:
            for i in range(self.reactants[reactant]):
                a *= parameters[reactant] - i
        if a < 0: a = 0  # prevent negative transition rates
        return a

    def react(self, parameters, bulk):
        """
        Simulate executing bulk reactions: i.e., updating species concentrations.
        """
        for reactant in self.reactants:
            parameters[reactant] -= self.reactants[reactant] * bulk
        for product in self.products:
            parameters[product] += self.products[product] * bulk
        return parameters
