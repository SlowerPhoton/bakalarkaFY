from collections import Counter
from math import factorial


class Reaction:
    reactants = {}
    products = {}
    rate_fun = None

    table_name = None
    ratio = 1

    def __init__(self, reactants, products, rate_fun=None):
        # convert list to a dictionary, example:
        # ["Ar", "e", "e"] -> {"Ar" : 1, "e" : 2}
        self.reactants = Counter(reactants)
        self.products = Counter(products)
        self.rate_fun = rate_fun

    def compute_a(self, parameters):
        a = self.rate_fun(parameters) * self.ratio
        for reactant in self.reactants:
            for i in range(self.reactants[reactant]):
                a *= parameters[reactant] - i
            a *= 1/factorial(self.reactants[reactant])
        return a

    def scale(self, parameters, overwrite_ratio=True):
        if self.rate_fun is None:
            return
        scale = (sum(self.reactants.values()) - 1) * 3
        if self.table_name in parameters:
            scale = parameters[self.table_name]
        if overwrite_ratio:
            self.ratio = 1
        self.ratio *= parameters["ratio"] ** scale

    def react(self, parameters):
        for reactant in self.reactants:
            parameters[reactant] -= self.reactants[reactant]
        for product in self.products:
            parameters[product] += self.products[product]
        return parameters
