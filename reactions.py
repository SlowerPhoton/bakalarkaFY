from collections import Counter


class Reaction:
    reactants = {}
    products = {}
    rate_fun = None

    rate_spec = ""

    def __init__(self, reactants, products, rate_fun):
        # convert list to a dictionary, example:
        # ["Ar", "e", "e"] -> {"Ar" : 1, "e" : 2}
        self.reactants = Counter(reactants)
        self.products = Counter(products)
        self.rate_fun = rate_fun

    def compute_a(self, parameters):
        # TODO: must be updated for more complicated instances
        a = self.rate_fun(None)  # this needs to be updated
        for reactant in self.reactants:
            a *= parameters[reactant]  # doesn't cover non-singular reactants
        return a

    def react(self, parameters):
        for reactant in self.reactants:
            parameters[reactant] -= self.reactants[reactant]
        for product in self.products:
            parameters[product] += self.products[product]
        return parameters
