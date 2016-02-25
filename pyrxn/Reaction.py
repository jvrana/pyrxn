from collections import defaultdict
import re

class Reaction:
    def __init__(self, reactants, products, rate_constant):
        self.elements = []
        self.reactant_elements = defaultdict(int)
        self.product_elements = defaultdict(int)
        self.rate_constant = rate_constant
        for r in reactants:
            s, e = self.convert_element(r)
            self.reactant_elements[e] = s
            self.elements.append(e)
        for p in products:
            s, e = self.convert_element(p)
            self.product_elements[e] = s
            self.elements.append(e)
        self.elements = list(set(self.elements))

    def convert_element(self, e):
        g = re.match("(?P<stoich>\d?)(?P<ele>[a-zA-Z0-9\*:]+)", e)
        stoich = g.group('stoich')
        if stoich == '':
            stoich = 1
        else:
            stoich = int(stoich)
        element = g.group('ele')
        return stoich, element

    def __repr__(self):
        s = "Reactant(" + str(self.reactant_elements) + " > " + str(self.product_elements) +")"
        return s