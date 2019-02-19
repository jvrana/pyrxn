from collections import defaultdict
import re
from copy import copy

class Reaction:
    def __init__(self, reactants, products, rate_constant, name=''):
        """Construct a new reaction"""
        self._elements = list()
        self._relements = defaultdict(int)
        self._pelements = defaultdict(int)
        self.rate = rate_constant
        self.name = name

        for r in reactants:
            s, e = self._convert_element(r)
            self._relements[e] = s
            if e not in self._elements:
                self._elements.append(e)
        for p in products:
            s, e = self._convert_element(p)
            self._pelements[e] = s
            if e not in self._elements:
                self._elements.append(e)
        sorted(self._elements)

    @property
    def elements(self):
        """List of all elements in reaction"""
        return list(self._elements)

    @property
    def reactants(self):
        """List of all reactants"""
        return copy(self._relements)

    @property
    def products(self):
        """List of all products"""
        return copy(self._pelements)

    def _convert_element(self, e):
        g = re.match("(?P<stoich>\d?)(?P<ele>[a-zA-Z0-9\*:]+)", e)
        stoich = g.group('stoich')
        if stoich == '':
            stoich = 1
        else:
            stoich = int(stoich)
        element = g.group('ele')
        return stoich, element

    def copy(self, suffix=None):
        reaction = Reaction.__new__(Reaction)

        relements = defaultdict(int)
        pelements = defaultdict(int)

        # rename elements if i is provided
        for k, v in self._relements.items():
            if suffix:
                k += str(suffix)
            relements[k] = v

        for k, v in self._pelements.items():
            if suffix:
                k += str(suffix)
            pelements[k] = v

        reaction._relements = relements
        reaction._pelements = pelements
        reaction.rate = self.rate
        name = self.name
        if suffix:
            name += "_" + str(suffix)
        reaction.name = name
        return reaction


    def __repr__(self):
        r = [str(self.reactants[k]) + str(k) for k in self.reactants]
        p = [str(self.products[k]) + str(k) for k in self.products]
        s = " + ".join(r) + " > " + " + ".join(p)
        return s