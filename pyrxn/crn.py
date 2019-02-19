'''

CRN.py
author: Justin Vrana
A Chemical Reaction Network Builder

'''

from pyrxn.exceptions import CRNException
import numpy as np
import pandas as pd
import re
import itertools
from scipy.integrate import odeint
from collections import defaultdict
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

class CRN:
    def __init__(self):
        self._reactions = []
        self._reactant_matrix = None
        self._product_matrix = None
        self._elements = None
        self._x = None
        self._init = None

    @property
    def reactions(self):
        return self._reactions[:]

    def _add_reaction(self, reaction):
        self._reactions.append(reaction)

    def r(self, eqn, rate):
        """
        Adds a new reaction from a string similar to the form `3A + 4B + 2D <> C`

        Direction: ['<', '>', '<>']

        :param eqn: reaction string
        :type eqn: basestring
        :param rates: list of rates. If direction is '<>', there must be two rates.
        :type rates: float
        :return: None
        :rtype: None
        """
        reactants, direction, products = re.split("\s*([<>=]+)\s*", eqn)
        reactants = re.split("\s*\+\s*", reactants)
        products = re.split("\s*\+\s*", products)
        while '' in reactants:
            reactants.remove('')
        while '' in products:
            products.remove('')
        if direction in ["<", ">"]:
            if not (isinstance(rate, int) or isinstance(rate, float)):
                raise CRNException("Rate must by a number for directional reactions. Instead found a {}".format(type(rate)))
            elif rate == 0:
                    raise CRNException("Rate cannot be zero.")
            if direction == ">":
                self._add_reaction(Reaction(reactants, products, rate))
            else:
                self._add_reaction(Reaction(products, reactants, rate))
        elif direction == "<>" or direction == "=":
            if not (isinstance(rate, list) or isinstance(rate, tuple)):
                raise CRNException("Bidirectional reactions include a list or tuple of rates.")
            elif 0 in rate:
                raise CRNException("Rate cannot be zero. {}".format(rate))
            self._add_reaction(Reaction(reactants, products, rate[0]))
            self._add_reaction(Reaction(products, reactants, rate[1]))
        else:
            raise CRNException("Direction {} not recognized".format(direction))
        self._stoich()

    def _all_elements(self):
        elements = [r.elements for r in self.reactions]
        elements = list(set(itertools.chain(*elements)))
        elements.sort()
        return elements

    def _stoich(self):
        elements = self._all_elements()
        stoich_matrix = np.zeros((len(elements), len(self.reactions)))
        reactant_matrix = stoich_matrix.copy()
        product_matrix = stoich_matrix.copy()
        for i, e in enumerate(elements):
            for j, r in enumerate(self.reactions):
                reactant_matrix[i,j] = r.reactants[e]
                product_matrix[i,j] = r.products[e]
                stoich_matrix[i,j] = r.products[e] - r.reactants[e]
        self._reactant_matrix = np.array(reactant_matrix)
        self._product_matrix = np.array(product_matrix)
        self._elements = elements
        return reactant_matrix, product_matrix, elements

    @property
    def R(self):
        return self._reactant_matrix

    @property
    def P(self):
        return self._product_matrix

    @property
    def N(self):
        return self.P - self.R

    @property
    def E(self):
        return list(self._elements)

    @property
    def x(self):
        return self._x

    @property
    def init(self):
        return self._init

    def get(self, matrix, element):
        index = self.E.index(element)
        return matrix[index, 0]

    def v(self, x):
        vmatrix = np.zeros((len(self.reactions), 1))
        for i, r in enumerate(self.reactions):
            val = r.rate
            for e in r.reactants:
                val *= self.get(x, e)**r.reactants[e]
            vmatrix[i, 0] = val
        return vmatrix

    def initialize(self, init):
        """
        Initialize CRN with starting concentrations

        :param init: dictioanry of starting concentrations
        :type init: dict
        :return: None
        :rtype: None
        """
        self._x = np.zeros((len(self.E), 1))
        for k, v in init.items():
            i = self.E.index(k)
            self._x[i, 0] = v
        self._init = self._x.copy()

    def get_rxn_flux_vector(self):
        """
        Returns the reaction flux vector

        :return: reaction flux vector
        :rtype: np.array
        """
        rxn_flux = []
        for j, reaction in enumerate(self.reactions):
            rate = 1
            for i, element in enumerate(self.elements):
                stoich_ele_i = reaction.reactant_elements[element]
                rate = rate * self.state[i]**stoich_ele_i
            rxn_flux.append(rate*reaction.rate_constant)
        self.reaction_flux_vector = np.array(rxn_flux)
        return self.reaction_flux_vector

    def apply_flux_vector(self, state):
        """
        Get the instantan
        :param state:
        :type state:
        :return:
        :rtype:
        """
        rxn_flux = []
        for j, reaction in enumerate(self.reactions):
            rate = 1
            for i, element in enumerate(self.elements):
                stoich_ele_i = reaction.reactant_elements[element]
                rate = rate * state[i]**stoich_ele_i
            rxn_flux.append(rate*reaction.rate_constant)
        return np.array(rxn_flux)

    def dX(self):
        self.get_rxn_flux_vector()
        return np.dot(self.N, self.v(self.x))

    def run(self, dt, t_final, init=None):
        if init is None:
            init = self.init.copy()

        def func(y, t):
            return np.dot(self.N, self.v(y.reshape(len(y), 1))).flatten()

        y0 = init.flatten()
        t = np.linspace(0, t_final, t_final/dt)
        y = odeint(func, y0, t)

        y = pd.DataFrame(y)
        y.index = y.index * dt
        y.columns = self._elements
        y.index.name = 'time'
        return y

    def __add__(self, other):
        c = self.__class__()
        c += self
        c += other
        return c

    def __iadd__(self, other):
        for r in other.reactions:
            self._add_reaction(r.copy())
        return self

    def dose_response(self, element, values, dt, t_final):
        index = list(self.elements).index(element)
        array = []
        for i in values:
            init_copy = self.init.copy()
            init_copy[index] = i
            d = self.run(dt, t_final, init=init_copy)
            array.append(d.iloc[-1])
        array = pd.DataFrame(array)
        array.index = values
        array.index.name = "{}_T".format(element)
        return array

