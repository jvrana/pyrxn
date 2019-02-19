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
from scipy.integrate import solve_ivp
from collections import defaultdict
from copy import copy


class Reaction:
    def __init__(self, reactants: dict, products: dict, rate_constant: float, name=''):
        """Construct a new reaction"""
        self._elements = list()
        self._relements = reactants
        self._pelements = products
        self.rate = float(rate_constant)
        self.name = name

        for e in reactants:
            self._relements = dict(reactants)
            if e not in self._elements:
                self._elements.append(e)
        for e in products:
            self._pelements = dict(products)
            if e not in self._elements:
                self._elements.append(e)
        sorted(self._elements)

    @classmethod
    def from_element_list(cls, reactants, products, rate_contant, name=''):
        """
        Construct reaction from a list of reactants and products.

        :param reactants: e.g ["1A", "2B"]
        :type reactants: str
        :param products: e.g ["1A", "2B"]
        :type products: str
        :param rate_contant: rate constant for this reaction
        :type rate_contant: float
        :param name: optional name of reaction
        :type name: str
        :return: new reaction
        :rtype: Reaction
        """
        rcoeff = dict()
        pcoeff = dict()
        for r in reactants:
            s, e = cls._convert_element(r)
            rcoeff[e] = s
        for p in products:
            s, e = cls._convert_element(p)
            pcoeff[e] = s
        return cls(rcoeff, pcoeff, rate_contant, name)

    @classmethod
    def from_string(cls, eqn, rate):
        """
        Construct reactions from unidirectional or bidirectional reactions.

        :param eqn: e.g. " A + B <> 2C"
        :type eqn: string
        :param rate: rate constant(s). Tuple must be provided if bidirectional reaction.
        :type rate: list or float
        :return:
        :rtype:
        """
        reactants, direction, products = re.split("\s*([<>=]+)\s*", eqn)
        reactants = re.split("\s*\+\s*", reactants)
        products = re.split("\s*\+\s*", products)
        reactions = []
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
                reactions.append(cls.from_element_list(reactants, products, rate))
            else:
                reactions.append(cls.from_element_list(products, reactants, rate))
        elif direction == "<>" or direction == "=":
            if not (isinstance(rate, list) or isinstance(rate, tuple)):
                raise CRNException("Bidirectional reactions include a list or tuple of rates.")
            elif 0 in rate:
                raise CRNException("Rate cannot be zero. {}".format(rate))
            reactions.append(cls.from_element_list(reactants, products, rate[0]))
            reactions.append(cls.from_element_list(products, reactants, rate[1]))
        else:
            raise CRNException("Direction {} not recognized".format(direction))
        return reactions

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

    @staticmethod
    def _convert_element(e):
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
    def __init__(self, reactions=None):
        self._reactions = []
        self._reactant_matrix = None
        self._product_matrix = None
        self._elements = list()
        self._x = None
        if reactions:
            for r in reactions:
                self.r(r)

    @property
    def reactions(self):
        return self._reactions[:]

    def _add_reaction(self, reaction):
        self._reactions.append(reaction)

    def add_reactions(self, reactions):
        """
        Add reactions from a list of Reactions

        :param reactions: list of Reactions
        :type reactions: list
        :return:
        :rtype:
        """
        elements = self._elements[:]
        for r in reactions:
            print(r)
            self._add_reaction(r)
            elements += r.elements
        self._elements = sorted(elements)
        self._stoich()
        self.set_state()
        return self

    def r(self, eqn, rate):
        """
        Adds a new reaction from a string similar to the form `3A + 4B + 2D <> C`

        :param eqn: reaction string
        :type eqn: basestring
        :param rates: list of rates. If direction is '<>', there must be two rates.
        :type rates: float
        :return: None
        :rtype: None
        """
        return self.add_reactions(Reaction.from_string(eqn, rate))

    def _stoich(self):
        """
        Initialize the stoichiometry matrix

        :return:
        :rtype:
        """
        elements = self.elements
        stoich_matrix = np.zeros((len(elements), len(self.reactions)))
        reactant_matrix = stoich_matrix.copy()
        product_matrix = stoich_matrix.copy()
        print(elements)
        for i, e in enumerate(elements):
            for j, r in enumerate(self.reactions):
                reactant_matrix[i,j] = r.reactants.get(e, 0.0)
                product_matrix[i,j] = r.products.get(e, 0.0)
                stoich_matrix[i,j] = r.products.get(e, 0.0) - r.reactants.get(e, 0.0)
        self._reactant_matrix = np.array(reactant_matrix)
        self._product_matrix = np.array(product_matrix)
        self._elements = elements
        return reactant_matrix, product_matrix, elements

    @property
    def R(self):
        """
        The reactant component of the stoichiometry matrix. A m x n matrix where m is the
        number of elements and n is the number of reactions.

        :return:
        :rtype:
        """
        return self._reactant_matrix

    @property
    def P(self):
        """
        The product component of the stoichiometry matrix. A m x n matrix where m is the
        number of elements and n is the number of reactions.

        :return:
        :rtype:
        """
        return self._product_matrix

    @property
    def N(self):
        """
        The stoichiometry matrix. A m x n matrix where m is the
        number of elements and n is the number of reactions.

        :return:
        :rtype:
        """
        return self.P - self.R

    @property
    def E(self):
        """
        Returns list of elements in the CRN

        :return:
        :rtype:
        """
        return list(self._elements)

    @property
    def elements(self):
        return list(self._elements)

    @property
    def x(self):
        """
        Current concentrations in the CRN. Indices correspond to the indices in the elements.

        :return:
        :rtype:
        """
        return self._x

    def __getitem__(self, element):
        """
        Gets the concentration of the element from the state matrix.

        :param element: element
        :type element: str
        :return: concentration
        :rtype: float
        """
        return self.get(element, self.x)

    def get(self, element, matrix, default=None):
        """
        Gets the concetration of element from a matrix.

        :param matrix:
        :type matrix:
        :param element:
        :type element:
        :return:
        :rtype:
        """
        try:
            index = self.E.index(element)
        except ValueError:
            return default
        return matrix[index, 0]

    def v(self, x=None):
        if x is None:
            x = self.x
        vmatrix = np.zeros((len(self.reactions), 1))
        for i, r in enumerate(self.reactions):
            val = r.rate
            for e in r.reactants:
                val *= self.get(e, x)**r.reactants[e]
            vmatrix[i, 0] = val
        return vmatrix

    def set_state(self, init=None):
        """
        Initialize CRN with starting concentrations

        :param init: dictioanry of starting concentrations
        :type init: dict
        :return: None
        :rtype: None
        """
        self._x = np.zeros((len(self.E), 1))
        if init is not None:
            self._x = self._state_matrix(init)

    def dX(self, x=None):
        """
        Compute the instantaneous rate.

        :return:
        :rtype:
        """
        if x is None:
            x = self.x
        self.get_rxn_flux_vector()
        return np.dot(self.N, self.v(x))

    def _state_matrix(self, init=None):
        init_vector = self.x
        size = len(self.E)
        if isinstance(init, dict):
            for k, v in init.items():
                i = self.E.index(k)
                init_vector[i, 0] = v
        elif isinstance(init, list):
            init_vector = np.array(init).reshape((size, 1))
        elif isinstance(init, np.ndarray):
            if init.ndim == 1:
                init_vector = init
                init_vector.reshape((size, 1))
            elif init.ndim == 2:
                init_vector = init
        if len(init_vector) != size:
            raise CRNException("Initialization vector size must be the same size as num elements ({})".format(size))
        return init_vector

    def solve(self, t_final, t_init=0, init=None, *args, **kwargs):
        def func(t, y):
            return np.dot(self.N, self.v(y.reshape(len(y), 1))).flatten()

        if init is None:
            init = np.zeros(len(self.E))
        print(init)
        sol = solve_ivp(func, [t_init, t_final], init)
        return sol

    def run(self, t_final, t_init=0, init=None, *args, **kwargs):
        """
        Run ode using solve_ivp from scipy.integrate

        :param dt:
        :type dt:
        :param t_final:
        :type t_final:
        :param init:
        :type init:
        :return:
        :rtype:
        """
        init_vector = self._state_matrix(init).flatten()
        sol = self.solve(t_final, t_init, init_vector, *args, **kwargs)
        df = pd.DataFrame(sol.y.T, columns=self.E)
        df['time'] = sol.t
        return df

    def __add__(self, other):
        c = self.__class__()
        c += self
        c += other
        return c

    def __iadd__(self, other):
        self.add_reactions(other.reactions)
        return self

    # def dose_response(self, element, values, dt, t_final):
    #     index = list(self.elements).index(element)
    #     array = []
    #     for i in values:
    #         init_copy = self.init.copy()
    #         init_copy[index] = i
    #         d = self.run(dt, t_final, init=init_copy)
    #         array.append(d.iloc[-1])
    #     array = pd.DataFrame(array)
    #     array.index = values
    #     array.index.name = "{}_T".format(element)
    #     return array

