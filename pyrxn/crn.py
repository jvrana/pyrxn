'''

CRN.py
author: Justin Vrana
A Chemical Reaction Network Builder

'''

from pyrxn.exceptions import CRNException
import numpy as np
import pandas as pd
import re
from scipy.integrate import solve_ivp
from collections import defaultdict
from copy import copy


class Reaction:
    def __init__(self, reactants: dict, products: dict, rate_constant: float, name=''):
        """Construct a new reaction"""
        self._relements = reactants
        self._pelements = products
        self.rate = float(rate_constant)
        self.name = name

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
    def from_string(cls, eqn, rate, name=''):
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
                reactions.append(cls.from_element_list(reactants, products, rate, name))
            else:
                reactions.append(cls.from_element_list(products, reactants, rate, name))
        elif direction == "<>" or direction == "=":
            if not (isinstance(rate, list) or isinstance(rate, tuple)):
                raise CRNException("Bidirectional reactions include a list or tuple of rates.")
            elif 0 in rate:
                raise CRNException("Rate cannot be zero. {}".format(rate))
            reactions.append(cls.from_element_list(reactants, products, rate[0], name+"_f"))
            reactions.append(cls.from_element_list(products, reactants, rate[1], name+"_r"))
        else:
            raise CRNException("Direction {} not recognized".format(direction))
        return reactions

    @property
    def elements(self):
        """List of all elements in reaction"""
        elements = []
        for e in self.reactants:
            self._relements = dict(self.reactants)
            if e not in elements:
                elements.append(e)
        for e in self.products:
            self._pelements = dict(self.products)
            if e not in elements:
                elements.append(e)
        return sorted(elements)

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

    def _to_str(self):
        """Helper function used in __hash__ and __repr__"""
        r = [str(self.reactants[k]) + str(k) for k in self.reactants]
        p = [str(self.products[k]) + str(k) for k in self.products]
        s = " + ".join(r) + " > " + " + ".join(p)
        return s

    def __eq__(self, other):
        return self.reactants == other.reactants and \
               self.products == other.products and \
               self.rate == other.rate

    def __hash__(self):
        return hash((self._to_str(), self.rate))

    def __str__(self):
        return self._to_str()

    def __repr__(self):
        return self._to_str()

class CRN:
    def __init__(self, reactions=None, name=''):
        self._reactions = []
        self._reactant_matrix = np.array([])
        self._product_matrix = np.array([])
        self._elements = list()
        self._x = np.array([])
        self.name = name
        if reactions:
            for r in reactions:
                self.r(r)

    @property
    def reactions(self):
        return self._reactions[:]

    def _add_reaction(self, reaction):
        self._reactions.append(reaction)

    def _group_reactions(self):
        # group similar reactions
        all_rxns = defaultdict(list)
        for r in self.reactions:
            all_rxns[r._to_str()].append(r)
        return all_rxns

    def reduce(self):
        grouped = self._group_reactions()
        if any([len(x) > 1 for x in grouped.values()]):
            reduced = []
            for k, rxns in grouped.items():
                copied = rxns[0].copy()
                copied.rate = sum([_r.rate for _r in rxns])
                reduced.append(copied)
            self.set_reactions(reduced)
            self._update()
        return self

    def remove_duplicates(self):
        """
        Remove duplicated reactions.
        :return:
        :rtype:
        """
        self._reactions = list(set(self.reactions))
        self._update()
        return self

    def clear_all(self):
        self._reactions = []
        self._reactant_matrix = np.array([])
        self._product_matrix = np.array([])
        self._elements = list()
        self._x = np.array([])

    def set_reactions(self, reactions, copy=False):
        self.clear_all()
        self.add_reactions(reactions, copy)

    def add_reactions(self, reactions, copy=False):
        """
        Add reactions from a list of Reactions

        :param reactions: list of Reactions
        :type reactions: list
        :return:
        :rtype:
        """
        elements = self._elements[:]
        for r in reactions:
            if copy: r = r.copy()
            self._add_reaction(r)
            elements += r.elements
        self._elements = sorted(list(set(elements)))

        self._update()
        return self

    def r(self, eqn, rate, name=''):
        """
        Adds a new reaction from a string similar to the form `3A + 4B + 2D <> C`

        :param eqn: reaction string
        :type eqn: basestring
        :param rates: list of rates. If direction is '<>', there must be two rates.
        :type rates: float
        :return: None
        :rtype: None
        """
        return self.add_reactions(Reaction.from_string(eqn, rate, name=name))

    def _update(self):
        """
        Updates the state matrices.

        :return:
        :rtype: tuple
        """
        elements = self.E
        stoich_matrix = np.zeros((len(elements), len(self.reactions)))
        reactant_matrix = stoich_matrix.copy()
        product_matrix = stoich_matrix.copy()
        for i, e in enumerate(elements):
            for j, r in enumerate(self.reactions):
                reactant_matrix[i, j] = r.reactants.get(e, 0.0)
                product_matrix[i, j] = r.products.get(e, 0.0)
                stoich_matrix[i, j] = r.products.get(e, 0.0) - r.reactants.get(e, 0.0)
        self._reactant_matrix = np.array(reactant_matrix)
        self._product_matrix = np.array(product_matrix)
        self._elements = elements

        prior_state = self.species
        self.clear_state()
        self.update_state(prior_state)

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
    def x(self):
        """
        Current concentrations in the CRN. Indices correspond to the indices in the elements.

        :return:
        :rtype:
        """
        return self._x

    @property
    def species(self):
        return dict(zip(self.E, self.x.flatten()))

    def clear_state(self):
        self._x = np.zeros((len(self.E), 1))
        return self

    def set_state(self, conc_dict, func=None):
        self.clear_state()
        return self.update_state(conc_dict, func=func)

    def set_state_from_matrix(self, matrix, func=None):
        self.clear_state()
        return self.update_state_from_matrix(matrix, func)

    def update_state_from_matrix(self, matrix, func=None, inplace=True):
        if matrix is None:
            return self.x
        x = self.x
        if matrix.shape != x.shape:
            raise CRNException("Cannot update. Expected a {} but found a {}".format(
                matrix.shape, x.shape))

        for i in range(len(matrix)):
            v = matrix[i, 0]
            if func:
                x[i, 0] = func(x[i, 0], v)
            else:
                x[i, 0] = v
        if inplace:
            self._x = x
            return x.copy()
        return x

    def update_state(self, conc_dict, func=None, inplace=True):
        """
        Updates the concentration matrix from a dictionary.

        :param conc_dict:
        :type conc_dict:
        :param func:
        :type func:
        :return:
        :rtype:
        """
        if conc_dict is None:
            return self.x
        x = self.x
        for k, v in conc_dict.items():
            i = self.E.index(k)
            if func:
                x[i, 0] = func(x[i,0], v)
            else:
                x[i, 0] = v
        if inplace:
            self._x = x
            return x.copy()
        return x

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
        """
        Rate matrix

        :param x:
        :type x:
        :return:
        :rtype:
        """
        if x is None:
            x = self.x
        vmatrix = np.zeros((len(self.reactions), 1))
        for i, r in enumerate(self.reactions):
            val = r.rate
            for e in r.reactants:
                val *= self.get(e, x)**r.reactants[e]
            vmatrix[i, 0] = val
        return vmatrix

    def dx(self, x=None):
        """
        Compute the instantaneous rate.

        :return:
        :rtype:
        """
        if x is None:
            x = self.x
        return np.dot(self.N, self.v(x))



    def solve(self, t_final, t_init=0, init=None, hold_constant=None, *args, **kwargs):
        e_to_i = dict(zip(self.E, range(len(self.E))))

        def dydt(t, y):
            rate = self.v(y.reshape(len(y), 1))
            dydt = np.dot(self.N, rate).flatten()
            if hold_constant:
                for _e in hold_constant:
                    dydt[e_to_i[_e]] = 0
            return dydt

        if init is None:
            init = np.zeros(len(self.E))
        sol = solve_ivp(dydt, [t_init, t_final], init)
        return sol

    def run(self, t_final, t_init=0, init=None, hold_constant=None, *args, **kwargs):
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
        init_vector = self.update_state(init, inplace=False).flatten()
        sol = self.solve(t_final, t_init, init_vector, hold_constant=hold_constant, *args, **kwargs)
        df = pd.DataFrame(sol.y.T, columns=self.E)
        df['time'] = sol.t
        return df

    def to_dataframe(self):
        d = []
        for r in self.reactions:
            d.append([str(r), r.rate, r.name])
        df = pd.DataFrame(d, columns=['reaction', 'rate', 'name'])
        return df

    @classmethod
    def merge_all(cls, crns):
        """
        Merge a list of CRNs

        :param crns:
        :type crns:
        :return:
        :rtype:
        """
        c = CRN()
        for other in crns:
            c.merge(other)
        return c

    def merge(self, other, func=None):
        species = self.species
        self += other
        self.update_state(species)
        self.update_state(other.species, func=func)
        return self

    def copy(self):
        c = self.__class__()
        c.merge(self)
        return c

    def __copy__(self):
        return self.copy()

    def __getitem__(self, element):
        """
        Gets the concentration of the element from the state matrix.

        :param element: element
        :type element: str
        :return: concentration
        :rtype: float
        """
        return self.get(element, self.x)

    def __setitem__(self, e, c):
        """
        Sets the concentration of a species.

        :param e:
        :type e:
        :param c:
        :type c:
        :return:
        :rtype:
        """
        return self.update_state({e: c})

    def __str__(self):
        s = "Chemical Reaction Network:\n"
        s += self.to_dataframe().to_string()
        return s

    def __add__(self, other):
        c = self.__class__()
        c += self
        c += other
        return c

    def __iadd__(self, other):
        self.add_reactions(other.reactions, copy=True)
        return self


