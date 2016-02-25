'''

CRN.py
author: Justin Vrana
A Chemical Reaction Network Builder

'''

from .Reaction import Reaction
from collections import defaultdict
import numpy as np
import pandas as pd
import re
import itertools
from scipy.integrate import odeint

class CRN:
    def __init__(self):
        self.reactions = []

    def r(self, eqn, *rates):
        ''' Adds reaction based on string format A + B <> C '''
        reactants, direction, products = re.split("\s*([<>]+)\s*", eqn)
        reactants = re.split("\s*\+\s*", reactants)
        products = re.split("\s*\+\s*", products)
        while '' in reactants:
            reactants.remove('')
        while '' in products:
            products.remove('')
        if re.match("[<>]+", direction) is not None:
            if ">" in direction:
                self.reactions.append(Reaction(reactants, products, rates[0]))
            if "<" in direction:
                self.reactions.append(Reaction(products, reactants, rates[1]))
        else:
            print "direction not recognized"

    def initialize(self, init):
        elements_array = [r.elements for r in self.reactions]
        elements_array = list(set(itertools.chain(*elements_array)))
        elements_array.sort()
        self.elements = np.array(elements_array)
        self.elements.setflags(write=False) # Read only

        init = defaultdict(float, init)
        self.state = np.array([init[e] for e in self.elements])
        self.init = self.state.copy()
        self.init.setflags(write=False) # Read only
        self.state.setflags(write=False) # Read only
        reactant_matrix = []
        product_matrix = []
        for r in self.reactions:
            reactant_array = []
            product_array = []
            for e in self.elements:
                reactant_array.append(r.reactant_elements[e])
                product_array.append(r.product_elements[e])
            reactant_matrix.append(reactant_array)
            product_matrix.append(product_array)
        self.reactant_matrix = np.array(reactant_matrix)
        self.product_matrix = np.array(product_matrix)

    def get_rxn_flux_vector(self):
        rxn_flux = []
        for j, reaction in enumerate(self.reactions):
            rate = 1
            for i, element in enumerate(self.elements):
                stoich_ele_i = reaction.reactant_elements[element]
                rate = rate * self.state[i]**stoich_ele_i
            rxn_flux.append(rate*reaction.rate_constant)
        self.reaction_flux_vector = np.array(rxn_flux)

    def apply_flux_vector(self, state):
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
        return np.dot(self.reaction_flux_vector, c.product_matrix - c.reactant_matrix)

    def _get_index(self, element):
        return list(self.elements).index(element)

    def run(self, dt, t_final, init=None):
        if init is None:
            init = self.init.copy()

        def func(y, t):
            flux_vector = self.apply_flux_vector(y)
            dx = np.dot(flux_vector, c.product_matrix - c.reactant_matrix)
            return dx

        y0 = init
        t = np.linspace(0,t_final,t_final/dt)
        y = odeint(func, y0, t)
        y = pd.DataFrame(y)
        y.index = y.index * dt
        y.columns = self.elements
        y.index.name = 'time'
        return y

    def dose_response(self, element, values, dt, t_final):
        index = list(c.elements).index(element)
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