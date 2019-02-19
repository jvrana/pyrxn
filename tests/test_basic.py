from pyrxn import CRN
from pyrxn.exceptions import CRNException
import pytest
import numpy as np

def test_basic():
    c = CRN()
    c.r("A + 2B <> C", [1.0, 2.0])

    assert(2 == len(c.reactions))

def test_basic2():
    c = CRN()
    c.r("A + 2B > 3C", 1.0)

    assert(1 ==len(c.reactions))

@pytest.mark.parametrize("reactants", ["A", "A + B", "A + B + C"])
@pytest.mark.parametrize("products", ["A", "A + 3B", "A + B + 2C"])
@pytest.mark.parametrize("direction", [">", "<"])
@pytest.mark.parametrize("rate", [0.1, 0.0, 1.0])
@pytest.mark.parametrize("num_rxns", [1, 2, 10])
def test_unidirectional(reactants, products, direction, rate, num_rxns):
    rxn = reactants + direction + products
    c = CRN()

    if rate == 0:
        with pytest.raises(CRNException):
            c.r(rxn, rate)
    else:
        for i in range(num_rxns):
            c.r(rxn, rate)
        assert(num_rxns == len(c.reactions))


@pytest.mark.parametrize("reactants", ["A", "A + B", "A + B + C"])
@pytest.mark.parametrize("products", ["A", "A + B", "A + B + C"])
@pytest.mark.parametrize("direction", ["=", "<>"])
@pytest.mark.parametrize("rates", [[1.0, 1.0], [2.0, 10], [0, 1]])
@pytest.mark.parametrize("num_rxns", [1, 2, 10])
def test_bidirectional(reactants, products, direction, rates, num_rxns):
    rxn = reactants + direction + products
    c = CRN()

    if 0 in rates:
        with pytest.raises(CRNException):
            c.r(rxn, rates)
    else:
        for i in range(num_rxns):
            c.r(rxn, rates)
        assert(2*num_rxns == len(c.reactions))

def test_N():
    c = CRN()
    c.r("A + 2B <> C", [1.0, 2.0])
    print(c.E)
    print(c.N)
    print(c.R)
    print(c.P)

    c.initialize(dict(C=55))
    print(c.x)
    print("v(x)")
    x = np.ones((3,1))*5
    print(c.v(x))

    d = c.run(0.001, 20)
    print(d)


def test_copy_reaction():
    c = CRN()
    c.r("2A + B > C", 1.0)

    reaction = c.reactions[0]
    r2 = reaction.copy(1)

    print(r2)