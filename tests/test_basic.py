from pyrxn import CRN
from pyrxn.exceptions import CRNException
import pytest
import numpy as np
from pyrxn.utilities import plot, gradient


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


def test_magic_add():
    c1 = CRN()
    c1.r("2A + B > C", 1.0)

    c2 = CRN()
    c1.r("3A + B <> C", [1.0, 3.0])

    c3 = c1 + c2
    assert(len(c3.reactions) == 3)

    for r in c3.reactions:
        print(r)


def test_run():
    c = CRN()

    c.r("P1 > m1 + P1", 0.01)
    c.r("m1 > ", 0.01)
    c.r("P1 + m2 <> PR1", [0.001, 0.001])

    c.initialize({"P1": 1})

    df = c.run(1, 750)
    df['time'] = df.index

    gradient(df, 'time', 1)
    gradient(df, 'time', 2)

    plot(df, 'time', c.E)



