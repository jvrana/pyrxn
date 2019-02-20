from pyrxn import CRN, Reaction
from pyrxn.exceptions import CRNException
import pytest
import numpy as np
from pyrxn.utilities import plot, gradient
from copy import copy

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

    c.set_state(dict(C=55))
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


def test_reaction__eq__():
    r1 = Reaction.from_string("A + B > C", 1.0, name="Reaction1")
    r2 = Reaction.from_string("A + B > C", 1.0, name="Reaction2")
    r3 = Reaction.from_string("A + B > C", 0.3)
    r4 = Reaction.from_string("A > C", 1.0)
    r5 = Reaction.from_string("A + B > D", 1.0)

    assert(r1 == r2), "Should be equivalent"
    assert(r1 != r3), "Reactions are not equal due to different rates"
    assert(r1 != r4), "Reactants are different"
    assert(r1 != r5), "Products are different"


def test_set_operator():
    c = CRN()
    c.r("A + B <> 2C", [1, 1])
    c['A'] = 1.1

    assert(list(c.x.flatten()) == [1.1, 0, 0])


def test_add_operator():
    """Expect to remove duplicates"""

    eqn1 = "A + B > C"
    eqn2 = "C <> A + B"

    c1 = CRN()
    c1.r(eqn1, 0.01)
    c2 = CRN()
    c2.r(eqn2, [1, 0.01])
    c3 = c1 + c2

    assert(len(c3.reactions) == 3)
    c3.remove_duplicates()
    assert(len(c3.reactions) == 2)


def test_reduce():
    """Expect to remove duplicates"""

    eqn1 = "A + B > C"
    eqn2 = "C <> A + B"

    c1 = CRN()
    c1.r(eqn1, 0.01)
    c2 = CRN()
    c2.r(eqn2, [1, 0.01])
    c3 = c1 + c2

    assert(len(c3.reactions) == 3)
    c3.reduce()
    assert(len(c3.reactions) == 2)
    assert(c3.reactions[0].rate == 0.02)


def test_copy():
    c = CRN()
    c.r("A + B <> 2C", [1, 1])
    c2 = c.copy()
    assert(id(c2) != id(c))
    assert(id(c2.reactions[0]) != id(c.reactions[0]))
    assert(c2.reactions[0] == c.reactions[0])

def test_magic_copy():
    c = CRN()
    c.r("A + B <> 2C", [1, 1])
    c2 = copy(c)
    assert(id(c2) != id(c))
    assert(id(c2.reactions[0]) != id(c.reactions[0]))
    assert(c2.reactions[0] == c.reactions[0])


def test_merge():
    c1 = CRN()
    c1.r("A + B <> 2C", [1, 1])
    c1['A'] = 2

    c2 = copy(c1)
    c2['A'] = 4
    c2['B'] = 1

    assert(c1['A'] == 2)
    assert(c1['B'] == 0)
    assert(c1['C'] == 0)

    c1.merge(c2)

    assert(c1['A'] == 4)
    assert(c1['B'] == 1)
    assert(c1['C'] == 0)


def test_merge_maintains_part():
    c = CRN()
    c.r("A + B > 2C", 1, part='part1')

    c2 = CRN()
    c2.r("AB > BC", 1, part='part2')


    c3 = CRN.merge_all([c, c2])
    for r in c3.reactions:
        print(r)
    assert c3.reactions[0].part == 'part1'
    assert c3.reactions[1].part == 'part2'

    print(c3)


def test_get_part():
    c = CRN()
    c.r("A + B > 2C", 1, part='part1')

    c2 = CRN()
    c2.r("AB <> BC", [1,1], part='part2')


    c3 = CRN.merge_all([c, c2])

    c4 = c3.get_part('part1')
    c5 = c3.get_part('part2')
    c6 = c3.get_part('part1', 'part2')

    assert len(c3.reactions) == 3
    assert len(c4.reactions) == 1
    assert len(c5.reactions) == 2
    assert len(c6.reactions) == 3


def test_merge_with_add():
    c1 = CRN()
    c1.r("A + B <> 2C", [1, 1])
    c1['A'] = 2

    c2 = copy(c1)
    c2['A'] = 4
    c2['B'] = 1

    assert(c1['A'] == 2)
    assert(c1['B'] == 0)
    assert(c1['C'] == 0)
    print(c1.x)
    c1.merge(c2, lambda a, b: a+b)

    assert(c1['A'] == 6)
    assert(c1['B'] == 1)
    assert(c1['C'] == 0)

def test_merge_all():
    c1 = CRN()
    c1.r("A + B <> 2C", [1, 1])
    c1['A'] = 2

    c2 = copy(c1)
    c2['A'] = 4
    c2['B'] = 1

    assert(c1['A'] == 2)
    assert(c1['B'] == 0)
    assert(c1['C'] == 0)
    c3 = CRN.merge_all([c1, c2])
    assert(c3['A'] == 4)
    assert(c3['B'] == 1)
    assert(c3['C'] == 0)


def test_add_reaction_shoud_maintain_state():

    c1 = CRN()
    c1.r("A + B > C", 1)
    c1['A'] = 1.0

    c1.r('B > C', 2)

    assert(c1['A'] == 1)

def test_run():
    c = CRN()

    c.r("P1 > m1 + P1", 0.01)
    c.r("m1 > ", 0.01)
    c.r("P1 + m2 <> PR1", [0.001, 0.001])

    c.set_state({"P1": 1})

    df = c.run(750)

    gradient(df, 'time', 1)
    gradient(df, 'time', 2)

    plot(df, 'time', c.E)



