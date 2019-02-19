from pyrxn.utilities.plotting import plot, gradient

def _op_dict(d1, d2, func):
    d3 = dict()

    allkeys = set(d1.keys()).union(d2.keys())

    for k in allkeys:
        d3[k] = func(d1.get(k, 0), d2.get(k, 0))
    return d3


def sum_dicts(d1, d2):
    return _op_dict(d1, d2, lambda a, b: a + b)


def subtract_dicts(d1, d2):
    return _op_dict(d1, d2, lambda a, b: a - b)

