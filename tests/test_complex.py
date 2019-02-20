import pytest
from pyrxn import CRN

def test_complex():
    def nor_gate(i, j, k):
        c = CRN()

        prom = "P{}{}".format(i, j)

        txn = "{0} > {0} + {1}".format(prom, k)
        deg = "{0} > ".format(k)
        r1 = "{0} + {1} <> {0}{1}".format(prom, i)
        r2 = "{0} + {1} <> {0}{1}".format(prom, j)

        c.r(txn, 1.0, 'txn')
        c.r(deg, 0.1, 'deg')
        c.r(r1, [0.01, 0.05], 'repr')
        c.r(r2, [0.01, 0.0005], 'repr')
        c[prom] = 1
        return c

    c1 = nor_gate("r1", "r1", "r5")
    c2 = nor_gate("r2", "r2", "r6")
    c3 = nor_gate("r5", "r6", "r7")
    c4 = CRN.merge_all([c1, c2, c3])

    c4.r(" > r1", 1.0)
    c4.r("r1 > ", 0.1)

    c4.remove_duplicates()
    c4.reduce()

    df = c4.run(1000)