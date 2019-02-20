from pyrxn.crn import CRN


def three_ring():
    c = CRN()
    c.r("PA > PA + A", 1)
    c.r("PB > PB + B", 1)
    c.r("PC > PC + C", 1)
    c.r("PC + 3B <> PC*", [0.5, 0.1])
    c.r("PB + 3A <> PB*", [0.5, 0.1])
    c.r("PA + 3C <> PC*", [0.5, 0.1])
    c.r("A > ", 0.1)
    c.r("B > ", 0.1)
    c.r("C > ", 0.1)
    c.update_state({"PA": 1, "PB": 1, "PC": 1})
    return c