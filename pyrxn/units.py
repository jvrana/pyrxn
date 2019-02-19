# 2016/03/26

import quantities as pq

def getUnits():
    pq.AU = pq.UnitQuantity('arbitrary units', 1, symbol='AU')
    pq.uL = pq.UnitQuantity('microliter', 0.001 * pq.mL, symbol='uL')
    pq.nL = pq.UnitQuantity('nanoliter', 0.001 * pq.uL, symbol='nL')
    pq.pL = pq.UnitQuantity('picoliter', 0.001 * pq.nL, symbol='pL')
    pq.fL = pq.UnitQuantity('femtoliter', 0.001 * pq.pL, symbol='fL')
    pq.uM = pq.UnitQuantity('micromolar', 0.001 * pq.mM, symbol='uM')
    pq.nM = pq.UnitQuantity('nanomolar', 0.001 * pq.uM, symbol='nM')
    pq.Rgas = pq.UnitQuantity('gas constant', 1.987 * pq.cal / pq.Kelvin / pq.mol )
    pq.kb = pq.UnitQuantity('Boltzmann\'s constant', 1.38E-23 * pq.J / pq.K)
    pq.Cal = pq.UnitQuantity('Calorie', 1000*pq.cal )
    pq.plank = pq.UnitQuantity('Plank\'s constant', 6.626E-34 * pq.J * pq.s )
    pq.molecule = pq.UnitQuantity('molecule', pq.mole/6.0221409e+23, symbol='molecule')
    return pq

class Organelle():
    def __init__(self, name):
        self.name = name

class Yeast():

    def __init__(self):
        self.cell = Organelle("cell")
        self.nucleus = Organelle('nucleus')
        self.cell.diameter = pq.Quantity(5, 'um')
        self.cell.volume = pq.Quantity(50, 'um**3')
        self.nucleus.diameter = pq.Quantity(2, 'um')
        self.nucleus.volume = pq.UncertainQuantity(2.9, 'um**3', 0.2)
        pq.nucleur_volume = pq.UnitQuantity('nuclear volume', self.nucleus.volume, symbol='NucVol')
        pq.cell_volume = pq.UnitQuantity('cellular volume', self.cell.volume, symbol='CellVol')

class Cas9():
    def __init__(self):
        self.Kd = pq.UncertainQuantity(1.2, pq.nM, 0.1)# From literature
        self.koff = pq.Quantity(5.5e-5, 1/pq.s).rescale(1/pq.min)
        self.kon = pq.Quantity(4.6e4, 1/(pq.s*pq.M))
        self.kon = self.kon.rescale(1/(pq.min*pq.nM))