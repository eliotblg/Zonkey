import os
import imp
import numpy as np
from chemsys import Chemsys

class QMMMinterface(object):

    def __init__(self, qmmodule, mmmodule, name="qmmm", memory=1, nproc=1):
        self.qm = qmmodule
        self.mm = mmmodule
        self.emm = 0.0
        self.eqm = 0.0
        self.eqmmm = 0.0

    def energy(self, coords, printener=True):
        self.emm = self.mm.energy(coords)
        self.eqm = self.qm.energy(coords)
        self.eqmmm = self.emm + self.eqm
        coords.qmmenergy = self.eqmmm 
        if printener:
            coords.printenergies()
        return self.eqmmm

    def gradients(self, coords, printener=True):
        self.emm, gmm = self.mm.gradients(coords)
        self.eqm, gqm = self.qm.gradients(coords) 
        self.eqmmm = self.emm + self.eqm
        gqmmm = gqm + gmm
        coords.qmmmenergy = self.eqmmm
        if printener:
            coords.printenergies()
        return self.eqmmm, gqmmm

    def print__(self):
        print("QM/MM energies (Hartrees)")
        print("%12s = %12s + %12s")%('Eqm/mm', 'Eqm', 'Emm')
        print("%.12E   %.12E   %.12E")%(self.eqmmm,self.eqm,self.emm)

