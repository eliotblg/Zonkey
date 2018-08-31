import os
import imp
import numpy as np
from . chemsys import Chemsys

class QMMMinterface(object):

    def __init__(self, qmmodule, mmmodule, name="qmmm", memory=1, nproc=1):
        self.qm = qmmodule
        self.mm = mmmodule
        self.emm = 0.0
        self.eqm = 0.0
        self.eqmmm = 0.0

    def energy(self, system, printener=True):
        self.emm = self.mm.energy(system)
        self.eqm = self.qm.energy(system)
        self.eqmmm = self.emm + self.eqm
        system.qmmenergy = self.eqmmm 
        if printener:
            system.printenergies()
        return self.eqmmm

    def gradients(self, system, printener=True):
        self.emm, gmm = self.mm.gradients(system)
        self.eqm, gqm = self.qm.gradients(system) 
        self.eqmmm = self.emm + self.eqm
        gqmmm = gqm + gmm
        system.qmmmenergy = self.eqmmm
        if printener:
            system.printenergies()
        return self.eqmmm, gqmmm

    def clean(self):
        self.qm.clean()
        self.mm.clean()

    def print__(self):
        print("QM/MM energies (Hartrees)")
        print("%12s = %12s + %12s")%('Eqm/mm', 'Eqm', 'Emm')
        print("%.12E   %.12E   %.12E")%(self.eqmmm,self.eqm,self.emm)

