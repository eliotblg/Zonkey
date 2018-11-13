import os
import imp
import time
import numpy as np
from . chemsys import Chemsys


BOHR2ANG = 0.5291772109217
ANG2BOHR = 1.0 / BOHR2ANG

HARTREE2KCALMOL = 627.50947415
KCALMOL2HARTREE = 1.0 / HARTREE2KCALMOL

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

    def interactive(self, system, printener=True):
        f = open("INTERACTIVE", "w")
        f.write('START')
        f.close()
        mmijob = self.mm.interactive(system)
        while True:
            #TODO has to kill itself once waited for very too long
            time.sleep(0.1)
            f = open("INTERACTIVE", "r")
            w = f.readline()
            f.close()
            if w[0:2] == "QM":
                # update system coordinates for QM computations
                # intcoords.txt file is printed by the interactive MM module 
                # TODO MM module has to handle the convertion factor 
                with open("intcoords.txt") as fp:
                    try:
                        for i, line in enumerate(fp):
                            l = line.rstrip().split()
                            system.coords[i] = np.array([float(l[0])*ANG2BOHR, \
                                    float(l[1])*ANG2BOHR, float(l[2])*ANG2BOHR])
                        print(system.coords)
                    except:
                        print('Cannot read coordinates from MM interactive module.\n')
                self.eqm, gqm = self.qm.gradients(system)
                f = open("qmgradients.txt", "w")
                # convert gradients to forces in NAMD units (to change if other code used)
                conv = -1.0 * HARTREE2KCALMOL * BOHR2ANG
                for g in gqm:
                    f.write(str(g[0]*conv) + ' ' + str(g[1]*conv) + \
                            ' ' + str(g[2]*conv) + '\n')
                f.close()
                f = open("INTERACTIVE", "w")
                f.write('MM')
                f.close()
            elif len(w) > 2 and w[0:4] == "STOP":
                # kill the job if still running | it's a subprocess object
                if mmijob.poll() == None:
                    mmijob.kill()
                break

    def clean(self):
        for f in ["INTERACTIVE", "qmgradients.txt"]:
            if os.path.isfile(f):
                os.remove(f)
        self.qm.clean()
        self.mm.clean()

    def print__(self):
        print("QM/MM energies (Hartrees)")
        print("%12s = %12s + %12s")%('Eqm/mm', 'Eqm', 'Emm')
        print("%.12E   %.12E   %.12E")%(self.eqmmm,self.eqm,self.emm)

