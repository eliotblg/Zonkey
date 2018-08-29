import os
import imp
import numpy as np
from chemsys import Chemsys

listofqmcodes = ['gaussian']

#
# Generic interface to QM codes for QM/MM or pure QM computations 
# No code specific commands/data is implied here
# It loads and uses the selected qm module from the interface directory
#
class QMinterface(object):

    def __init__(self, qmcode, epath=None, name="qmjob", memory=1, nproc=1, \
                 spath = None):
        self.execpath = epath
        if epath == None:
            self.execpath = None #TODO Here load from files and check
        self.scratchpath = spath
        self.name = name
       
        self.memory=memory
        self.nproc=nproc

        # load QM module
        srcdir = os.path.dirname(os.path.abspath(__file__))
        self.qmmodule = imp.load_source(qmcode, \
                        srcdir + '/interfaces/' + qmcode + '.py')
        # Check if environement for the QM module is set
        if not self.qmmodule.checkenvironment():
            self.qmmodule.setupenvironment(epath, spath)

    def setmethod(self, hamiltonian, basis=None, extra=''):
        self.hamiltonian = hamiltonian
        self.basis = basis
        self.extra = extra

    def runjob(self, coords, jfile='qmjob', method='sp', memory=1, nproc=1):
        #TODO read self.memory and self.nproc if None is provided
        self.qmmodule.printjob(coords, jfile=jfile, nproc=nproc, mem=memory, \
                               ham=self.hamiltonian, basis=self.basis,       \
                               meth=method, charge=coords.charge,            \
                               mult=coords.mult, extra=self.extra)
        self.qmmodule.runjob(self.execpath, jfile=jfile)
        ##self.qmmodule.checkjob()

    def energy(self, coords):
        self.runjob(coords, jfile='energy-' + self.name, method='sp', \
                    memory=self.memory, nproc=self.nproc)
        ener = self.qmmodule.extractdata(coords, jfile='energy-' + self.name, \
                                         val='energy')
        coords.qmenergy = ener
        return ener

    # compute energy and gradients
    def gradients(self, coords):
        self.runjob(coords, jfile='gradients-' + self.name, method='grad', \
                    memory=self.memory, nproc=self.nproc)    
        ener, grad = self.qmmodule.extractdata(coords, jfile='gradients-' + \
                                               self.name, val='gradients')
        coords.qmenergy = ener
        return ener, grad

