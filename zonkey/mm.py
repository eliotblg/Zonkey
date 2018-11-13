import os
import imp
import numpy as np
from . chemsys import Chemsys

listofmmcodes = ['namd']

#
# Generic interface to MM codes for QM/MM or pure MM computations 
# No code specific commands/data is implied here
# It loads and uses the selected mm module from the interfaces directory
#
class MMinterface(object):

    def __init__(self, mmcode, epath=None, name="mmjob", memory=1, nproc=1, \
                 spath = None):

        self.execpath = epath
        if epath == None:
            self.execpath = None #TODO Here load from files and check
        self.scratchpath = spath
        self.name = name
        self.memory=memory
        self.nproc=nproc

        #Load MM module | TODO do this part properly
        srcdir = os.path.dirname(os.path.abspath(__file__))
        self.mmmodule = imp.load_source(mmcode, srcdir + '/interfaces/' \
                                      + mmcode + '.py')
        #Check if environement for the MM module is set
        if not self.mmmodule.checkenvironment():
            self.mmmodule.setupenvironment(epath, spath)

    def setstructure(self, coords, refcoords, structure, prm):
        if type(prm) is not list: prm = [ prm ]
        self.ref = refcoords
        self.struct = structure
        self.prm = prm 
        if coords.nqmatoms > 0:
            self.struct, self.prm = self.mmmodule.qmmmprepare(coords, \
                                    structure, prm)

    def runjob(self, coords, jfile='mmjob', method='sp', memory=1, nproc=1, \
               pbc = None, cutoff = None, extra = None, nsteps = None, \
               temperature = None):
        self.mmmodule.printjob(coords, jfile=jfile, nproc=nproc, mem=memory, \
                               structure = self.struct, prm = self.prm, \
                               refcoords = self.ref, meth=method, pbc = pbc, \
                               cutoff = cutoff, extra=extra, nsteps=nsteps, \
                               temperature = temperature)
        if method == 'interactive':
            mmijob = self.mmmodule.runinteractivejob(coords, self.execpath, \
                                                     jfile=jfile, nproc=nproc)
            return mmijob
        else:
            self.mmmodule.runjob(coords, self.execpath, jfile=jfile, nproc=nproc)

    def energy(self, coords, pbc = None, cutoff = None, extra = None):
        jobname = 'mmener-' + self.name
        self.runjob(coords, jfile = jobname, method = 'sp', \
                    memory = self.memory, nproc = self.nproc, pbc = pbc, \
                    cutoff = cutoff, extra = extra)  
        ener = self.mmmodule.extractdata(coords, jfile=jobname,val='energy')
        coords.mmenergy = ener
        return ener

    def gradients(self, coords, pbc = None, cutoff = None, extra = None):
        jobname = 'mmgrad-' + self.name
        self.runjob(coords, jfile = jobname, method = 'grad', \
                    memory = self.memory, nproc = self.nproc, pbc = pbc, \
                    cutoff = cutoff, extra = extra)
        e, g = self.mmmodule.extractdata(coords, jfile=jobname, val='gradient')
        coords.mmenergy = e
        return e, g

    def interactive(self, coords, pbc = None, cutoff = None, extra = None, \
                    nsteps = 1000, temperature = 30):
        jobname = 'mminter-' + self.name
        mmijob = self.runjob(coords, jfile = jobname, method = 'interactive', \
                    memory = self.memory, nproc = self.nproc, pbc = pbc, \
                    cutoff = cutoff, extra = extra, nsteps=nsteps, \
                    temperature = temperature)
        return mmijob 

    def clean(self):
#        self.mmmodule.clean(['mmgrad-' + self.name, 'mmener-' + self.name, \
#                             'mminter-' + self.name], self.struct)
        pass







