import os
from shutil import copyfile
import numpy as np

from . utils import *

HARTREE2KCALMOL = 627.50947415
KCALMOL2HARTREE = 1.0 / HARTREE2KCALMOL

#
# Chemical System class
# Contains everything about atoms, molecules, etc...
# It is the main class of zonkey
#
class Chemsys(object):
    def __init__(self, coordfile=None):
        self.charge = 0
        self.mult   = 1
        if coordfile is None:
            self.natoms = 0
            self.nbonds = 0
            self.infile = None
        else:
            zkfile = 'zk-' + coordfile
            copyfile(coordfile, zkfile)
            self.infile = coordfile
            self.zkfile = zkfile
            self.readcoordfile(coordfile)
        self.qmenergy  = None
        self.mmenergy  = None
        self.qmmmenergy = None

    def readcoordfile(self, coordfile):
        self.infile = coordfile
        
        # check type of file based on the extension
        fileext = os.path.splitext(coordfile)[-1][1:]
        print('Reading coordinates from ' + coordfile)
        if fileext == 'xyz':
            self.atypes, self.coords = readfile.readxyz(coordfile)
        elif fileext == 'pdb':
            self.atypes, self.coords = readfile.readpdb(coordfile)
        else:
            print('file extension ' + fileext + ' not supported')
            exit(1)

        self.natoms = len(self.coords)
        self.charges = np.zeros(self.natoms)

        if self.natoms == 0:
            print('No atom found in coord file')
            exit(1)

        self.coords = np.array(self.coords)
        self.grad   = np.zeros([len(self.coords),3])

    def getbonds(self, structure=None):
        # implementation only covers psf at this stage
        if structure != None:
#            srcdir = os.path.dirname(os.path.abspath(__file__))
#            mpsf = imp.load_source('mpsf', srcdir + '/utils/manipulatepsf.py')    
            self.bonds = manipulatepsf.getbonds(structure)
            self.nbonds = len(self.bonds)
        else:
            print('No other method than reading from psf exist yet for getting bonds')
            exit(1)
        
    # set QM and MM region
    def setregions(self, qmlist = None, mmlist = None):
        if self.natoms == 0:
            print('You should load coordinates before setting regions')
            exit(1)
        self.qmatoms = []
        self.mmatoms = []
        self.ghostatoms = []
        if qmlist != None:
            self.qmatoms = qmlist
        else:
            ##print(mmlist)
            if mmlist == None or len(mmlist) == 0:
                print('Error setting QM and MM regions: should give qmlist or mmlist')
                exit(1)
            for i in range(self.natoms):
                if i not in mmlist:
                    self.qmatoms.append(i)

        if mmlist != None and len(mmlist) > 0:
            self.mmatoms = mmlist
        else:
            for i in range(self.natoms):
                if i not in qmlist:
                    self.mmatoms.append(i)

        for i in range(self.natoms):
            if i not in self.qmatoms and i not in self.mmatoms:
                 self.ghostatoms.append(i)
                
        self.nqmatoms = len(self.qmatoms)
        self.nmmatoms = len(self.mmatoms) 
        self.nghostatoms = len(self.ghostatoms)

        print(str(self.nqmatoms) + ' atoms in QM region\n' + \
              str(self.nmmatoms) + ' atoms in MM region\n' + \
              str(self.nghostatoms) + ' ghost atoms')

    def setactivefrozen(self, alist=None, flist=None):
        if self.natoms == 0:
            print('You should load coordinates before setting active/frozen regions')
            exit(1)

        if alist != None:
            self.activelist = alist
            if flist != None:
                self.frozenlist = flist 
                if set(alist) & set(flist):
                    print('Overlaping active and frozen lists')
                    exit(1)
                else:
                    self.frozenlist = []
                    for i in range(self.natoms):
                        if i not in alist:
                            self.frozenlist.append(i)
        elif flist != None: 
            self.frozenlist = flist
            self.activelist = []
            for i in range(self.natoms):
                if i not in flist:
                    self.activelist.append(i)
        else:
            return 0
        self.nactives = len(self.activelist)
        self.nfrozens = len(self.frozenlist)
        print(str(self.nactives) + ' atoms in active region\n' + \
              str(self.nfrozens) + ' atoms in frozen region\n' + \
              str(self.natoms - self.nactives - self.nfrozens) + \
              ' atoms have no active/frozen region\n') 


    def printenergies(self, units='KCAL/MOL'):
        if units == 'KCAL/MOL':
            scale = HARTREE2KCALMOL
        else:
            scale = 1.0
            units = 'HARTREE'
        if self.qmenergy is not None:
            print('QM    energy (' + units.lower() + '): ' + \
                  str(self.qmenergy*scale))
        if self.mmenergy is not None:
            print('MM    energy (' + units.lower() + '): ' + \
                  str(self.mmenergy*scale))
        if self.qmmmenergy is not None:
            print('QM/MM energy (' + units.lower() + '): ' + \
                  str(self.qmmmenergy*scale))







