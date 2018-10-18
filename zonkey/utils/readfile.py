import os
import imp
import numpy as np

BOHR2ANG = 0.5291772109217
ANG2BOHR = 1.0 / BOHR2ANG

def readxyz(coordfile):
    anames = []
    coord = []
    with open(coordfile) as fp:
        coorf = fp.readlines()
    natoms = int(coorf[0].rstrip())
    for i in range(2, natoms+2):
        c = coorf[i].rstrip().split()
        anames.append(c[0])
        coord.append(np.array([ float(c[1])*ANG2BOHR, \
                                float(c[2])*ANG2BOHR, \
                                float(c[3])*ANG2BOHR ]))
    return anames, np.array(coord)


def readpdb(coordfile):
    anames = []
    coord  = []
    with open(coordfile) as fp:
        for line in fp:
            l = line.split()
            if len(l) and l[0] in ['ATOM', 'HETAM']:
                anames.append(getelement(l[2])) 
                m = line[28:56].split()
                coord.append(np.array([ float(m[0])*ANG2BOHR, \
                                        float(m[1])*ANG2BOHR, \
                                        float(m[2])*ANG2BOHR ]))
    return anames, np.array(coord) 

def getelement(aname):
    # not always obvious with force field definition
    digitaname2 = ['LI', 'BE', 'MG', 'AL', 'SI', 'CL', 'AR', 'FE']
    if len(aname) > 1:
        if aname[0:2] in digitaname2:
            return aname[0:2]
            if aname[0:3] == 'POT':
                return 'K' 
            if aname[0:3] == 'CLA':
                return 'CL'
    return aname[0]



