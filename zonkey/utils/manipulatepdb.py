import os
import imp
import shutil 
import numpy as np

BOHR2ANG = 0.5291772109217
ANG2BOHR = 1.0 / BOHR2ANG


def changecoordinates(coordfile, atoms, newcoords, newfile = None):
#    tomodify = dict(zip(atoms, [["% 3.3f"%(j*BOHR2ANG) for j in i] for i in newcoords]))
    tomodify = dict(zip(atoms, [['{0:> 8.3f}'.format(j*BOHR2ANG) \
                                 for j in i] for i in newcoords]))

    # create tmp pdb to avoid loosing data if mistake
    upcoords = open('tmp.pdb','w')
    with open(coordfile, 'r') as fp:
        n = 0
        for line in fp.readlines():
            if line[0:4] == 'ATOM':
                if n in tomodify:
                     line = line[0:30] + tomodify[n][0] + tomodify[n][1] + \
                           tomodify[n][2] + line[54:]
                n+=1
            upcoords.write(line)

    upcoords.close()
    if newfile == None:
        shutil.copy("tmp.pdb", coordfile)
    else:
        shutil.copy("tmp.pdb", newfile)
    os.remove("tmp.pdb")
