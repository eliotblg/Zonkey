import sys
import numpy as np
from readfile import readpdb

def selectaround(coords, r=np.array([0.0, 0.0, 0.0]), dist=20.0, \
                 pbc=None, pbcenter=None):
    dist2 = dist*dist
    print(coords)

names, coords = readpdb(str(sys.argv[1]))
selectaround(coords)
 
