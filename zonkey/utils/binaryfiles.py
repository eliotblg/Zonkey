import os
import struct

ANG2BOHR        = 0.5291772109217
BOHR2ANG        = 1.0 / ANG2BOHR

def readnamdbin(filename):
    with open(filename, 'rb') as fp: 
        natoms=struct.unpack('i', fp.read(4))[0]
        x = list(struct.unpack('d'*natoms*3, fp.read(8*natoms*3)))
    return [x[i:i+3]*ANG2BOHR for i in range(0, natoms*3, 3)] 

def printnamdbin(filename, x):
    with open(filename, 'wb') as fp:
        fp.write(struct.pack('i', len(x)))
        for xx in x:
            fp.write(struct.pack('ddd', \
            xx[0]*BOHR2ANG, xx[1]*BOHR2ANG, xx[2]*BOHR2ANG))

# TODO Have some "block" reading if the file size > available memory
def readdcd(filename, maxfsize=2000000000):
    ##if os.path.getsize(filename) < maxfsize:
    with open(filename, 'rb') as fp:
        natoms = struct.unpack('cccc', fp.read(4))[0]
        x = list(struct.unpack('d'*natoms*3, fp.read(8*natoms*3)))
    return [x[i:i+3]*ANG2BOHR for i in range(0, natoms*3, 3)]


#
