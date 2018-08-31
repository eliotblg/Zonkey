
ANG2BOHR        = 0.5291772109217
BOHR2ANG        = 1.0/ANG2BOHR

def printxyz(coords, atypes, filename):
    if filename[-4:] != '.xyz':
        filename = filename + '.xyz'
    fp = open(filename, 'w')
    fp.write(str(len(coords)) + '\n')
    for i, xyz in enumerate(coords):
        fp.write('\n' + atypes[i] + ' ' + \
                        str(xyz[0]*BOHR2ANG) + ' ' + \
                        str(xyz[1]*BOHR2ANG) + ' ' + \
                        str(xyz[2]*BOHR2ANG) + ' ')
    fp.close()

