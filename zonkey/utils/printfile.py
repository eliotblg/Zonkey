
BOHR2ANG = 0.5291772109217
ANG2BOHR = 1.0 / BOHR2ANG

def printxyz(coords, atypes, filename, append=False):
    if filename[-4:] != '.xyz':
        filename = filename + '.xyz'
    if append:
        fp = open(filename, 'a')
    else:
        fp = open(filename, 'w')
    fp.write(str(len(coords)) + '\n')

    for i, xyz in enumerate(coords):
        fp.write('\n' + atypes[i] + ' ' + \
                        str(xyz[0]*BOHR2ANG) + ' ' + \
                        str(xyz[1]*BOHR2ANG) + ' ' + \
                        str(xyz[2]*BOHR2ANG) + ' ')
    fp.write('\n')
    fp.close()

