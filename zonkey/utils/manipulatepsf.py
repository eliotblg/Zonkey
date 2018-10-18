#
# Set the charges of a list of atoms to the input value 
#
def setcharges(psf, atomlist, charge=0.0, psfout='mod.psf'):
    # convert charge to proper format for printing in the psf
    scharge = str(charge)
    if scharge[0] != '-':
        scharge = ' ' + scharge
    if len(scharge) > 9:
        scharge = scharge[0:9]
    if len(scharge) < 9:
        for i in range(len(scharge), 9):
            scharge = scharge + '0'

    # new psf file with modified charges
    npsf = open(psfout, 'w')

    check = False 
    acount = 0
    with open(psf) as fp:
        for line in fp:
            ll = line.split()
            if check == True and int(ll[0]) - 1 in atomlist:
                line = line[0:55] + scharge + line[64:]
                acount += 1
            if len(ll) and ll[1] == '!NATOM':
                check = True 
            if acount == len(atomlist):
                check = False
            npsf.write(line)

    npsf.close()
    return psfout

#
# Get charges from psf and output them as a list of natoms size
#
def getcharges(psf):
    charges = []
    check = False
    acount = 0
    natoms = -1
    with open(psf) as fp:
        for line in fp:
            ll = line.split()
            if check == True:
                charges.append(float(line[55:70]))
                acount += 1
            if len(ll) and ll[1] == '!NATOM':
                check = True
                natoms = int(ll[0])
            if acount == natoms:
                break
    return charges

#
# returns a list of bonds taken directly from the psf
#
def getbonds(psf):
    bondlist = []
    check = False
    with open(psf) as fp:
        for line in fp:
            ll = line.rstrip().split()
            if len(ll) >=2 :
                if ll[1] == '!NTHETA:':
                    break
                if check:
                    for i in range(0, len(ll), 2):
                        bondlist.append([int(ll[i])-1, int(ll[i+1])-1])     
                if ll[1] == '!NBOND:':
                    check = True 
    return bondlist

#
# Remove bonded terms from PSF. Take a list of atoms as input. 
# The terms removed are those listed in the nabond dictionary 
# which maps them to the number of atoms used to describe them
#
# !!! VDW or electrostatic 1-4 interaction definition
#
def removebonded(psf, atomlist , psfout='modbond.psf'):
    npsf = open(psfout, 'w')
    check = 0 
    nabond = {'!NBOND:':2, '!NTHETA:':3, '!NPHI:':4, '!NIMPHI:':4, '!NDON:':2, '!NACC:':2}

    if atomlist != None:
        satomlist = [str(i + 1) for i in atomlist]
    else:
        print('ERROR in removebonded in manipulate psf: no atomlist provided')
        exit(1)
    with open(psf) as fp:
        for line in fp:
            ll = line.rstrip().split()
            if check:
                if len(ll):
                    for i in range(0, len(ll), check):
                        a = [ll[j] for j in range(i,i+check)]
                        if set(a).isdisjoint(satomlist):                    
                            al.append(a)
                            nc += 1
                        nb += 1
            if len(ll) and ll[1] in nabond: 
                check = nabond[ll[1]]
                savedline = line
                na = int(ll[0])
                nb = 0
                nc = 0
                al = []
            if not check:
                npsf.write(line)
            elif nb >= na:
                tmp = '%10d'%nc + savedline[10:]
                j = 0
                while j < nc:
                    for i in range(0, int((min(8 + (check%2), (nc-j)*check))/check)):
                        for k in al[j]:
                            tmp = tmp + '%10s'%k
                        j += 1 
                    tmp = tmp + '\n'
                npsf.write(tmp)
                check = 0
            
    npsf.close()
    return psfout

#
# create new atom types for QM atoms
# they interact with MM atoms through vdw interactions at MM level
# they don't interact together with these interactions
#
def modatomtype(psf, prm, atomlist, psfout='newtypes.psf', \
                prmout='newtypes.prm', startnum=1, createnbfix=True):
    check = False
    acount = 0
    atommap = []
    npsf = open(psfout, 'w')

    with open(psf) as fp:
        for line in fp:
            ll = line.split()
            if check == True and int(ll[0]) - 1 in atomlist:
                oldname = line[47:55]
                newname = "{:<8}".format('QM' + str(startnum + acount))
                line = line[0:47] + newname +line[55:]
                atommap.append([oldname.strip(), newname.strip()]) 
                acount += 1

            if len(ll) and ll[1] == '!NATOM':
                check = True
            if acount == len(atomlist):
                check = False
            npsf.write(line)

    #atomparam = [[0 for x in range(4)] for y in range(len(atomlist))]
    if type(prm) is not list: prm = [ prm ]
    for prmfile in prm:
        check = False
        with open(prmfile) as fp:
            for line in fp:
                if line[0:4] == 'MASS':
                    ll = line.split()
                    for j in range(len(atommap)):
                        if atommap[j][0] == ll[2]:
                            atommap[j].append(ll[3])
                if check:
                    ll = line.split()
                    for j in range(len(atommap)):
                        if len(ll) and atommap[j][0] == ll[0]:
                            atommap[j].append(ll[2])
                            atommap[j].append(ll[3])
                if line[0:4] == 'NONB':
                    check = True    
    nprm = open(prmout, 'w')
    wprm = '* New atom types to use in QM/MM computation with Zonkey *\n*\n\nATOMS\n'
    for a in atommap:
        wprm = wprm + 'MASS    -1 ' + a[1] + '  ' + a[2] + '\n'
    wprm = wprm + '\nBONDS\n\nANGLES\n\nIMPROPERS\n'
    wprm = wprm + '\nNONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -\n'
    wprm = wprm + 'cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5\n\n'
    for a in atommap:
        wprm = wprm + a[1] + '   0.0   ' + a[3] + '  ' + a[4] + '\n'
    if createnbfix:
        wprm += '\nNBFIX\n'
        for i in range(len(atommap)+1):
            for j in range(i+1,len(atommap)):
                wprm += atommap[i][1] + ' ' + atommap[j][1] + '   0.0    0.0\n'


    wprm = wprm + '\nEND\n'
          
    nprm.write(wprm)
    npsf.close()
    nprm.close()

    return psfout, prmout

