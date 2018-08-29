import os
import re
import numpy as np
from chemsys import Chemsys
##import constant
# TODO use restart file with checkpoints

ANG2BOHR          = 0.5291772109217
BOHR2ANG          = 1.0 / ANG2BOHR

#MY_ANG2BOHR       = 0.5291772109217
#GAUSSIAN_ANG2BOHR = 0.52917721092
#GAUSSANG2MYANG = MY_ANG2BOHR / GAUSSIAN_ANG2BOHR

# Check if environment is set
def checkenvironment():
    if not 'GAUSS_EXEDIR' in os.environ:
        print('GAUSS_EXEDIR not set')
        return False
    if not 'GAUSS_SCRDIR' in os.environ:
        print('GAUSS_SCRDIR not set')
        return False
    return True

# Set environment; both variables are mandatory
def setupenvironment(executable=None, scratch=None):
    if executable != None:
        gaussiandir = os.path.dirname(executable)
        os.environ['GAUSS_EXEDIR'] = gaussiandir
    elif not 'GAUSS_EXEDIR' in os.environ:
        print('You have to provide a gaussian executable directory path')
        exit(1)
    if scratch != None:
        os.environ['GAUSS_SCRDIR'] = scratch 
    elif not 'GAUSS_SCRDIR' in os.environ:
        print('You have to provide a gaussian scratch directory path')
        exit(1)

# Print job file
def printjob(coords, jfile='jf', nproc=1, mem=1, ham='hf', basis='STO-3G', \
             meth='sp', charge=None, mult= None, extra=''):
    # if no charge or multiplicity is specified, load them from the chemsys
    if charge == None:
        charge = coords.charge
    if mult == None:
        mult = coords.mult
    method = meth
    hamiltonian = ham

    # open input file
    fp = open(jfile + '.com', 'w')
    # do some translation specific to gaussian
    if meth == 'sp':
        method = ''
    elif meth == 'grad':
        method = 'force'
    if ham == 'pbe0':
        hamiltonian = 'pbe1pbe'
    if basis != None and basis != '':
        hamiltonian = hamiltonian + '/' + basis
   
    # if QM/MM load extra keywords to include point charges
    if coords.nmmatoms > 0:
        extra = extra + ' Charge'
        if "nosymm" not in extra.lower():
            extra = extra + ' NoSymm'
        if meth == 'grad':
            extra = extra + ' Prop=Grid' #'SCF=verytight'

    if os.path.isfile(jfile + '.chk'):
        extra += ' Guess=Read'

    # write the actual job file
    # write the header
    fp.write('%nproc=' + str(nproc) + '\n')
    fp.write('%mem=' + str(int(round(mem))) + 'GB\n')
    fp.write('%Chk=' + jfile + '.chk\n')
    fp.write('# ' + method + ' ' + hamiltonian + ' ' + extra + '\n\n')
    fp.write('This is a comment\n\n')
    fp.write(str(charge) + ' ' + str(mult) + '\n')

    # if no MM atoms just run std QM computation
    #if coords.nmmatoms == 0:
    #    coords.nqmatoms = coords.natoms
    #    coords.qmatoms  = range(coords.natoms) 
    # if MM atoms are set but no QM atoms => error
    if coords.nqmatoms == 0:
        print('No QM atoms defined in QM computation using Gaussian')
        exit(1)
    # QM/MM 
    else:
        for i in coords.qmatoms:
            c = coords.coords[i]
            fp.write(coords.atypes[i] + ' '+ str(c[0]*BOHR2ANG) + ' ' + \
                     str(c[1]*BOHR2ANG) + ' ' + str(c[2]*BOHR2ANG) + '\n')
        fp.write('\n')
        # write MM cordinates and charge values to be included into QM computation
        for i in coords.mmatoms:
            c = coords.coords[i]
            fp.write(str(c[0]*BOHR2ANG) + ' ' + str(c[1]*BOHR2ANG) + ' ' + \
                     str(c[2]*BOHR2ANG) + ' ' + str(coords.charges[i]) + '\n')

        # add grid point at MM atoms to evaluate the potential
        if meth == 'grad' and coords.nmmatoms > 0:
            fp.write('\n' + str(coords.nmmatoms) + ',2,41,42\n')
            fortf = open('fort.41', 'w')
            for i in coords.mmatoms:
                c = coords.coords[i]
                fortf.write('%20f%20f%20f\n'%(c[0]*BOHR2ANG, c[1]*BOHR2ANG, c[2]*BOHR2ANG))
            fortf.close()
                
    # mandatory blank line at the end 
    fp.write('\n')
    fp.close()

def runjob(executable, jfile='jf'):
    os.system(executable + ' < ' + jfile + '.com > ' + jfile + '.log')

#TODO def checkjob()

def extractdata(coords, jfile='jf', val='energy'):
    with open(jfile + '.log', 'r') as f:
        fdata = f.read()
    if val in ['energy', 'gradients']:
        # look after the energy
        energy = None
        reg = '.*SCF Done:.*'
        match = re.search(reg,fdata)
        if match:
            energy = float(match.group().split()[4])
            # if just energy needed, return here
            if val == 'energy':
                return energy 
        if energy == None:
            print('Could not locate energy in ' + jfile + '.log')
            exit(1)
        # look after the QM gradients
        grad = np.zeros([coords.natoms,3])
        ngrad = -2 # skip 2 lines after finding indication of gradients 
        storegrad = False
        for line in fdata.split('\n'):
            if storegrad:
                d = line.split()
                if ngrad >= coords.nqmatoms:
                    break
                if ngrad >= 0:
                    grad[coords.qmatoms[ngrad]] = [float(d[2]),float(d[3]),float(d[4])]         
                ngrad += 1
            if line[0:7] == ' Center':
                d = line.split()
                if d[2] == 'Forces':
                    storegrad = True
        # look for the MM gradients
        if coords.nmmatoms:
            fortf = open('fort.42','r')
            for i in coords.mmatoms:
                fortf.readline()
                d = fortf.readline().split()
                q = coords.charges[i]
                grad[i] = [float(d[0])*q, float(d[1])*q, float(d[2])*q]                
    grad = np.array(grad) * -1.0
    return energy, grad

