from zonkey import *
from optparse import OptionParser

# parse and extract options from input file
parser = OptionParser()
parser.add_option('-i', '--input', dest='infile', \
                  help='Input file', metavar='FILE')

(options, args) = parser.parse_args()
if not options.infile:   # if filename is not given
    inout.printoptionhelp()
    parser.error( \
        'No input file specified\nusage: \'python zonkey.py -i myinputfile\'')

# read input file and obtain a dictionary of options
opts = inout.readinputfile(options.infile)

###############################

#TODO have coordinates potentialy read from a .xyz or alike file

# create a "Chemical system" based on the input coordinates
system = Chemsys(opts['coordinates'])

# set QM and MM region
system.setregions(qmlist=opts['qmlist'], mmlist=opts['mmlist'])

# set active and frozen region
system.setactivefrozen(alist=opts['activelist'], flist=opts['frozenlist'])

# get charges from psf file | here can input charge manualy too
system.getcharges(structure = opts['mmstruct'])

# create an interface to the given QM program
if system.nqmatoms > 0:
    qmjob = QMinterface(opts['qmcode'],epath='/home/eliot/software/g09E01/g09',\
                     memory=opts['memory'], nproc=opts['nproc'],spath='/tmp')
    qmjob.setmethod(opts['qmham'], basis=opts['qmbasis'], extra=opts['qmextra'])
    if system.nmmatoms <= 0:
        method = qmjob
# create MM
if system.nmmatoms > 0 :
    mmjob = MMinterface(opts['mmcode'], \
        epath='/home/eliot/software/NAMD_2.12_Linux-x86_64-multicore/namd2', \
        memory = opts['memory'], nproc = opts['nproc'])
    mmjob.setstructure(system, system.infile, structure = opts['mmstruct'], \
        prm = opts['mmparam'])
    if system.nqmatoms <= 0:
        method = mmjob

# if both QM and MM atoms, load QM/MM module
if system.nqmatoms > 0 and system.nmmatoms > 0:
    method = QMMMinterface(qmjob, mmjob)


############################################################################

optimizer = Geopt(method)
#optimizer.setoptimizer(system.nactives)
#optimizer.optimize(system)
optimizer.steepestd(system)

#method.clean()

