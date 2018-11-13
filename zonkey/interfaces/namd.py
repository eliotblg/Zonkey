import os
import re
import imp
import numpy as np
import subprocess

srcdir = os.path.dirname(os.path.abspath(__file__))
manipsf = imp.load_source('manipulatepsf', srcdir + \
                          '/../utils/manipulatepsf.py')
manipdb = imp.load_source('manipulatepdb', srcdir + \
                          '/../utils/manipulatepdb.py')
binfile = imp.load_source('binaryfiles', srcdir + \
                          '/../utils/binaryfiles.py')

BOHR2ANG = 0.5291772109217
ANG2BOHR = 1.0 / BOHR2ANG
HARTREE2KCALMOL = 627.50947415
KCALMOL2HARTREE = 1.0 / HARTREE2KCALMOL 

# dictionary of default namd keywords used
defnkey = {'structure': None, 'coordinates': None, 'outputname': None,       \
           'temperature': '0.0','exclude': 'scaled1-4', '1-4scaling': '1.0', \
           'cutoff': '999.9','switchdist': '997.9','pairlistdist': '1001.0', \
           'timestep': '1.0','nonbondedFreq': '1','fullElectFrequency':'1',  \
           'stepspercycle': '1', 'paraTypeCharmm': 'on'}

# dictionary of default namd keywords used for intractive method
defnkeyint = defnkey
defnkeyint.update({ 'temperature': '300', 'cutoff': '12.0', 'switchdist': '10.0', \
           'pairlistdist': '14.0', 'restartfreq': '100', 'dcdfreq': '1', 'xstfreq': '1', \
           'outputEnergies': '1', 'switching': 'on', 'vdwForceSwitching': 'yes', \
           'rigidBonds':'none' })

# dictionary of defaults for frozen atoms 
keyfrozen = {'fixedAtoms': 'on', 'fixedAtomsForces': 'off', 'fixedAtomsFile': 'fixed.pdb'}

# dictionary of defaults for periodic boundary conditions
keyperiod = {'wrapAll': 'on', 'PME': 'yes', 'PMEInterpOrder': 6.0, 'PMEGridSpacint': 1.0, \
             'useGroupPressure': 'yes', 'langevin': 'on', 'langevinDamping': 'on', \
             'langevinTemp': '300', 'langevinHydrogen': 'off', 'langevinPiston': 'on', \
             'langevinPistonTarget': '1.01325', 'langevinPistonPeriod': '50.0', \
             'langevinPistonDecay': '25.0', 'langevinPistonTemp': '300'}

# Check if environment is set
def checkenvironment():
#    if executable == None:
#        return False
    return True

# Set environment
def setupenvironment(executable=None, scratch=None):
    if executable == None:
        print('Please provide an executable for NAMD')
        exit(1)

# Modify psf and prm files for QM/MM computations
def qmmmprepare(coords, structure, prm):
    srcdir = os.path.dirname(os.path.abspath(__file__))
    psfname = structure
    # set charges of QM atoms to 0.0
    structure = manipsf.setcharges(structure, coords.qmatoms, \
                        charge = 0.0, psfout = 'mod-' + psfname)
    # remove every bonded terms implying a QM atom
    structure = manipsf.removebonded(structure, coords.qmatoms, \
                        psfout = 'mod2-' + psfname)
    # remove vdw interactions in between QM atoms but not QM/MMvdw
    structure, newprm = manipsf.modatomtype(structure, prm, coords.qmatoms, \
                        psfout = 'qmmm-' + psfname, prmout = 'newtypes.prm')
    os.remove('mod-' + psfname)
    os.remove('mod2-' + psfname)
    prm.append(newprm)
    return structure, prm


# Print job file
def printjob(coords, jfile='namdjob', nproc=1, mem=1, meth='sp', pbc = None,     \
             structure = None, refcoords = None, prm = None, cutoff = None, \
             extra = None, nsteps = None, temperature = None):

    if meth == 'interactive':
        namdkey = defnkeyint
    else:
        namdkey = defnkey

    if structure == None:
        print('You have to provide a xplor type psf file to run NAMD')
        exit(1)

    if refcoords == None:
        print('You have to provide a coordinate file to run NAMD')   
        exit(1)

    namdkey['structure'] = structure
    namdkey['coordinates'] = refcoords
    namdkey['outputname'] = jfile + '-namd.out'

    if cutoff != None:
        namdkey['cutoff'] = str(cutoff)
        namdkey['switchdist'] = str(cutoff - 2.0)
        namdkey['pairlistdist'] = str(cutoff + 2.0)

    namdkey['bincoordinates'] = 'zonkey.bin.coor'

    if temperature != None:
        namdkey['temperature'] = str(temperature)

    # if given, add additional keyword provided
    if extra != None:
        for ex in extra:
            namdkey[ex[0]] = ex[1]

    fp = open(jfile + '.conf', 'w')
    
    for ip in prm:
        fp.write('parameters ' + ip + '\n')

    # print all options in the dictionary to the NAMD input (.conf) file
    namdoptions = ''
    for key in namdkey:
        namdoptions += key + ' ' + namdkey[key] + '\n'
    fp.write(namdoptions)

    if meth == 'sp':
        fp.write('run 0\n')
    # if gradient print TCL force script which makes namd print forces in a file
    elif meth == 'grad':
        fp.write('tclforces on\n' + \
                 'tclforcesscript {\n' \
                 '  set natoms ' + str(coords.natoms) + '\n' + \
                 '  set ee 0\n' + \
                 '  for {set i 1} {$i <= $natoms} {incr i 1} {\n' + \
                 '    addatom $i\n' + \
                 '  }\n' + \
                 '  enabletotalforces\n' + \
                 '  proc calcforces {} {\n' + \
                 '    set forcesfile [open "' + jfile + '-namdforces.txt" "w"]\n'+ \
                 '    global natoms\n' + \
                 '    loadtotalforces f\n' + \
                 '    global ee\n' + \
                 '    for {set i 1} {$i <= $natoms} {incr i 1} {\n' + \
                 '      if {[array exists f]} {\n' + \
                 '        puts $forcesfile $f($i)\n' + \
                 '      }\n' + \
                 '    }\n' + \
                 '    close $forcesfile' + \
                 '  }\n'  + \
                 '}\n' )
              # Make NAMD crash to avoid computing second step 
#             'if {$ee == 1} {\nexec MAKENAMDCRASH\n} else {\nset ee 1\n}' + \
        # run one step (otherwise gradients not computed)
        fp.write('\nrun 1\n')    
    elif meth == 'interactive':
        fp.write('tclforces on\n' + \
                 'tclforcesscript {\n' \
                 '  set natoms ' + str(coords.natoms) + '\n' + \
                 '  set ee 0\n' + \
                 '  set atoms {}\n' + \
                 '  for {set i 1} {$i <= $natoms} {incr i 1} {\n' + \
                 '    addatom $i\n' + \
                 '    lappend atoms $i\n' + \
                 '  }\n' + \
                 '  enabletotalforces\n' + \
                 '  proc calcforces {} {\n' + \
                 '    global atoms' + \
                 '    global natoms\n' + \
		 '    set interfile [open "INTERACTIVE" "w"]\n' + \
                 '    puts $interfile "QM"\n' + \
                 '    close $interfile\n' + \
                 '    set coorfile [open "intcoords.txt" "w"]\n' + \
                 '    loadcoords c\n' + \
                 '    for {set i 1} {$i <= $natoms} {incr i 1} {\n' + \
                 '      puts $coorfile $c($i)\n' + \
#                 '      puts stdout $c($i)   \n' + \
                 '    }\n' + \
                 '    close $coorfile\n' + \
                 '    while { 1 > 0} {\n' + \
                 '      after 200\n' + \
                 '      set interfile [open "INTERACTIVE" "r"]\n' + \
                 '      set check [ gets $interfile line ]\n' + \
                 '      if { $line == "MM"} {\n' + \
                 '        set gradfile [open "qmgradients.txt" "r"]\n' + \
                 '        foreach atom $atoms {\n' + \
                 '          set check [gets $gradfile line ]\n' + \
                 '          set force [ split $line ]\n' + \
#                 '          puts stdout \"$atom $force\" \n' + \
                 '          addforce $atom $force\n' + \
                 '        }\n' + \
                 '        close $gradfile\n' + \
                 '        close $interfile\n' + \
                 '        break\n' + \
                 '      }\n' + \
                 '      close $interfile\n' + \
                 ' ' + \
                 '    }\n' + \
                 '    set tstep [ getstep ] \n' + \
                 '    if { $tstep == ' + str(nsteps) + ' } {\n' + \
                 '      set interfile [open "INTERACTIVE" "w"]\n' + \
                 '      puts $interfile "STOP"\n' + \
                 '      close $interfile\n' + \
                 '    }\n' + \
                 '  }\n' + \
                 '}\n' )
        fp.write('\nrun ' + str(nsteps) + '\n')
    fp.close()

def runjob(coords, executable, refcoords=None, jfile='jf', nproc=1, \
           updatecoordinates=True):
    # use bin file because pdb is only with a few digits
    if updatecoordinates:
        binfile.printnamdbin('zonkey.bin.coor', coords.coords) 
    # run the actual job
    os.system(executable + ' +p ' + str(nproc) + ' ' + jfile + '.conf > ' + \
              jfile + '.log')

def runinteractivejob(coords, executable, refcoords=None, jfile='jf', nproc=1, \
           updatecoordinates=True):

    # use bin file because pdb coordinates only use a few digits
    if updatecoordinates:
        binfile.printnamdbin('zonkey.bin.coor', coords.coords)

    # run the actual job
    namdlog = open( jfile + '.log', 'w')
    if nproc > 1:
        # TODO modify this line
        namdijob = subprocess.Popen([executable, '+p', str(nproc), \
                         jfile + '.conf', '>', jfile + '.log'])
    else:
#        namdijob = subprocess.Popen([executable, jfile + '.conf', \
#                         '>', jfile + '.log'])
#        namdijob = subprocess.Popen([executable, jfile + '.conf > ' + jfile + '.log'])
        namdijob = subprocess.Popen(executable + ' ' + jfile + '.conf', stdout=namdlog, stderr=namdlog, shell=True)
        # if the job is no longuer runing close the log file
        if namdijob.poll() != None:
            namdlog.close()
    return namdijob

#TODO def checkjob()

def extractdata(coords, jfile='jf', val='energy'):

    if val in ['energy', 'gradient']:

        # look for the energy
        for line in open(jfile + '.log').readlines():
            if line[0:4] == 'ENER':
                break
        e = float(line.split()[13]) * KCALMOL2HARTREE
        # if single point return as no need for gradients
        if val == 'energy':
            return e
        
        # go for gradients 
        g = np.zeros((coords.natoms, 3))
        i = 0
        # -1.0 becauses its force and we use gradients
        conv = -1.0 * KCALMOL2HARTREE * ANG2BOHR
        # open the text file printed by the tclforce script 
        with open(jfile + '-namdforces.txt') as fp:
            for line in fp.readlines():
                l = line.split()
                g[i] = [float(l[0])*conv, float(l[1])*conv, float(l[2])*conv]
                i += 1
        return e, g

def clean(nametodelete, structurefile):
    extensions = ['.log', '-namdforces.txt', '-namd.out.xsc', '-namd.out.vel', '-namd.out.coor', \
                  '-namd.out.xsc.BAK', '-namd.out.vel.BAK', '-namd.out.coor.BAK', '.conf']
    filestodelete = ['zonkey.bin.coor', 'newtypes.prm']

    if structurefile[0:5] == 'qmmm-':
        filestodelete.append(structurefile)
    for f in nametodelete:
        for e in extensions:
            filestodelete.append(f + e)
    for f in filestodelete:
        if os.path.isfile(f): 
            os.remove(f)









