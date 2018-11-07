import os
import nlopt
import numpy as np
import imp
from . utils import manipulatepdb
from . utils import printfile

class Geopt(object):
    def __init__(self, ofunction, algorithm='LD_LBFGS', historyf='geopt-history.xyz'):
        self.ofunction = ofunction
        self.algorithm = 'nlopt.' + algorithm.upper()
        # if there is already a history file from previous run move it to BCKP-
        if os.path.isfile(historyf):
            os.rename(historyf, 'BCKP-' + historyf)
        self.histfile = historyf 

    def objectivefunction(self, x, grad):

        print('*** Geometry optimization with ' + self.algorithm + \
              ' | Step: ' + str(self.nstep) + '/' + str(self.opt.get_maxeval()))  

        self.setcoordinates(self.c, x)
        # if gradient optimization (normaly the case)
        if grad.size > 0:
            e, fullgrad = self.ofunction.gradients(self.c)
            grad[:] = self.getgradientsarray(self.c.activelist, fullgrad)
        # if optimization algorithm just needs energy
        else:
            e = self.ofunction.energy(self.c)
              
        de  = e - self.previouse
        dde = de - self.previousde
        print('Delta  E(Hartree) : ' + "% .2E"%(de) + \
              ' | ' + "% .2E"%(self.opt.get_ftol_abs()))
        print('Delta DE(Hartree) : ' + "% .2E"%(dde) + \
              ' | ' + "% .2E"%(self.opt.get_ftol_rel()))

        if grad.size > 0:
            madg, maxg = self.getmadmaxgrad(grad)
            print('Max grad compo(au): ' + "% .2E"%(maxg) + \
                  ' | ' + "% .2E"%(self.gmaxtol))
            print('Mad grad compo(au): ' + "% .2E"%(madg) + \
                  ' | ' + "% .2E"%(self.gmadtol))
            if maxg < self.gmaxtol and madg < self.gmadtol:
                print('Max and Mad Gradients converged at step ' + str(self.nstep))
                self.opt.force_stop() 

        printfile.printxyz(self.c.coords, self.c.atypes, \
                           self.histfile, append = True)

        self.previouse  = e
        self.previousde = de
        self.nstep += 1
        return e

    def setoptimizer(self, n, maxeval = 200, step = 0.1, rtol = 1.0e-10, \
                     atol = 1.0e-8, gmaxtol = 4.0e-4, gmadtol = 5.0e-4,  \
                     algorithm = None):
        self.dim = n * 3
        if algorithm == None:
            self.opt = nlopt.opt(eval(self.algorithm), self.dim)
        else:
            self.opt = nlopt.opt(eval('nlopt.' + algorithm.upper()), self.dim)
        self.opt.set_initial_step(step)
        self.opt.set_min_objective(self.objectivefunction)
        self.opt.set_maxeval(maxeval)
        self.opt.set_ftol_rel(rtol) 
        self.opt.set_ftol_abs(atol)
        self.gmaxtol = gmaxtol
        self.gmadtol = gmadtol
        self.nstep = 1

    def optimize(self, coords):
        # load data
        x = self.getactivearray(coords)
        self.c = coords
        self.previouse  = 0.0
        self.previousde = 0.0

        # run optimization; ForcedStop signal is sent when converged
        try:
            newx = self.opt.optimize(x)
        except nlopt.ForcedStop:
            newx = self.opt.last_optimize_result()
        # extract data
        e = self.opt.last_optimum_value()
        print('Energy after optimization: ' + "%20.10f"%(e) + ' Hartrees')
        coords = self.c
        manipulatepdb.changecoordinates(coords.infile, list(range(coords.natoms)),\
                                coords.coords, newfile=coords.infile[0:-4] + '-opt.pdb')
        printfile.printxyz(coords.coords, coords.atypes, \
                           coords.infile[0:-4] + '-opt.xyz')

    def getactivearray(self, coords):
        x = np.array([])
        for i in coords.activelist:
            x = np.append(x, coords.coords[i])
        return x

    def setcoordinates(self, coords, x):
        j = 0
        for i in coords.activelist:
            for ii in range(3):
                coords.coords[i][ii] = x[j]
                j += 1

    def getgradientsarray(self, indexlist, fullgrad):
        g = np.array([])
        for i in indexlist:
            g = np.append(g, fullgrad[i])
        return g

    def getmadmaxgrad(self, grad):
        t = np.absolute(grad)
        return np.mean(t), np.amax(t)  

    # simple steepest (gradient) descend just for testing
    def steepestd(self, coords, step = 0.1, maxeval = 200):
        self.c = coords
        pc = coords.coords
        de = 1E9
        pe = 1E9
        for i in range(maxeval):
            self.c.coords = pc
            e, g = self.ofunction.gradients(self.c)
            de = e - pe
            pe = e
            print('Energy and DE at step ' + str(i+1) + ': ' + str(e) +' \ ' + str(de) )
            if abs(de) < 1.0E-8:
                print('Converged after ' + str(i+1) + ' steepest descent steps\n')
                break
            c = pc - step * g           
#            if i > 0: 
#                dg = np.ravel(g-pg)
#                #ndg = np.sum(np.absolute(dg))
#                ndg = np.linalg.norm(dg)
#                #print(dg)
#                #print(np.transpose(np.array(np.ravel(c-pc))))
#                step = np.matmul(np.transpose(np.ravel(c - pc)), dg) / (ndg*ndg)
#                print('step size corrected to: ' + str(step))
            printfile.printxyz(c, coords.atypes, self.histfile, append = True)
            pc = c
            pe = e
            pg = g 
            

















