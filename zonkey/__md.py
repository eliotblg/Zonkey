
class Moleculardynamics(object):

    def __init__(self, options = None):
        self.opts = options

        self.execpath = epath
        if epath == None:
            self.execpath = None #TODO Here load from files and check
        self.scratchpath = spath
        self.name = name
        self.memory=memory
        self.nproc=nproc

    
