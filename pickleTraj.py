def pickleTraj(TrajName, LogFile = None, LogFileToken = None, Verbose = False):
    import os, cPickle as pickle
    import sim
    
    pickleName = TrajName + '.pickle'
    if os.path.isfile(pickleName):
        of = open(pickleName, 'r')
        if Verbose: print 'Loading from pickle...'
        Trj = pickle.load(of)
        of.close()
    else:
        if Verbose: print 'Pickling trajectory...'
        Trj = sim.traj.lammps.Lammps(TrajName, LogFile = LogFile, LogFileToken = LogFileToken)
        of = open(pickleName, 'w')
        pickle.dump(Trj, of)
        of.close()
    
    init = Trj[0] # needed to enable framedata parsing
    return Trj

import sys
sys.modules[__name__] = pickleTraj
