### Helper classes and functions to create an optimizer
import os
import sim
from collections import namedtuple

from dotdict import dotdict
import forcefield as ff
import numpy as np

### wishlist
# 1) extended ensemble
# 2) Hessian recalculation

### Create optimizer
def create_optimizer(Sys, traj, md_engine, steps_equil, steps_prod, steps_stride, ElecSys=None):
    # === create a mapping object (1:1, since actual mapped traj created earlier) ===
    Map = sim.atommap.PosMap()
    print(Sys.Name)
    print('NMol: {}'.format(Sys.NMol))
    print('NAtom: {}'.format(Sys.NAtom))
    print('NDOF: {}'.format(Sys.NDOF))
 
    for (i, a) in enumerate(Sys.Atom):
        Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]

    # === Read in the Trajectory ===
    # set box-len same as AA traj
    Sys.BoxL = traj.FrameData['BoxL']
    print("Box: {}".format(Sys.BoxL))

    # === Set up Optimizers ===
    if md_engine.lower() in ['omm','openmm']:
        print('... using OpenMM engine')
        OptClass = sim.srel.optimizetrajomm.OptimizeTrajOpenMMClass
    elif md_engine.lower() in ['lammps']:
        print('... using LAMMPS engine')
        OptClass = sim.srel.optimizetrajlammps.OptimizeTrajLammpsClass
    else:
        print('... default using sim engine')
        OptClass = sim.srel.optimizetraj.OptimizeTrajClass
        fobj = open('tmp_measures{}.dat'.format(Sys.Name), 'w')
        Sys.Measures.VerboseOutput(fobj = fobj, StepFreq = StepsStride)
        Int = Sys.Int
        Trj = sim.traj.lammps.LammpsWrite("tmp_traj_{}.lammpstrj".format(Sys.Name))
        Trj.AddAction(Int, StepFreq = steps_stride)
 
    # === Create Optimizer ===
    ewaldPotentials = [ P for P in Sys.ForceField if P.Names[0] == 'ewald' ]
    UseTarHists = False
    if not ewaldPotentials: #if no ewald
        UseTarHists = True
    print("... use target histograms: {}".format(UseTarHists))

    Opt = OptClass(Sys, Map, Beta=1./Sys.TempSet, Traj=traj, FilePrefix=Sys.Name, ElecSys=ElecSys,
            SaveLoadArgData=True, TempFileDir=os.getcwd(), UseTarHists=UseTarHists)

    Opt.ConstrainNeutralCharge()

    # md iterations
    Opt.StepsEquil = int(steps_equil) #int(steps_equil/Sys.Int.Method.TimeStep)
    Opt.StepsProd  = int(steps_prod) #int(steps_prod/Sys.Int.Method.TimeStep)
    Opt.StepsStride = int(steps_stride) #int(steps_stride/Sys.Int.Method.TimeStep)

    # extra prep for the optimizer object
    sim.srel.optimizetraj.PlotFmt = 'png'
    Opt.MinReweightFrames = 10
    #Opt.MinReweightFrac = 0.15

    # === Diagonstic, check whether parameters are fixed ===
    print('...check fixed parameters status:')
    print "i,Fixed"
    for potential in Sys.ForceField:
        print("--- {} ---".format(potential.Label))
        for (i,Fixed) in enumerate(potential.Param.Fixed):
            print i,Fixed

    return Opt
"""
    if recalc:
        #only manually unfix external potential. fixing/unfixing other variables should be done in main or recalc.py
        PExtSin = [ force for force in Sys.ForceField if force.Names[0]=='external_sinusoid' ]
        PExtSp = [ force for force in Sys.ForceField if force.Names[0]=='external_spline' ]
        for pext in PExtSin:
            pext.UConst.Fixed = False #pext.Fixed = False
        for pext in PExtSp:
            pext.Fixed = False #pext.Fixed = False

        PExtSin = [ force for force in ElecSys.ForceField if force.Names[0]=='external_sinusoid' ]
        PExtSp = [ force for force in ElecSys.ForceField if force.Names[0]=='external_spline' ]
        for pext in PExtSin:
            pext.UConst.Fixed = False #pext.Fixed = False
        for pext in PExtSp:
            pext.Fixed = False #pext.Fixed = False

        #other optimizer set up
        Opt.CheckReady()
        Opt.StartIter = Opt.Iter
        SearchDirIter = 0
        Opt.Mode = "INIT"
        Opt.Backtracking = False
        Opt.UseMaxChange = False
        Opt.PrevLineSearch = False
        Opt.CGRestartCount = 0
        Opt.Initdx()

        Opt.InitConstraints(Verbose = Opt.Verbose)
        Opt.SetParam(Opt.Param)
        if Opt.Verbose:
            print Opt.Output0()
        Opt.OutputTarHistFile()
        Opt.ReweightTar = False #i.e. make a new model trajectory
"""




def run_opt(Opts, Weights, Name, UseWPenalty, MaxIter, SteepestIter, StageCoefs):
    Optimizer = sim.srel.OptimizeMultiTrajClass(Opts, Weights=Weights)
    Optimizer.FilePrefix = ("{}".format(Prefix))
    if not UseWPenalty:
        Optimizer.RunConjugateGradient(MaxIter=MaxIter, SteepestIter=SteepestIter)
    else:
        Optimizer.RunStages(StageCoefs = StageCoefs)

def recalc(Opt,prefix_append='recalc'):
    import time
    StartTime = time.time()
    #Opt = Opts[0] 
    #Opt = sim.srel.OptimizeMultiTrajClass(Opts) #creating a new object would require not just appending, but giving new file prefix name
    Opt.FilePrefix += '_' + prefix_append
    print('saving recalculation with prefix {}'.format(Opt.FilePrefix))

    # set recalculation options 
    Opt.CheckReady()
    Opt.StartIter = Opt.Iter
    SearchDirIter = 0
    Opt.Mode = "INIT"
    Opt.Backtracking = False
    Opt.UseMaxChange = False
    Opt.PrevLineSearch = False
    Opt.CGRestartCount = 0
    Opt.Initdx()

    Opt.InitConstraints(Verbose = Opt.Verbose)
    Opt.SetParam(Opt.Param)
    if Opt.Verbose:
        print Opt.Output0()
    Opt.OutputTarHistFile()
    Opt.ReweightTar = False #i.e. make a new model trajectory

    # === begin calculation ===
    print('===== Begin recalculating trajectory and parametric derivatives =====')   
    if not Opt.ReweightTar:
        print("OldModTraj is: {}".format(Opt.OldModTraj))
        print("ReweightOldModTraj was: {}".format(Opt.ReweightOldModTraj))
        print("...setting to False")
        Opt.ReweightOldModTraj =  False
        print("...ReweightOldModTraj is now: {}".format(Opt.ReweightOldModTraj))
        Opt.UpdateModTraj()
        Opt.OutputModHistFile()
        Opt.OutputPlot()

    NotFixed = np.logical_not(Opt.Fixed)
    Opt.CalcObj(CalcDeriv = True)

    print("=== Hessian/DDSrel ===")
    np.savetxt("{}_Hessian.dat".format(Opt.FilePrefix),Opt.DDObj)

    Ignore = Opt.Fixed.copy()
    for i in Opt.Constrained:
        Ignore[i] = True

    for (i, ThisDObj) in enumerate(Opt.DDObj):
        if all(Opt.DDObj[i,:] == 0):
            Ignore[i] = True

    Use = np.logical_not(Ignore)

    DObj = Opt.DObj[Use]
    DDObj = Opt.DDObj[np.ix_(Use, Use)] + Opt.HessianPad * np.identity(len(DObj))

    print('Used variables: {}'.format(Use))
    print('Non-zero Hessian:\n{}'.format(DDObj))
    np.savetxt('{}_Hessian_masked.dat'.format(Opt.FilePrefix),DDObj)

    print("=== Gradient/DSrel ===")
    np.savetxt('{}_grad.dat'.format(Opt.FilePrefix),Opt.DObj)
    print('Non-zero Grad:\n{}'.format(DObj))
    np.savetxt('{}_grad_masked.dat'.format(Opt.FilePrefix),DObj)

    print("=== Misc. ===")
    print("Total recalculation runtime: {}min".format( (time.time()-StartTime)/60. ))


