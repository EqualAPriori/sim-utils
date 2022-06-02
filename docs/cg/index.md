# Coarse Graining a System

## I. Phenomenological models
### Ia) Model 1: Helfand compressibility with $\chi$ parameters
This is the typical model used in polymer field theory models.

### Ib) Model 2: Asymmetric bead volumes, the general excluded volume model
This model is a generalization of the Helfand compressibility model.

## II. Systematic coarse graining with the relative entropy
We provide standard tools for rapidly setting up coarse-grained simulations using the `mdfts` software package.

Core steps of coarse graining involve:

1. mapping
2. running relative entropy minimization
3. converting to field theory inputs

### IIa) Using the `sim` coarse graining package
Below is a short example showing all of the essential steps of setting up and running the `sim` package. The `sim-utils` package has utilities for automating a lot of the system definition steps via input files.

There are 7 main steps:

1. creating atom types
2. defining molecule types from the atom types
3. create a `world` container for the molecule types, and a `system` object containing the world and the potentials.
4. define potentials and add them to the system
5. set up molecular dynamics integrators, conditions, and initialize the system
6. compile the system
7. run a relative entropy optimization

```python
import sim

### ===== 1. Create atom aypes =====
atom_type_A = sim.chem.AtomType("A", Mass = 1.0, Color = (1,0,0))

### ===== 2. Define molecule types and bonding =====
mol_type_A = sim.chem.MolType("MA", [atom_type_a,atom_type_a,atom_type_a])
mol_type_A.Bond(0,1) #harmonic bond

### ===== 3. create world, system, and add molecules =====
world = sim.chem.World([mol_type_A], Dim = 3, Units = sim.units.DimensionlessUnits) #stores molecules
sys_name = "example_system"
sys = sim.system.System(world, Name = sys_name)
n_mols_A = 10
for i in range(n_mols_A):
    sys += mol_type_A.New()
sys.BoxL = 20

### ===== 4. make filters and create potentials =====
filter_aa_nonbond = sim.atomselect.PolyFilter([atom_type_a,atom_type_a], Bonded=False)
P = sim.potential.LJGaussian(sys, Cut = 3.0,
                     Filter = filter_aa_nonbond,
                     Epsilon = 0.0, Sigma = 1.0, B = 1.0, Kappa=1.0, Dist0=0.0,
                     Shift = True, Label = "LJG")

bondl,bondf = 1.0,10.0
filter_aa_bond = sim.atomselect.PolyFilter([atom_type_a,atom_type_a], Bonded=True)
p_bond = sim.potential.Bond(sys, Filter = filter_aa_bond, Dist0 = bondl, FConst = bondf, Label="Bonded")

sys.ForceField.extend([p_ljg,p_bond])

### ===== 5. Set up integrators and histogramming =====
integrator = sys.Int
integrator.Method = sys.Int.Methods.VVIntegrate
integrator.Method.Thermostat = Int.Method.ThermostatLangevin
integrator.Method.LangevinGamma = 1.0
integrator.Method.TimeStep = 0.01

sys.Measures.PEnergy.SetupHist(-1200, -500, 300)

### ===== 6. Compile system and give initial settings =====
sys.Load()
temperature = 1.0
sys.TempSet = temperature
sim.system.positions.CubicLattice(sys, Random = 0.1)
sim.system.velocities.Canonical(sys, Temp = TempSet)

### ===== 7. Set up optimizer and run Srel =====
trj = sim.traj.Simple("../sampletraj/ljtraj.trj.bz2") #load a zipped trajectory
trj.BoxL = sys.BoxL.copy()
opt = sim.srel.OptimizeTrajClass(sys, Beta = 1., Traj = trj, FilePrefix = "test")

opt.StepsEquil = 10000
opt.StepsProd = 500000
opt.StepsStride = 100
opt.Run()

### ===== 8. Optional: run the system =====
prefix = "test_md"
n_steps_minimization = 100
ret = sim.export.omm.MakeOpenMMTraj(sys, Verbose = True,
                                       NStepsMin = n_steps_minimization,
                                       NStepsEquil = 1000,
                                       NStepsProd = 10000,
                                       WriteFreq = 1000,
                                       Prefix = prefix,
                                       DelTempFiles = False)
traj, traj_file, dcd_file = ret
```
