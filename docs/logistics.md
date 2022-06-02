# Software installation, dependencies, and logistics
The tutorials and software we provide rely heavily on the following software packages. The functions of each piece of software are described, so that researchers may make substitutions where appropriate.


## Package Management
It is highly recommended to use a package manager to install the relevant code. Conda is recommended, and some packages are available only through pip (which comes with conda).

1. [conda](https://docs.conda.io/en/latest/)
2. [pip](https://pip.pypa.io/en/stable/)

## Molecular Dynamics
1. [OpenMM](https://openmm.org/)
    
    Particle simulations, used for both all-atom and coarse-grained simulations.

    ```
    conda install -c omnia openmm
    ```

2. [Packmol](http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml)
    
    For packing and initializing large molecular systems.

3. [MDTraj](https://mdtraj.org/1.9.4/api/generated/mdtraj.Trajectory.html)

    For manipulating and analyzing molecular dynamics trajectories

    `conda install -c conda-forge mdtraj`

## Coarse Graining
1. [Sim-Utils](https://github.com/EqualAPriori/sim-utils)

    Interface and utilities for running and driving the coarse-graining package used at BioPACIFIC. While the interface for manipulating and working with system definitions and force fields is public, for access to the underlying coarse-graining engine, contact BioPACIFIC.

## Field Theory
1. [PolyFTS](https://www.mrl.ucsb.edu/~fredrickson/)

    In-house field theory package maintained by Professor Glenn Fredrickson's group. For access, contact BioPACIFIC. 

## Other related tools
1. [LAMMPS](https://www.lammps.org/) - particle based simulation
2. [cgOpenMM](https://github.com/shirtsgroup/cg_openmm) - coarse grained simulations built on top of the OpenMM engine
3. [HOOMD-blue](http://glotzerlab.engin.umich.edu/hoomd-blue/) - GPU accelerated coarse grained particle simulations
4. [Freud](https://freud.readthedocs.io/en/latest/index.html) - analysis tools of self assembling systems
5. [VOTCA](https://www.votca.org/) - another open source coarse graining package
6. [MoSDeF](https://mosdef.org/) - a package of tools for setting up simulations.