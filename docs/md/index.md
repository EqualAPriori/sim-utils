# Particle-based Molecular Dynamics

## All atom simulations
For now we defer to extensive guides provided [elsewhere](http://docs.openmm.org/latest/userguide/application/02_running_sims.html). However, the key features that one needs to be familiar with include:

1. **Creating force fields and partial charges**. This can either be done with a software package like [AmberTools](https://ambermd.org/AmberTools.php), downloaded from papers that provide force fields, or obtained from force field repositories like [OpenKim](https://openkim.org/model-developer-directory/) or the [Open Force Field Initiative](https://openforcefield.org/).
2. **Create Initial Configurations**. This can be done, e.g. with tools like [packmol](http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml). Standard simulation tools like [AmberTools](https://ambermd.org/AmberTools.php) often have utility functions that aid in this process as well.
3. **Create a system topology**. In essence, this is an enumeration of bead type names and bonding topologies for each molecule type of interest. This can be done through manual curation with tools like [Avogadro](https://avogadro.cc/) or through full stack tools like [AmberTools](https://ambermd.org/AmberTools.php).
4. **Run**. With the above data in a suitable format, one can use any all atom simulation package of preference, such as [OpenMM](https://openmm.org/) or [Gromacs](https://www.gromacs.org/). Tools like [ParmEd](https://github.com/ParmEd/ParmEd) and [MoSDeF](https://mosdef.org/) can help with conversion between file formats.

## Coarse-grained simulations
We provide standard tools for rapidly setting up coarse-grained simulations using the `mdfts` software package and format for defining force fields and topologies.

