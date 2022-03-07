# NaCl
A Nacl salt solution provides a quick case study of what it takes to include [electrostatics](electrostatics.md) in a coarse grained model.

This case study assumes one already has a coarse-grained [water](../water/index.md) model constructed. It then represents each ion as individual beads. For simplicity, the ion smearing radii are taken to be the same as that of water.

The key is to use a smeared Coulombic interaction, with the Bjerrum length $l_b$ chosen to be the solvent's dielectric constant (at the desired temperature). For example, at room temperature the measured dielectric constant for water using the OPC water model is $\varepsilon_r = 78.4$. This leads to:

$$l_b = \frac{e^2}{4\pi\varepsilon_0\varepsilon_r k_BT} \approx 0.715\text{nm}$$

Theoretical details are discussed on the [electrostatics theory](../electrostatics.md) page, code and files can be found in the [repository](https://github.com/EqualAPriori/sim-utils/tree/scout/docs/tutorials/salt) for all the files, and an example [jupyter notebook](NaCl.ipynb) is also available. The novel feature is the addition of the smeared electrostatic potential. In the jupyter notebook this is set via the configuration files. To do it manually, one needs to first make sure that ion atom type definitions are defined with charges:

```python
Na_type = sim.chem.AtomType("Na", Mass=1.0, Charge = 1.0)
Cl_type = sim.chem.AtomType("Cl", Mass=1.0, Charge =-1.0)
```

And subsequently, to add both the bare Coulomb correction and the smeared Coulomb correction *to* the Coulomb interaction.
```python
p = sim.potential.Ewald( Sys, Cut=2.5, Shift=True, Coef=0.715, FixedCoef=True, Label='ewald' )
pcorr = sim.potential.smearedcoulombEwCorr( Sys, Coef = 0.715, Cut = 2.5, BornA = 0.55, FixedCoef = True, FixedBornA = True, Label="smearedCoulCorr" )
Sys.ForceField.extend([p,pcorr])
```

The jupyter notebook further shows how if the charge parameters are held fixed, the simulation can be further accelerated by using a proxy system without electrostatic interactions to calculate energy differences with variations in non-bonded interactions. (since the electrostatic energy parameters are held constant).
