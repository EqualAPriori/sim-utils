# Water
Water provides a case study in understanding the [main interaction model](../explicit-solvent-pref.md) underlying simulations suitable for running field theory.

We coarse grain water using a 1 bead model:

![water](water.png){: style="height:111px;width:305px"}

See [this notebook](water.ipynb) for example code coarse graining water at $25^o$C. It fixes the water smearing radius to $a_w=0.3107$nm $=\rho_w^{-1/3}$, and determines a suitable repulsion strength $u_w$.