# NaCl
A Nacl salt solution provides a quick case study of what it takes to include [electrostatics](electrostatics.md) in a coarse grained model.

This case study assumes one already has a coarse-grained [water](../water/index.md) model constructed. It then represents each ion as individual beads. For simplicity, the ion smearing radii are taken to be the same as that of water.

The key is to use a smeared Coulombic interaction, with the Bjerrum length $l_b$ chosen to be the solvent's dielectric constant (at the desired temperature). For example, at room temperature the measured dielectric constant for water using the OPC water model is $\varepsilon_r = 78.4$. This leads to:

$$l_b = \frac{e^2}{4\pi\varepsilon_0\varepsilon_r k_BT} \approx 0.715\text{nm}$$

Details are discussed on the [electrostatics theory](../electrostatics.md) page, and an example [jupyter notebook](NaCl.ipynb) is also available.