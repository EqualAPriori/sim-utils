# Electrostatic interactions
Electrostatic interactions are important for maintaining the aqueous solubility of many polymeric and surfactant molecules, and are also intimately involved in many self assembly processes, driving electrostatic complexation.

## Smeared Electrostatics
The simplest electrostatic interaction compatible with both particle-based molecular simulation and field theory is a *smeared* electrostatic interaction:

$$
U_{elec}    = \frac{l_b}{2} \sum_{ij} \frac{1}{r_{ij}}erf\left( \frac{r_{ij}}{2\sqrt{a_i^2/2+a_j^2/2}}\right)
$$

Where $a_i$ and $a_j$ are the electrostatic smearing length. It is conventional to take this smearing length to be the same as the smearing length in the [nonbonded interactions](explicit-solvent-pref.md). This is assumed in the field theory models, although it is not strictly necessary. The main advantage of using a smeared electrostatics model are many: 1) it avoids singularities in the energy when two charged species overlap, and 2) it allows a field theory model to physically capture a combination of both [a Born solvation and correlation energies](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.81.021501). Roughly speaking, the solvation energy (will come out to $l_b/2\pi^{1/2}a$ ) represents the energy to solvate a single ion in solvent, while the correlation energy is a many-body effect describing the interaction of ions with other ions.

The $l_b$ factor is the Bjerrum length

$$
l_b = \frac{e^2}{4\pi\epsilon_0\epsilon_r k_BT}
$$

which characterizes the electrostatic interaction strength. The dielectric constant $\epsilon_r\approx80$ in water at room temperature, leading to $l_b\approx 0.7$nm. For implicit solvent water models it is natural to choose $\epsilon_r\approx80$ to model the electrostatic effects of water. For the explicit solvent models, if we do not describe the solvent molecules with dipoles (see following section), the dielectric screening effect must be taken into account by again choosing an appropriate dielectric constant $\epsilon_r$ or $l_b$.

The Bjerrum length should be adjusted for every new temperature that a model is parameterized at. It is customary to choose $l_b$ to be the solvent's dielectric constant, as this exactly reproduces the dilute-limit electrostatic interaction between two charges. This approximation is good in the limit where a solvent is dominant. There are no standard procedures for interpolating between two environments with very different dielectric constants, although approximations can be made by [adding additional $\chi$ interactions to model the solvation penalty](https://pubs.acs.org/doi/abs/10.1021/acs.macromol.1c00095) of moving from a high dielectric to low dielectric environment. In a similar vein, the relative entropy minimization will naturally try to account for these electrostatic effects in effective nonbonded interactions.

This smeared electrostatic interaction has been used in [field theory](https://www.pnas.org/doi/10.1073/pnas.1900435116) and also in [dissipative particle dynamics models](https://arxiv.org/abs/1608.00626).

## Dipolar interactions
The field theory also allows for [dipolar interactions](https://aip.scitation.org/doi/full/10.1063/1.4964680). Such field theoretic descriptions allow for explicit dipoles on solvent molecules to naturally generate solvation effects and dielectric screening effects. However, these can not be efficiently and exactly implemented in particle-based simulations, so we do not consider them in the context of the all atom to field theory work flow.

