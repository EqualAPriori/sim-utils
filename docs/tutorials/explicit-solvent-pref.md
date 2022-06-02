# Explicit Solvents and (Coarse-Grained) Reference Pressures
In the spirit of coarse graining, it is recommended to coarse grain to soft interactions that smear out short length scale behaviors. The motivation is that, if possible, one does not want to spend computational resources resolving short length scale molecular structuring. Note that in principle this does not preclude one from retaining thermodynamically appropriate behavior on longer length scales. For example, the virial pressure is calculated as a weighted integral over a pair distribution function, and the same pressure can be reproduced by many different pair correlation functions. In this spirit, coarse grained models with soft short-range interactions only seek to preserve longer-length scale structural correlations and thermodynamic behaviors.

## The Gaussian Interaction Model
The interaction potential of choice is a Gaussian repulsive potential:
$$
u_{ij}(r) = \frac{u_{ij}}{(2\pi(a_i^2+a_j^2))^{3/2}}e^{-r^2/(a_i^2+a_j^2)}
$$
where $u_{ij}$ is the interaction repulsion strength and $a_i,a_j$ are the smearing interaction widths of beads $i$ and $j$, respectively.

In the limit where all smearing lengths are the same and all $u_ii$ are the same, this model can be exactly mapped to smeared molecular models that use the Helfand compressibility and a $\chi$ parameter.

As illustrated in [Shen et al. 2020](https://aip.scitation.org/doi/full/10.1063/5.0022808) and [Sherck et al 2021](https://pubs.acs.org/doi/10.1021/acsmacrolett.1c00013), using such potentials greatly attenuates the repulsive core in pair distributions. This allows for particle overlap and larger time steps in particle simulations. Further, with appropriately chosen repulsive potential strength $u_0$, miscibilities and macroscopic phase behavior can be reproduced, sometimes quantitatively [(Shen et al. 2020)](https://aip.scitation.org/doi/full/10.1063/5.0022808).

In practice, families of coarse grained models can be calibrated to what is termed in the DPD literature as a "reference coarse-grained pressure" (Fraije) that is not necessarily equal to that of the all-atom system. The intuition is that in nearly incompressible systems (i.e. most solution formulations and solid materials), the actual pressure of the system is usually not of great importance. Rather, it is more important to simply reproduce this nearly incompressible behavior. For example, if modeling water, choosing $a_w=0.3nm$ and $u_{ww}=20kT*a_w^3$ leads to a reference coarse-grained pressure of $P_{cg,ref}\approx10.31kT/a_w^3$ and compressibility $0.062_w/kT\approx4.45\cdot 10^{-10}Pa^{-1}$ that is near that of water at room temperature. Again, for these repulsive Gaussian models, most likely $P_{cg,ref}$ does *not* equal that of the all-atom system. As long as the pressure-dependences are not of interest, coarse grained simulations run at $P_{cg}$ model simulations of the all-atom system at constant pressure (e.g. $1$atm for many aqueous solution studies). Further, if pressure changes *are* of interest, as long as the compressibility is reproduced, it still makes sense to compare AA simulations studied at pressure $P=P_{AA,ref}+\Delta P$ to CG simulations studied at $P=P_{cg,ref}+\Delta P$.

Models coarse grained to this $P_{ref}$ are said to be consistent and calibrated. This follows the spirit of [dissipative particle dynamics (DPD)](https://aip.scitation.org/doi/10.1063/1.474784), which use soft, truncated quadratic potentials instead of the soft Gaussian potentials employed in field theory. DPD also requires choosing a characteristic $P_ref$ to operate at.


## Estimating Initial and Choosing Parameters
With sufficiently large smearing lengths, mean field analytical approximations give really good estimates of the pressure of dense, nearly incompressible systems like liquids.

For example, for a homopolymer with length $N$, the mean-field pressure is:

$$
\beta P = \frac{\rho}{N} + \frac{1}{2}u_{ii}\rho^2
$$

where $\rho$ is the bead density. The compressibility in this model is:
$$
\kappa = \left( \frac{\rho}{N} + u_{ii}\rho^2 \right)^{-1}.
$$

For water, choosing $u_{ii}=15v_w$ and $\rho= 1/v_w$, where $v_w$ is the water volume, gives $P_{ref}=8.5kT/v_w$ and a compressibility $k_BT/16v_w\approx4.5\times10^{-10}$Pa, which roughly reproduces the compressibility water at room temperature. Note that since water is quite incompressible at room temeprature, it can actually be safe to [reduce the compressibility of the coarse grained water by even a factor of 3](https://pubs.acs.org/doi/pdf/10.1021/acsmacrolett.1c00013).