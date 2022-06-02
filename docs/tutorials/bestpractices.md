# Best practices, in a nutshell
This section collects best practices and tutorials, organized by topic. It is mean to be a quick guide for finding relevant example molecular systems and theoretical considerations

Most realistic systems are highly multicomponent, and involve several subtle coarse graining steps and decisions. 


## I. Choosing a coarse-grained interaction model
The following tutorials emphasize the theoretical considerations when constructing models amenable for running field theory. 

### 1. Use explicit solvent models with soft interactions

For a sample system, see this the [pure water](water/index.md) case study.

For theoretical considerations, see [here](explicit-solvent-pref.md).

In the spirit of coarse graining, it is recommended to use soft interactions that smear out short length scale structure and behaviors, and obtain coarse grained parameters that still reproduce larger scale structural and thermodynamic properties.

The interaction potential of choice is a Gaussian repulsive potential:
$$
u_{ij}(r) = \frac{u_0}{(2\pi(a_i^2+a_j^2))^{3/2}}e^{-r^2/(a_i^2+a_j^2)}
$$
where $u_0$ is the interaction repulsion strength and $a_i,a_j$ are the smearing interaction widths of beads $i$ and $j$, respectively.


### 2. Use soft bonded interactions
For a sample system, see the [pure dodecane](dodecane/index.md) case study.

For theoretical considerations, see [here](bonded.md)

A typical bonded interaction is of the harmonic form:
$$
u(r) = \frac{k}{2}(r-r_0)^2
$$
where $k$ is the bond strength and $r_0$ is a reference bond length. 

### 3. Charged interactions
For a sample system, see the [NaCl](salt/index.md) case study.

For theoretical considerations, see [electrostatics](electrostatics.md). It is recommended to use 

$$
U_{elec}    = \frac{l_b}{2} \sum_{ij} \frac{1}{r_{ij}}erf\left( \frac{r_{ij}}{2\sqrt{a_i^2/2+a_j^2/2}}\right)
$$

Where $a_i$ and $a_j$ are the electrostatic smearing length and the $l_b$ factor is the Bjerrum length

$$
l_b = \frac{e^2}{4\pi\epsilon_0\epsilon_r k_BT}
$$
with $\epsilon_r$ the dielectric constant of the solvent.

## II. Collecting appropriate all-atom simulation data
The collection of appropriate all-atom molecular trajectories is of paramount importance in obtaining reasonable coarse grained models. It is well-documented that systematic coarse graining can be [highly dependent on the all-atom simulations used to derive the models](https://aip.scitation.org/doi/full/10.1063/5.0022808). The tutorials below help guide one through that process.

### 1. Homogeneous, bulk systems
For a sample system, see the [pure water](water/index.md) case study.

As described in the [theoretical considerations for explicit solvent](explicit-solvent-pref.md), it is highly recommended that bulk systems be coarse grained in the isobaric (NPT) ensemble at an appropriately chosen coarse grained reference $P_{cg}$. Usually, $P_{cg}$ does not need to match the pressure of the all atom system, as long as the compressibility is approximately right (what matters are volume changes with changing pressure and composition, not the absolute pressure). For typical applications, $P_{AA}=1$bar and $P_{cg}=8.5 k_BT/v_w$ is a good choice if one wants to use fairly incompressible water as a coarse grained reference.

### 2. Mixtures and multiple bead types

Mixtures, which involve multiple bead types, needs some extra care. [Recent work](https://aip.scitation.org/doi/full/10.1063/5.0022808) has shown that naive coarse graining in the NVT ensemble can lead to very poor predictions of solution thermodynamics like activity coefficients. The key is to simulate from molecular configurations and ensembles where thermodynamically relevant collective variables are sensitive to the underlying molecular interactions.

### 3. Miscible systems
For a sample system, see the [methanol-water](methanol_water/index.md) case study.

For a miscible system, the density is a relevant collective variable conjugate to the activity coefficient. The density response can in turn be probed by coarse graining using an [external potential ensemble](https://aip.scitation.org/doi/full/10.1063/5.0022808).


### 4. Immiscible systems
For a sample system, see the [dodecane](dodecane/index.md) case study.

For an immiscible system like dodecane-water, the system already naturally forms two segregated phases. In practice these systems are coarse-grained from a bulk-phase-segregated state. The appropriate ensemble simply barostats the axis normal to the interface, while fixing the dimensions parallel to the interface. In OpenMM, the interface-normal should be the $z-axis$.

### 5. Structured systems
For a sample system, see the [SDS](SDS/index.md) case study.

For a structured system like one involving a surfactant, one can experiment with coarse graining from the various structures formed by the molecule. It is suggested to coarse grain from well-reproducible structures, for example micelles at the experimentally measured aggregation number, or surfactants deposited at a planar interface between two fluids.


## Running Field Theory
### 1. Phase separation estimation: [Mean Field Theory](theory/meanfield.md)
For documentation, see [here](https://pubs.acs.org/doi/pdf/10.1021/acsmacrolett.1c00013).

For a sample system, see the [PEO](PEO/index.md) case study.

Especially for homogeneous phases, the mean field theory can provide quick analytical estimates of molecular compatibility. 

### 2. Phase separation estimation of charged systems: coacervation and [RPA theory](theory/rpa.md)
Systems undergoing phase separation induced by charge-complexation (as opposed to other causes like hydrophobicity) can be described using RPA theory, which describes attractions due to field fluctuations, or charge charge correlations.

For theoretical considerations, see this notebook studying [PAA-PAH](PE/index.md).

### 3. Comparing structured phase stability with self consistent field theory
For a sample system, see this notebook studying [SDS](SDS/index.md).


