# Bonded interactions
## Harmonic Bonds
A typical bonded interaction is of the form:
$$
u(r) = \frac{k}{2}(r-r_0)^2
$$
where $k$ is the bond strength and $r_0$ is a reference bond length. There are two helpful limits: 

1.) $k$ large, which approximates a **freely jointed chain** model with bond length $r_0$.

2.) a zero-centered $r_0=0$ model, which is the **discrete Gaussian chain** model, with average bond length $b=(3/k)^{1/2}$

Zero-centered bonds are especially recommended because 1) they are softer and allow for larger coarse grained MD time steps, and 2) the optimization problem for zero-centered bonds is faster, less stiff, and easily allows for co-optimization with other target properties like chain end-to-end distance.

## Angular Bonds
While angular bonded potentials can be obtained easily with the coarse-graining procedure, they are expensive to compute in the field theory and are thus usually ignored.

