# Core Simulations for Tiling Implementation on the Extended Hubbard Model

This is a set of simplified self-contained computer programs for simulations and computations on the extended Hubbard model by implementing the [newly developed](https://iopscience.iop.org/article/10.1088/1367-2630/ae36c8) **tiling auxiliary model approach**. The programs present in the [project repository](https://github.com/abhirup-m/Tiling-Embedded-Esiam) are more complicated because they involve additional code for plotting, parallelisation and other housekeeping tasks. This repository strips away all additional code and presents the minimal version necessary to implement the physics, making them pedagogically more suitable.

## How to read this
1. *Unitary renormalisation group (RG) analysis of the lattice-embedded impurity model*
Start from here, by looking through the `RGFlow.jl` file. This file sources the `Helpers.jl` file and calls its functions whenever necessary, so look through it as and when necessary.

2. *RG Phase diagram calculation*
After the RG flow, the next step is generating a phase diagram from the fixed point couplings. This is carried out in `PhaseDiagram.jl`, and invokes the `RgFlow.jl` script.

- Impurity model $k-$space static correlation calculation (spin, charge)
- Impurity model local spectral function and quasiparticle residue
- Impurity model local self-energy and optical conductivity
- Tiled real-space static correlation calculation (spin, entanglement measures, QFI)
- Tiled $k-$space spectral function and self-energy
