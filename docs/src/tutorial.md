# Tutorial

The basic usage is probably best illustrated with a brief example. In the following, we solve a simple model for a photon pulse scattered by a quantum dot in a chiral waveguide. 

We start by loading the package, defining some symbolic parameters and the photonic annihilation operator `a` as well as the atomic transition operator `σ`, which denotes a transition from level `j` to level `i` as `σ(i,j)`. This allows us to quickly write down the Hamiltonian and the collapse operators of the system with their corresponding decay rates.

```@example tutorial
using Latexify # hide
set_default(double_linebreak=true) # hide
using QuantumInputOutput
```
