# [Theoretical Background](@id theory)

This section summarizes the theoretical ideas behind **QuantumInputOutput.jl**, i.e. input-output theory for quantum pulses, its formulation via SLH networks, and the use of virtual cavities to handle traveling wave packets. We also explain how temporal output-mode bases are determined using two-time correlation functions.

## Input-Output Theory with Quantum Pulses

The central goal of the input-output theory with quantum pulses is to describe the dynamics of a photon pulse, which interacts with a system. For such a system one is usually interested in the pulse shape (temporal mode) and the quantum state of the pulse, as well as the quantum state of the system during and after the interaction. 
We consider a local quantum system with Hamiltonian $H_s$ coupled to a traveling bosonic field $\hat b_{in}(t)$
through a system operator $\hat c$. Under standard Markov and dispersionless-propagation assumptions, the
interaction can be written as

```math
\hat V_{SB}(t) = i \sqrt{\gamma}\,\big(\hat c\,\hat b_{in}^\dagger(t) - \hat c^\dagger \hat b_{in}(t)\big).
```

The outgoing field satisfies the input-output relation

```math
\hat b_{out}(t) = \hat b_{in}(t) + \sqrt{\gamma}\,\hat c(t).
```

Quantum pulses are not single-frequency modes but wave packets occupying a continuum of modes. A normalized
temporal mode $u(t)$ defines a creation operator

```math
\hat b_u^\dagger = \int dt\, u(t)\,\hat b^\dagger(t).
```

When the system interacts with such a pulse, the outgoing radiation can become multimode. The goal is to
describe the full system–field dynamics in a tractable way, without explicitly discretizing the entire continuum.

## Virtual Cavities and Pulse Modes

The key idea is to replace the traveling input and output pulses by **virtual cavities**, each with a time-dependent, complex coupling to the system determined by the corresponding temperal mode.  

For a normalized input mode $u(t)$, the coupling

```math
g_u(t) = \frac{u(t)}{\sqrt{1 - \int_0^t dt' \,|u(t')|^2}}
```

emits the initial intracavity state into the traveling mode $u(t)$. Conversely, for a chosen output mode
``v(t)``, the coupling

```math
g_v(t) = -\frac{v^*(t)}{\sqrt{\int_0^t dt' \,|v(t')|^2}}
```

absorbs that mode into the virtual output cavity. These expressions make it possible to treat input and output pulses as **single oscillator modes** in a cascaded network. Note, that this formalism can be extended to multiple input and output modes with modified couplings for each additional mode, see e.g. [A. Kiilerich, et al., Phys. Rev. A 102, 023717 (2020)](https://doi.org/10.1103/PhysRevA.102.023717). 

## SLH Networks for Cascaded Pulses

In the SLH formalism, each component is specified by a triple ``(S, L, H)``:

- ``S``cattering matrix
- ``L``Lindblad operators
- ``H``amiltonian

Networks are built using composition rules:

1. **Cascade** (feed-forward): $G_1 \triangleright G_2$ connects the output of $G_1$ into the input of $G_2$.
2. **Concatenation**: $G_1 \boxplus G_2$ stacks channels side-by-side.

By modeling the input and output pulses as virtual cavities and cascading them with the physical system, we
obtain an effective SLH triple for the full problem. This describes a master equation involving the system and and of auxiliary modes, which can be solved with standard Lindblad solvers. 

## Output Modes and the Correlation Function

The output field is generally multimode. To determine a **basis of temporal output modes** and their
occupations, we compute the first-order correlation function

```math
g^{(1)}(t_1, t_2) = \langle \hat L_s^\dagger(t_1)\, \hat L_s(t_2) \rangle,
```

where $\hat L_s$ is the output operator (e.g., $\hat L_s = \hat a$ for a cavity).
The eigen-decomposition

```math
g^{(1)}(t_1, t_2) = \sum_i n_i\, v_i^*(t_1)\, v_i(t_2)
```

defines an orthonormal set of temporal modes $v_i(t)$ with mean occupations $n_i$. The most populated
output modes correspond to the largest eigenvalue and are used as the physically relevant pulse shapes for
further analysis or for feedback into a second virtual cavity.

This procedure yields a principled way to extract a small number of relevant modes from a continuum and is
central to the workflows implemented in **QuantumInputOutput.jl**.

### References

The above summary follows the virtual-cavity and SLH-based pulse formalism developed in the input-output
literature on quantum pulses and cascaded systems, in particular the works by Kiilerich, Mølmer, and collaborators,
as well as subsequent extensions by Christiansen and coworkers.

- Input-output theory with quantum pulses
  - [A. Kiilerich, et. al., Phys. Rev. Lett. 123, 123604 (2019)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.123604)
  - [A. Kiilerich, et al., Phys. Rev. A 102, 023717 (2020)](https://doi.org/10.1103/PhysRevA.102.023717)
- The SLH framework for modeling quantum inputoutput networks
  - [J. Combes, et. al. Advances in Physics: X, 2:3, 784-888 (2017)](https://doi.org/10.1080/23746149.2017.1343097) 