# Theory

This section summarizes the theoretical ideas behind **QuantumInputOutput.jl**, i.e. input-output theory for quantum pulses, its formulation via SLH networks, and the use of virtual cavities to handle traveling wave packets. We also explain how temporal output-mode bases are determined using two-time correlation functions and describe the interaction-picture, as well as pulse delay for traveling pulses. 

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

The key idea is to replace the traveling input and output pulses by **virtual cavities**, each with a time-dependent, complex coupling to the system determined by the corresponding temporal mode.  

For a normalized input mode $u(t)$, the coupling

```math
g_u(t) = \frac{u(t)}{\sqrt{1 - \int_0^t dt' \,|u(t')|^2}}
```

emits the initial intracavity state into the traveling mode $u(t)$. Conversely, for a chosen output mode
$v(t)$, the coupling

```math
g_v(t) = -\frac{v^*(t)}{\sqrt{\int_0^t dt' \,|v(t')|^2}}
```

absorbs that mode into the virtual output cavity. These expressions make it possible to treat input and output pulses as **single oscillator modes** in a cascaded network. 

For multiple input and output modes, one first constructs effective modes
``u_i^{\mathrm{eff}}(t)`` and ``v_i^{\mathrm{eff}}(t)`` (which include distortions from
the other virtual cavities). In recursive form (for ``j<i``):

```math
\dot\alpha_{v_i}^{(j)}(t) =
-g_{v_j}(t)\,v_i^{(j-1)}(t)
+ \frac{1}{2}|g_{v_j}(t)|^2\alpha_{v_i}^{(j)}(t),
\qquad
v_i^{(j)}(t) = v_i^{(j-1)}(t) + g_{v_j}^*(t)\alpha_{v_i}^{(j)}(t),
```

```math
\dot\alpha_{u_i}^{(j)}(t) =
-g_{u_j}(t)\,u_i^{(j-1)}(t)
+ \frac{1}{2}|g_{u_j}(t)|^2\alpha_{u_i}^{(j)}(t),
\qquad
u_i^{(j)}(t) = u_i^{(j-1)}(t) - g_{u_j}^*(t)\alpha_{u_i}^{(j)}(t).
```

With ``u_i^{(0)}(t)=u_i(t)`` and ``v_i^{(0)}(t)=v_i(t)``, the effective modes are
``u_i^{\mathrm{eff}}(t)=u_i^{(i-1)}(t)`` and ``v_i^{\mathrm{eff}}(t)=v_i^{(i-1)}(t)``.
The corresponding couplings then follow from the single-mode expressions:

```math
g_{u_i}^{\mathrm{eff}}(t) =
\frac{\big(u_i^{\mathrm{eff}}(t)\big)^*}{\sqrt{1 - \int_0^t dt' \,|u_i^{\mathrm{eff}}(t')|^2}},
\qquad
g_{v_i}^{\mathrm{eff}}(t) =
-\frac{\big(v_i^{\mathrm{eff}}(t)\big)^*}{\sqrt{\int_0^t dt' \,|v_i^{\mathrm{eff}}(t')|^2}}.
```

This is the multi-mode extension used by [`effective_input_mode`](@ref) and
[`effective_output_mode`](@ref), see e.g.
[A. Kiilerich, et al., Phys. Rev. A 102, 023717 (2020)](https://doi.org/10.1103/PhysRevA.102.023717).

## SLH Networks for Cascaded Pulses

In the SLH formalism, each component is specified by a triple `(S, J, H)`:

- `S`cattering matrix
- `J`ump (Lindblad) operators
- `H`amiltonian

Networks are built using composition rules:

1. **Cascade**: $G_1 \triangleright G_2$ connects the output of $G_1$ into the input of $G_2$ 

```math
G_1 \triangleright G_2 = \left(
S_2 S_1,\;
J_2 + S_2 J_1,\;
H_1 + H_2 + \frac{1}{2i}\left(J_2^\dagger S_2 J_1 - J_1^\dagger S_2^\dagger J_2\right)
\right).
```

2. **Concatenation**: $G_1 \boxplus G_2$ stacks channels side-by-side 

```math
G_1 \boxplus G_2 = \left(
\begin{bmatrix}
S_1 & 0 \\
0 & S_2
\end{bmatrix},\;
\begin{bmatrix}
J_1 \\
J_2
\end{bmatrix},\;
H_1 + H_2
\right).
```

3. **Feedback reduction**: connecting output port $x$ back into input port $y$ of an $n$-port component
$G = (S, J, H)$ yields a reduced $(n-1)$-port model $[(G)]_{x \to y} = (\tilde S, \tilde J, \tilde H)$ with

```math
\tilde S = S_{\bar x,\bar y} + S_{\bar x,y}(1 - S_{x,y})^{-1}S_{x,\bar y},
```

```math
\tilde J = J_{\bar x} + S_{\bar x,y}(1 - S_{x,y})^{-1}J_x,
```

```math
\tilde H = H + \frac{1}{2i}\left[
\left(\sum_{j=1}^n J_j^\dagger S_{j,y}\right)(1 - S_{x,y})^{-1}J_x - \mathrm{h.c.}
\right].
```

Here $\bar x$ and $\bar y$ denote the remaining output and input ports after removing $x$ and $y$. This rule closes internal loops in an SLH network and is used whenever a network contains a direct coherent feedback path.

By modeling the input and output pulses as virtual cavities and cascading them with the physical system, we
obtain an effective SLH triple for the full problem. This describes a master equation involving the system and auxiliary modes, which can be solved with standard Lindblad solvers. 

## Common SLH Elements

Appendix 1 of [Combes et al. (2017)](https://doi.org/10.1080/23746149.2017.1343097) lists SLH triples for commonly used components. The most frequently used elements in this package are summarized below using the standard ``(S, J, H)`` ordering.

Phase shifter (single input/output):

```math
G_{\mathrm{PS}} = \left(e^{i\phi},\, 0,\, 0\right).
```

Beam splitter (two inputs/outputs, with unitary $S$):

```math
G_{\mathrm{BS}} = \left(
\begin{bmatrix}
r_{11} & t_{12} \\
t_{21} & r_{22}
\end{bmatrix},\,
0,\,
0
\right),
\quad S^\dagger S = I.
```

Coherent drive (displaces the input field by $\alpha(t)$):

```math
G_{\mathrm{coh}} = \left(1,\, \alpha(t),\, 0\right).
```

One-sided cavity (decay rate $\kappa$, detuning $\Delta_c$, mode operator $a$):

```math
G_{\mathrm{1s}} = \left(I,\, \sqrt{\kappa}\,a,\, \Delta_c\,a^\dagger a\right).
```

Two-sided cavity / Fabry–Perot (decay rates $\kappa_1,\kappa_2$):

```math
G_{\mathrm{2s}} = \left(
I_2,\,
\begin{bmatrix}
\sqrt{\kappa_1}\,a \\
\sqrt{\kappa_2}\,a
\end{bmatrix},\,
\Delta_c\,a^\dagger a
\right).
```

Two-level system side-coupled to a waveguide (guided decay $\kappa_g$ and unguided decay $\kappa_\perp$):

```math
G_{\mathrm{TLS}} = \left(
I_2,\,
\begin{bmatrix}
\sqrt{\kappa_g}\,\sigma_- \\
\sqrt{\kappa_\perp}\,\sigma_-
\end{bmatrix},\,
\frac{\omega}{2}\sigma_z
\right).
```

## Output Modes and the Correlation Function

The output field is generally multimode. To determine a proper **basis of temporal output modes** and their
occupations, we compute the first-order correlation function

```math
g^{(1)}(t_1, t_2) = \langle \hat J_s^\dagger(t_1)\, \hat J_s(t_2) \rangle,
```

where $\hat J_s$ is the output operator (e.g., $\hat J_s = g_u(t) \hat{a}_u + \sqrt{\kappa} \hat c$ for a single sided cavity interacting with a pulse$u(t)$).
The eigen-decomposition

```math
g^{(1)}(t_1, t_2) = \sum_i n_i\, v_i^*(t_1)\, v_i(t_2)
```

defines an orthonormal set of temporal modes $v_i(t)$ with mean occupations $n_i$. The most populated
output modes correspond to the largest eigenvalue and are used as the physically relevant pulse shapes for
further analysis. This procedure yields a principled way to extract a small number of relevant modes from a continuum and is
central to the workflows implemented in **QuantumInputOutput.jl**.

## Interaction Picture for Traveling Pulses

For pulse scattering problems (with large photon numbers) it is often advantageous to transform to an interaction picture that removes the pure mode-transfer dynamics between the virtual input and output cavities. Following [Christiansen et al. (2023)](https://doi.org/10.1103/PhysRevA.107.013706), 
one splits the total Hamiltonian into the system part and a coupling part that only exchanges excitations between
the upstream and downstream virtual modes. For the two-mode case this is

```math
\hat H_{uv}(t) = \frac{i}{2}\left[g_u(t) g_v^*(t)\,\hat a_u^\dagger \hat a_v - g_v(t) g_u^*(t)\,\hat a_v^\dagger \hat a_u\right].
```

The corresponding interaction-picture transformation affects only the field-mode operators. Writing
$(\hat a_u(t), \hat a_v(t))^T = M(t)\,(\hat a_u(0), \hat a_v(0))^T$, the coefficient matrix solves

```math
\frac{d}{dt} M(t) = A(t)\,M(t), \qquad
A(t) = \frac{1}{2}\begin{bmatrix}
0 & g_u(t) g_v^*(t) \\
-g_u^*(t) g_v(t) & 0
\end{bmatrix},
```

with $M(0)=I$. In this picture the remaining Hamiltonian and dissipators are obtained by substituting the
interaction-picture mode operators into the full cascaded model and subtracting $\hat H_{uv}(t)$ from the original
Hamiltonian. The key practical benefit is that the virtual modes no longer undergo full emptying and refilling
during the pulse; instead, only a narrow band of Fock states is populated at any time, which allows much smaller
Hilbert-space truncations for numerical simulations. 

The interaction picture can be generalized to multiple cascaded modes with

```math
A(t) = \frac{1}{2}\begin{bmatrix}
0 & g_1 g_2^* & \cdots & g_1 g_N^* \\
-g_1^* g_2 & 0 & \ddots & \vdots \\
\vdots & \ddots & 0 & g_3 g_4^* \\
-g_1^* g_N & \cdots & -g_{N-1}^* g_N & 0
\end{bmatrix}
```

which corresponds to the interaction of the Hamiltonian obtained from $G_1 \triangleright ... \triangleright G_N$, where $G_i = (0, g_i \hat a_i, 0)$. 

## Pulse Delay

Propagation delays become essential when pulses traverse different path lengths as for example in a Mach-Zehnder interferometer. Such delays can be modeled without discretizing the entire field continuum by introducing a **virtual delay cavity** that first absorbs a chosen pulse and later re-emits it with a controlled delay.
For long delays, the pulse can first be captured and then re-emitted. For short delays, however, capture and emission must occur simultaneously. In that case, the couplings are modified to account for both incoming and outgoing photons, yielding

```math
\tilde g_{\mathrm{out},u,v}(t) = \frac{u^*(t)}{\sqrt{\int_0^t dt' \,|v(t')|^2 - \int_0^t dt' \,|u(t')|^2}}
\qquad
\tilde g_{\mathrm{in},v,u}(t) = -\frac{v^*(t)}{\sqrt{\int_0^t dt' \,|v(t')|^2 - \int_0^t dt' \,|u(t')|^2}}
```

The same formalism extends to interferometric settings by introducing one delay cavity per path and adjusting only the relative delays, which usually is the physically relevant quantity. 

## References

- Input-output theory with quantum pulses
  - [A. Kiilerich, et al., Phys. Rev. Lett. 123, 123604 (2019)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.123604)
  - [A. Kiilerich, et al., Phys. Rev. A 102, 023717 (2020)](https://doi.org/10.1103/PhysRevA.102.023717)
- The SLH framework for modeling quantum inputoutput networks
  - [J. Combes, et al. Advances in Physics: X, 2:3, 784-888 (2017)](https://doi.org/10.1080/23746149.2017.1343097) 
- Interaction picture and pulse delay for quantum pulses
  - [Christiansen et al., Phys. Rev. A 107, 013706 (2023)](https://doi.org/10.1103/PhysRevA.107.013706)
  - [Christiansen et al., Phys. Rev. A 113, 013730 (2026)](https://doi.org/10.1103/3f3w-jmj8)
