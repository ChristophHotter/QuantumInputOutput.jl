# Implementation

This section explains how the package moves from a symbolic SLH model to a numerical time evolution, and how the pulse-specific tools fit into that pipeline. 

## Symbolic expressions

The symbolic layer is built on [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl). You define Hilbert spaces and operators explicitly and keep the SLH model usualy analytic. A typical setup looks like this:

```julia
hu = FockSpace(:u)
hs = NLevelSpace(:s, 2)
hv = FockSpace(:v)
h = hu ⊗ hs ⊗ hv

au = Destroy(h, :a_u, 1)
σ = Transition(h, :σ, 1, 2, 2)
av = Destroy(h, :a_v, 3)

γ, gu, gv = rnumbers("γ g_u g_v")
```

An SLH component is represented as `(S, L, H)` by the [`SLH`](@ref) type. The cascade [`▷`](@ref) and concatenation [`⊞`](@ref) rules implement the standard network composition from the SLH framework (Combes et al., 2017). The resulting effective operators are accessed by [`get_hamiltonian`](@ref) and [`get_lindblad`](@ref) and remain symbolic until translation. This is especially useful when you want to further manipulate the expressions, e.g. to transform into the interaction picture. 

If you directly want to use [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) operators and functions, the [`SLHqo`](@ref) type skips the symbolic layer and allows time-dependent `L` or `H` as callables while still using the same cascade and concatenate rules.

## Translate to numerics

To solve the dynamics of the master equation we first need to create the corresponding numeric operators for the Hamilton and the Lindblad terms.  
This can be done with [`translate`](@ref), which converts symbolic operators into [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) operators on a chosen basis. It accepts two parameter substitution dictionaries:

- `parameter`: numeric parameters used in algebraic substitution
- `time_parameter`: time-dependent parameters, given as functions of `t` 

If `time_parameter` is non-empty, [`translate`](@ref) returns a callable `t -> op(t)` so that the Hamiltonian and jump operators can be supplied to `timeevolution.master_dynamic`. 

In some cases it can be useful to define your own set of numeric operators which should replace the symbolic expressions, e.g. to reduce the Hilbert space if the output cavities are not analyzed but they are already included in the symbolic derivation. Such a list of operators can be provide with the dictionary `operators`.

If the dynamics of the system should be solved with a higher-order mean-field approximation, the symbolic Hamiltonian and Lindblad terms can be directly used in [QuantumCumulants.jl](https://github.com/qojulia/QuantumCumulants.jl)


## Field coupling terms

Quantum pulses are encoded through virtual cavities. Given a normalized input mode `u(t)` and output mode `v(t)`, the package constructs time-dependent couplings

``g_u(t) = u^*(t) / \sqrt{1 - \int_0^t |u(t')|^2 dt'}`` and
``g_v(t) = -v^*(t) / \sqrt{\int_0^t |v(t')|^2 dt'}``. 

The implementation uses cumulative numerical integration on a time grid `T` and a linear interpolation. 

- [`u_to_gu`](@ref) and [`v_to_gv`](@ref) build interpolated couplings from sampled modes
- [`u_to_gu_Gauss`](@ref) and [`v_to_gv_Gauss`](@ref) provide analytic expressions for Gaussian pulses 

For multiple input/output modes the distortion of the pulse due to the subsequent/preceding virtual cavities needs to be taken into account. 

For multiple pulses, the effective input mode ``u_i^{eff}(t)`` and output mode ``v_i^{eff}(t)`` for the virtual cavity `i` are constructed via [`u_eff`](@ref) and [`v_eff`](@ref).  
TODO: explain more; show equation; 


## Output modes and the correlation function

The dominant output modes are extracted by computing the two-time correlation matrix

``g^{(1)}(t_1, t_2) = \langle L_s^\dagger(t_1) L_s(t_2) \rangle``

and diagonalizing it. In this package, [`two_time_corr_matrix`](@ref) builds that matrix from a previously computed trajectory $\rho(t)$ and a chosen output operator $L_s(t)$, using the quantum regression theorem. This means, for each time point $t_1$ we calculate $L_s(t_1) \rho(t_1)$ and use this as the initial "state" for the propagation of $t_2$, with the same Hamiltonian and Lindblad terms. 

The eigenvectors of the matrix $g^{(1)}(t_1, t_2)$ correspond to temporal modes and the eigenvalues to their mean photon-number weights. The full procedure is illustrated in the [Tutorial](@ref). 

## Interaction picture

For networks with multiple virtual modes, an interaction-picture transformation can simplify the dynamics by removing the pure mode-transfer dynamics between the virtual input and output cavities. 

For the interaction picture of $N$ modes, corresponding to the Hamiltonian obtained from $G_1 \triangleright ...  G_N$, where $G_i = (0, g_i \hat a_i, 0)$, we obtain the transformation for the operators 
$(\hat a_1(t),  ..., \hat a_N(t))^T = M(t)\,(\hat a_1(0), ..., \hat a_N(0))^T$, where the coefficient matrix $N(t)$ is determined by ``\dot M(t) = A(t) M(t)``, with ``M(t_0) = I`` and the coupling matrix

```math
A(t) = \frac{1}{2}\begin{bmatrix}
0 & g_1 g_2^* & \cdots & g_1 g_N^* \\
-g_1^* g_2 & 0 & \ddots & \vdots \\
\vdots & \ddots & 0 & g_3 g_4^* \\
-g_1^* g_N & \cdots & -g_{N-1}^* g_N & 0
\end{bmatrix}
```

The functions [`interaction_picture_A_2modes`](@ref), [`interaction_picture_A_3modes`](@ref), and [`interaction_picture_A_4modes`](@ref) build the coupling matrices for two, three and four modes. For two equal modes ($u(t) = v(t)$) the analytic solution of $M(t)$ is provided by [`interaction_picture_M_2modes_equal`](@ref).

## Pulse delay

Pulse propagation delays are modeled by a virtual delay cavity that absorbs a pulse `v(t)` while emitting a target pulse `u(t)` with a controlled delay. The functions

- [`uv_to_gin`](@ref)
- [`uv_to_gout`](@ref)

compute the in-coupling and out-coupling strengths `g_in(t)` and `g_out(t)` for this delay cavity, according to

```math
\tilde g_{\mathrm{out},u,v}(t) = \frac{u^*(t)}{\sqrt{\int_0^t dt' \,|v(t')|^2 - \int_0^t dt' \,|u(t')|^2}},
\qquad
\tilde g_{\mathrm{in},v,u}(t) = -\frac{v^*(t)}{\sqrt{\int_0^t dt' \,|v(t')|^2 - \int_0^t dt' \,|u(t')|^2}}.
```

The implementation mirrors [`u_to_gu`](@ref) and [`v_to_gv`](@ref): compute cumulative integrals on a time grid and return interpolated time-dependent couplings. 

A simple single-pulse case is demonstrated in the example [Simple Pulse Delay with a Virtual Cavity](examples/08-1_pulse-delay__simple.md), where an input pulse is emitted, delayed by a virtual cavity, and captured into an output cavity. 