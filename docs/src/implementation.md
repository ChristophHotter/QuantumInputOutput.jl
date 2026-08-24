# Implementation

This section explains how the package moves from a symbolic SLH model to a numerical time evolution, and how the pulse-specific tools fit into that pipeline. 

## Symbolic expressions

The symbolic layer is built on [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl). You define Hilbert spaces and operators explicitly and usually derive the SLH model. A typical setup looks like this:

```julia
hu = FockSpace(:u)
hs = NLevelSpace(:s, 2)
hv = FockSpace(:v)
h = hu ⊗ hs ⊗ hv

au = Destroy(h, :a_u, 1)
σ = Transition(h, :σ, 1, 2, 2)
av = Destroy(h, :a_v, 3)

@variables γ::Real g_u::Complex g_v::Complex
```

An SLH component is represented as `(S, J, H)` by the [`SLH`](@ref) type. The cascade [`▷`](@ref), concatenation [`⊞`](@ref), and feedback reduction [`feedback`](@ref) rules implement the standard network composition from the SLH framework.

The resulting effective operators are accessed by [`hamiltonian`](@ref) and [`jump_operator`](@ref) and remain symbolic until translation. This is especially useful when you want to further manipulate the expressions, e.g. to transform into the interaction picture.

```julia
G_cas = ▷(G_u, G_s, G_v)

H = hamiltonian(G_cas)
J = jump_operator(G_cas)
```

For networks with internal loops, the symbolic model can be reduced directly with [`feedback`](@ref), which applies the SLH feedback reduction rule before translation. This keeps the symbolic workflow consistent: build a network from cascades and concatenations, eliminate internal connections symbolically, and only then translate the reduced Hamiltonian and Lindblad operators to numerics.

If you directly want to use [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) operators and functions, you can pass numeric operators and callables directly to [`SLH`](@ref), which supports both symbolic and numeric operator types as well as time-dependent `J` or `H` while still using the same cascade, concatenate, and feedback rules. This can be much faster. 

## Translate to numerics

To solve the dynamics of the master equation we first need to create the corresponding numeric operators for the Hamilton and the Lindblad terms.  
This can be done with [`to_numeric`](@ref), which converts symbolic operators into [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) operators on a chosen basis. It accepts two parameter substitution dictionaries:

- `parameter`: numeric parameters used in algebraic substitution
- `time_parameter`: time-dependent parameters, given as functions of `t` 

If `time_parameter` is non-empty, [`to_numeric`](@ref) returns a callable `t -> op(t)` so that the Hamiltonian and jump operators can be supplied to `timeevolution.master_dynamic`. 

```julia
bu = FockBasis(2)
bs = NLevelBasis(2)
bv = FockBasis(2)
b = bu ⊗ bs ⊗ bv

dict_p = Dict(γ => 1.0, g_v => 0.0)
gu_t = coupling_input(t -> exp(-t^2), 0:0.01:5)
dict_p_t = Dict(g_u => gu_t)

H_QO = to_numeric(H, b; parameter=dict_p, time_parameter=dict_p_t)
J_QO = to_numeric(J, b; parameter=dict_p, time_parameter=dict_p_t)
```

In some cases it can be useful to define your own set of numeric operators which should replace the symbolic expressions, e.g. to reduce the Hilbert space if the output cavities are not analyzed but they are already included in the symbolic derivation. Such a list of operators can be provide with the dictionary `operators`.

If the dynamics of the system should be solved with a higher-order mean-field approximation, the symbolic Hamiltonian and Lindblad terms can be directly used in [QuantumCumulants.jl](https://github.com/qojulia/QuantumCumulants.jl)

## Field coupling terms

Quantum pulses are encoded through virtual cavities. Given a normalized input mode `u(t)` and output mode `v(t)`, the package constructs time-dependent couplings

``g_u(t) = u^*(t) / \sqrt{1 - \int_0^t |u(t')|^2 dt'}`` and
``g_v(t) = -v^*(t) / \sqrt{\int_0^t |v(t')|^2 dt'}``. 

The implementation uses cumulative numerical integration on a time grid `T` and a linear interpolation. 

- [`coupling_input`](@ref) and [`coupling_output`](@ref) build interpolated couplings from sampled modes
- For Gaussian pulses, pass a [`Gaussian`](@ref) descriptor to [`coupling_input`](@ref) or [`coupling_output`](@ref) for analytic expressions

```julia
T = [0:0.002:1;]*12
σ = 1.0
τ = 4.0
u(t) = 1/(√(σ)*π^(1/4)) * exp(-0.5 * ((t - τ) / σ)^2)

gu_t = coupling_input(u, T)
gv_t = coupling_output(u, T)
```

For multiple input/output modes the distortion of the pulse due to the subsequent/preceding virtual cavities needs to be taken into account. The effective input mode ``u_i^{\mathrm{eff}}(t)`` and output mode ``v_i^{\mathrm{eff}}(t)`` for the virtual cavity `i` are constructed via [`effective_input_mode`](@ref) and [`effective_output_mode`](@ref). We first compute `u_i_eff` / `v_i_eff` and then pass them to [`coupling_input`](@ref) / [`coupling_output`](@ref).  

```julia
# Two input modes and two output modes on the same time grid T
u_modes = [u1_mode, u2_mode]
v_modes = [v1_mode, v2_mode]

# Effective modes for the second virtual cavities: u_2^(1), v_2^(1)
u2_eff = effective_input_mode(u_modes, T, 2)   # = u_2^eff = u_2^(2-1)
v2_eff = effective_output_mode(v_modes, T, 2)  # = v_2^eff = v_2^(2-1)

# Couplings from effective modes
gu2_eff = coupling_input(u2_eff, T)
gv2_eff = coupling_output(v2_eff, T)
```

## Output modes and the correlation function

The dominant output modes are extracted by computing the two-time correlation matrix

```math
g^{(1)}(t_1, t_2) = \langle J_s^\dagger(t_1) J_s(t_2) \rangle
```

and diagonalizing it. In this package, [`correlation_matrix`](@ref) builds that matrix from a previously computed trajectory $\rho(t)$ and a chosen output operator $J_s(t)$, using the quantum regression theorem. This means, for each time point $t_1$ we calculate $J_s(t_1) \rho(t_1)$ and use this as the initial "state" for the propagation of $t_2$, with the same Hamiltonian and Lindblad terms. 

The eigenvectors of the matrix $g^{(1)}(t_1, t_2)$ correspond to temporal modes and the eigenvalues to their mean photon-number weights. The full procedure is illustrated in the [Tutorial](@ref). 

```julia
Js(t) = gu_t(t) * au_qo + √(1.0) * c_qo
g1 = correlation_matrix(T, ρt, input_output_1, Js)
F = eigen(g1)
v_mode = F.vectors[:, end] / sqrt(T[2] - T[1])
```

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

The function [`coupling_matrix`](@ref) builds the coupling matrix for an arbitrary number of modes. The ODE for the coefficient matrix $M(t)$ is solved by [`solve_mode_evolution`](@ref). For two equal modes ($u(t) = v(t)$) the analytic solution of $M(t)$ is provided by [`solve_mode_evolution_symmetric`](@ref).

Using the interaction picture for scattering with a Fock state $| n=20 \rangle$ on a two-level system, is  demonstrated in the example [Interaction Picture Scattering with a Quantum Pulse](examples/06-1_interaction-picture__PRA2023_107-013706_fig2.md). 

## Pulse delay

Pulse propagation delays are modeled by a virtual delay cavity that absorbs a pulse `v(t)` while emitting a target pulse `u(t)` with a controlled delay. The functions

- [`coupling_delay_in`](@ref)
- [`coupling_delay_out`](@ref)

compute the in-coupling and out-coupling strengths `g_in(t)` and `g_out(t)` for this delay cavity, according to

```math
\tilde g_{\mathrm{out},u,v}(t) = \frac{u^*(t)}{\sqrt{\int_0^t dt' \,|v(t')|^2 - \int_0^t dt' \,|u(t')|^2}},
\qquad
\tilde g_{\mathrm{in},v,u}(t) = -\frac{v^*(t)}{\sqrt{\int_0^t dt' \,|v(t')|^2 - \int_0^t dt' \,|u(t')|^2}}.
```

The implementation mirrors [`coupling_input`](@ref) and [`coupling_output`](@ref): compute cumulative integrals on a time grid and return interpolated time-dependent couplings. 

A simple single-pulse case is demonstrated in the example [Simple Pulse Delay with a Virtual Cavity](examples/08-1_pulse-delay__simple.md), where an input pulse is emitted, delayed by a virtual cavity, and captured into an output cavity. 
