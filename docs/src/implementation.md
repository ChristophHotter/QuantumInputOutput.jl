# Implementation

This section explains how the package moves from a symbolic SLH model to a numerical time evolution, and how the pulse-specific tools fit into that pipeline. 

## Symbolic expressions

The symbolic layer is built on `SecondQuantizedAlgebra.jl`. You define Hilbert spaces and operators explicitly and keep the SLH model usualy analytic. A typical setup looks like this:

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

An SLH component is represented as `(S, L, H)` by the `SLH` type. The cascade `▷` and concatenation `⊞` rules implement the standard network composition from the SLH framework (Combes et al., 2017). The resulting effective operators are accessed by `get_hamiltonian` and `get_lindblad` and remain symbolic until translation. This is especially useful when you want to further manipulate the expressions, e.g. to transform into the interaction picture. 

If you directly want to use `QuantumOptics.jl` operators and functions, the `SLHqo` type skips the symbolic layer and allows time-dependent `L` or `H` as callables while still using the same cascade and concatenate rules.

## Translate to numerics

To solve the dynamics of the master equation we first need to create the corresponding numeric operators for the Hamilton and the Lindblad terms.  
This can be done with `translate`, which converts symbolic operators into `QuantumOptics.jl` operators on a chosen basis. It accepts two parameter substitution dictionaries:

- `parameter`: numeric parameters used in algebraic substitution
- `time_parameter`: time-dependent parameters, given as functions of `t` 

If `time_parameter` is non-empty, `translate` returns a callable `t -> op(t)` so that the Hamiltonian and jump operators can be supplied to `timeevolution.master_dynamic`. 

In some cases it can be useful to define your own set of numeric operators which should replace the symbolic expressions, e.g. to reduce the Hilbert space if the output cavities are not analyzed but they are already included in the symbolic derivation. Such a list of operators can be provide with the dictionary `operators`.

If the dynamics of the system should be solved with a higher-order mean-field approximation, the symbolic Hamiltonian and Lindblad terms can be directly used in `QuantumCumulants.jl`


## Field coupling terms

Quantum pulses are encoded through virtual cavities. Given a normalized input mode `u(t)` and output mode `v(t)`, the package constructs time-dependent couplings

``g_u(t) = u^*(t) / \sqrt{1 - \int_0^t |u(t')|^2 dt'}`` and
``g_v(t) = -v^*(t) / \sqrt{\int_0^t |v(t')|^2 dt'}``.

The implementation uses cumulative numerical integration on a time grid `T` and a linear interpolation. 

- `u_to_gu(u, T)` and `v_to_gv(v, T)` build interpolated couplings from sampled modes
- `u_to_gu_Gauss(τ, σ; δ=0)` and `v_to_gv_Gauss(τ, σ; δ=0)` provide analytic expressions for Gaussian pulses 

For multiple input/output modes the distortion of the pulse due to the subsequent/preceding virtual cavities needs to be taken into account. 

For multiple pulses, the effective input mode `u_i^{eff}(t)` and output mode `v_i^{eff}(t)` for the virtual cavity `i` are constructed via `u_eff` and `v_eff`.  
TODO: explain more; show equation; 


## Output modes and the correlation function

The dominant output modes are extracted by computing the two-time correlation matrix

``g^{(1)}(t_1, t_2) = \langle L_s^\dagger(t_1) L_s(t_2) \rangle``

and diagonalizing it. In this package, `two_time_corr_matrix` builds that matrix from a previously computed trajectory `ρ(t)` and a chosen output operator `L_s(t)`. The eigenvectors correspond to temporal modes and the eigenvalues to their photon-number weights, consistent with the pulse-mode treatment of Kiilerich & Mølmer (2019) and the cascaded master-equation construction in Christiansen et al. (2023).

The full end-to-end procedure is illustrated in the `Tutorial` and in the cavity-scattering example `examples/01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3.md`.

## Interaction picture

For networks with multiple virtual modes, an interaction-picture transformation can simplify the dynamics by factoring out known mode mixing. The functions in `src/interaction_picture.jl` implement

``\dot M(t) = A(t) M(t)``, with ``M(t_0) = I``,

where `A(t)` depends on the couplings `g_u(t), g_v(t), ...`. The numerical solver is exposed as `interaction_picture_M`, and analytic shortcuts are provided in special cases, e.g. `interaction_picture_M_2modes_equal` for two equal pulses (Christiansen et al., 2023, PRA 107, 013706).

The `interaction_picture_A_2modes`, `interaction_picture_A_3modes`, and `interaction_picture_A_4modes` helpers build the coupling matrices used in those transformations. These tools are used in the interaction-picture examples and are the preferred route when you want stable dynamics for long pulses or weak couplings.

## Pulse delay

Pulse propagation delays are modeled by a virtual delay cavity that absorbs a pulse `v(t)` while emitting a target pulse `u(t)` with a controlled delay. The functions

- `uv_to_gin(u, v, T)`
- `uv_to_gout(u, v, T)`

compute the in-coupling and out-coupling strengths `g_in(t)` and `g_out(t)` for this delay cavity. This construction follows the propagation-delay treatment in Christiansen & Mølmer (2026) and builds on the pulse-shaping ideas of Kiilerich & Mølmer (2019).

For a concrete implementation, see the pulse-delay example `examples/08-1_pulse-delay__simple.md` and the related draft in `.other/08-2_pulse-propagation-delay__PRA2026_113-013730_fig4.jl`.
