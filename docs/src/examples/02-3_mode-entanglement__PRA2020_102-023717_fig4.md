```@meta
EditURL = "../../../examples/02-3_mode-entanglement__PRA2020_102-023717_fig4.jl"
```

# Photon Number and Mode Entanglement with a Quantum Emitter

We simulate the decay of a three-level Λ emitter through a cavity into two dominant temporal
output modes. We identify these modes from the field autocorrelation
function, capture them with two virtual output cavities, and visualize the
resulting mode populations and the final atom-mode density matrix.
This example, reproduces Fig. 4 of
[A. Kiilerich and K. Molmer, Phys. Rev. A 102, 023717 (2020)](https://doi.org/10.1103/PhysRevA.102.023717).

````@example 02-3_mode-entanglement__PRA2020_102-023717_fig4
using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using Plots
using LaTeXStrings
using LinearAlgebra
using DataInterpolations
````

````@example 02-3_mode-entanglement__PRA2020_102-023717_fig4
# symbolic Hilbert space
hc = FockSpace(:c)
hs = NLevelSpace(:atom, 3)
hv1 = FockSpace(:v1)
hv2 = FockSpace(:v2)
h = hc ⊗ hs ⊗ hv1 ⊗ hv2

# symbolic operators
a = Destroy(h, :a, 1)
σ(i, j) = Transition(h, :σ, i, j, 2)
av1 = Destroy(h, :a_v1, 3)
av2 = Destroy(h, :a_v2, 4)

# symbolic parameters
@variables γ::Real g::Real ω12::Real
gv1, gv2 = complex_vars("g_v1 g_v2")
nothing # hide
````

The localized system consists of a cavity mode coupled to a three-level
Λ emitter. The transition |g₁⟩ ↔ |e⟩ is resonant with the cavity and
|g₂⟩ ↔ |e⟩ is detuned by ω₁₂.

````@example 02-3_mode-entanglement__PRA2020_102-023717_fig4
H_s = g * (a' * σ(1, 3) + a * σ(3, 1) + a' * σ(2, 3) + a * σ(3, 2)) + ω12 * σ(2, 2)
G_s = SLH(1, √(γ) * a, H_s)

G_v1 = SLH(1, gv1' * av1, 0)
G_v2 = SLH(1, gv2' * av2, 0)
G = cascade(G_s, G_v1, G_v2)

H = hamiltonian(G)
L = lindblad(G)[1]
nothing # hide
````

We use the parameters quoted in the paper and initialize the emitter in the excited state.

````@example 02-3_mode-entanglement__PRA2020_102-023717_fig4
γ_ = 1.0
g_ = 0.1γ_
ω12_ = 0.5γ_

T = [0:0.005:1;]*100/γ_
ΔT = T[2] - T[1]

dict_p_1 = Dict([γ, g, ω12, gv1, gv2] .=> [γ_, g_, ω12_, 0.0, 0.0])
nothing # hide

# numeric bases and operators
bc = FockBasis(1)
bs = NLevelBasis(3)
bv1 = FockBasis(1)
bv2 = FockBasis(1)
b = bc ⊗ bs ⊗ bv1 ⊗ bv2

a_qo = destroy(bc) ⊗ one(bs) ⊗ one(bv1) ⊗ one(bv2)
σ_qo(i, j) = one(bc) ⊗ transition(bs, i, j) ⊗ one(bv1) ⊗ one(bv2)
av1_qo = one(bc) ⊗ one(bs) ⊗ destroy(bv1) ⊗ one(bv2)
av2_qo = one(bc) ⊗ one(bs) ⊗ one(bv1) ⊗ destroy(bv2)
nothing # hide
````

We first solve the decay dynamics without explicit output cavities and
determine the two dominant temporal output modes from the first-order
correlation matrix g⁽¹⁾(t₁, t₂).

````@example 02-3_mode-entanglement__PRA2020_102-023717_fig4
H_QO_1 = translate_qo(H, b; parameter = dict_p_1)
L_QO_1 = translate_qo(L, b; parameter = dict_p_1)

function input_output_1(t, ρ)
    J = [L_QO_1]
    return H_QO_1, J, dagger.(J)
end

ψ0 = fockstate(bc, 0) ⊗ nlevelstate(bs, 3) ⊗ fockstate(bv1, 0) ⊗ fockstate(bv2, 0)
t_1, ρt_1 = timeevolution.master_dynamic(T, ψ0, input_output_1)
nothing # hide

Ls = √(γ_) * a_qo
g1_m = correlation_matrix(T, ρt_1, input_output_1, Ls)

F = eigen(g1_m)
n_avg = real.(F.values) * ΔT

v1_mode = F.vectors[:, end] / √(ΔT)
v2_mode = F.vectors[:, end-1] / √(ΔT)
n1 = n_avg[end]
n2 = n_avg[end-1]
nothing # hide
````

After identifying the two modes, we add two cascaded virtual output cavities.
The coupling of the second cavity must be corrected for the reshaping caused by
the first output cavity, which is done with [`effective_output_mode`](@ref).

````@example 02-3_mode-entanglement__PRA2020_102-023717_fig4
gv1_t = coupling_output(v1_mode, T)
v2_eff = effective_output_mode([v1_mode, v2_mode], T, 2)
gv2_t = coupling_output(v2_eff, T)

dict_p_2 = Dict([γ, g, ω12] .=> [γ_, g_, ω12_])
dict_p_t_2 = Dict(gv1 => gv1_t, gv2 => gv2_t)

H_QO_2 = translate_qo(H, b; parameter = dict_p_2, time_parameter = dict_p_t_2)
L_QO_2 = translate_qo(L, b; parameter = dict_p_2, time_parameter = dict_p_t_2)

function input_output_2(t, ρ)
    J = [L_QO_2(t)]
    return H_QO_2(t), J, dagger.(J)
end

t_2, ρt_2 = timeevolution.master_dynamic(T, ψ0, input_output_2)
nothing # hide
````

We monitor the excited-state population, the cavity population, and the
populations transferred to the two output modes.

````@example 02-3_mode-entanglement__PRA2020_102-023717_fig4
P_e_t = real.(expect(σ_qo(3, 3), ρt_2))
n_cavity_t = real.(expect(a_qo' * a_qo, ρt_2))
n1_t = real.(expect(av1_qo' * av1_qo, ρt_2))
n2_t = real.(expect(av2_qo' * av2_qo, ρt_2))
nothing # hide
````

````@example 02-3_mode-entanglement__PRA2020_102-023717_fig4
p_a = heatmap(
    T,
    T,
    real.(g1_m);
    c = :inferno,
    xlabel = L"\gamma t_2",
    ylabel = L"\gamma t_1",
    colorbar_title = L"g^{(1)}(t_1,t_2)",
    xlims = (0, 100),
    ylims = (0, 100),
    aspect_ratio = 1,
)

p_b = plot(
    T,
    real.(v1_mode);
    color = :blue,
    ls = :dash,
    lw = 2,
    label = "n₁ = $(round(n1; digits = 2))",
)
plot!(p_b, T, real.(v2_mode); color = :red, lw = 2, label = "n₂ = $(round(n2; digits = 2))")
plot!(p_b; xlabel = L"\gamma t", ylabel = L"\Re[v_i(t)]", legend = :topright)

p_c = plot(T, n_cavity_t; color = :green, ls = :dot, lw = 2, label = L"n_\mathrm{cavity}")
plot!(p_c, T, P_e_t; color = :black, ls = :dashdot, lw = 2, label = L"P(|e\rangle)")
plot!(p_c, T, n1_t; color = :blue, lw = 2, label = L"n_1")
plot!(p_c, T, n2_t; color = :red, ls = :dash, lw = 2, label = L"n_2")
plot!(
    p_c;
    xlabel = L"\gamma t",
    ylabel = "mean excitation",
    ylims = (0, 1),
    legend = :topright,
)

plot(p_a, p_b, p_c, layout = (1, 3), size = (1200, 360))
````

Note that the plotted real part of the output-mode functions are not the same
as in the paper due the arbitrary global phase.

## Package versions

These results were obtained using the following versions:

````@example 02-3_mode-entanglement__PRA2020_102-023717_fig4
using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    [
        "QuantumInputOutput",
        "SecondQuantizedAlgebra",
        "QuantumOptics",
        "Plots",
        "DataInterpolations",
        "LaTeXStrings",
    ],
    mode = PKGMODE_MANIFEST,
)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

