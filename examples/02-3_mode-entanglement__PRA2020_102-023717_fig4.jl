# # Photon Number and Mode Entanglement with a Quantum Emitter

# In this example, we reproduce Fig. 4 of
# [A. Kiilerich and K. Molmer, Phys. Rev. A 102, 023717 (2020)](https://doi.org/10.1103/PhysRevA.102.023717),
# following the model described in Sec. III.B.
# A three-level Λ emitter decays through a cavity into two dominant temporal
# output modes. We identify these modes from the field autocorrelation
# function, capture them with two virtual output cavities, and visualize the
# resulting mode populations and the final atom-mode density matrix.

using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using Plots
using LaTeXStrings
using LinearAlgebra
using DataInterpolations

#

## symbolic Hilbert space
hc = FockSpace(:c)
hs = NLevelSpace(:atom, 3)
hv1 = FockSpace(:v1)
hv2 = FockSpace(:v2)
h = hc ⊗ hs ⊗ hv1 ⊗ hv2

## symbolic operators
a = Destroy(h, :a, 1)
σ(i, j) = Transition(h, :σ, i, j, 2)
av1 = Destroy(h, :a_v1, 3)
av2 = Destroy(h, :a_v2, 4)

## symbolic parameters
@rnumbers γ g ω12
gv1, gv2 = cnumbers("g_v1 g_v2")
nothing # hide

# The localized system consists of a cavity mode coupled to a three-level
# Λ emitter. The transition |g₁⟩ ↔ |e⟩ is resonant with the cavity and
# |g₂⟩ ↔ |e⟩ is detuned by ω₁₂.

H_s = g * (a' * σ(1, 3) + a * σ(3, 1) + a' * σ(2, 3) + a * σ(3, 2)) + ω12 * σ(2, 2)
G_s = SLH(1, √(γ) * a, H_s)

G_v1 = SLH(1, gv1' * av1, 0)
G_v2 = SLH(1, gv2' * av2, 0)
G = cascade(G_s, G_v1, G_v2)

H = get_hamiltonian(G)
L = get_lindblad(G)[1]
nothing # hide

# We use the parameters quoted in the paper and initialize the emitter in the
# excited state.

γ_ = 1.0
g_ = 0.1γ_
ω12_ = 0.5γ_

T = [0:0.005:1;]*100/γ_
ΔT = T[2] - T[1]

dict_p_1 = Dict([γ, g, ω12, gv1, gv2] .=> [γ_, g_, ω12_, 0.0, 0.0])
nothing # hide

## numeric bases and operators
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

# We first solve the decay dynamics without explicit output cavities and
# determine the two dominant temporal output modes from the first-order
# correlation matrix g⁽¹⁾(t₁, t₂).

H_QO_1 = translate_qo(H, b; parameter = dict_p_1)
L_QO_1 = translate_qo(L, b; parameter = dict_p_1)

function input_output_1(t, ρ)
    J = [L_QO_1]
    return H_QO_1, J, dagger.(J)
end

ψ0 = fockstate(bc, 0) ⊗ nlevelstate(bs, 3) ⊗ fockstate(bv1, 0) ⊗ fockstate(bv2, 0)
t_1, ρt_1 = timeevolution.master_dynamic(T, ψ0, input_output_1)
nothing # hide

Ls(t) = √(γ_) * a_qo
g1_m = two_time_corr_matrix(T, ρt_1, input_output_1, Ls)

F = eigen(g1_m)
n_avg = real.(F.values) * ΔT

v1_mode = F.vectors[:, end] / √(ΔT)
v2_mode = F.vectors[:, end - 1] / √(ΔT)
n1 = n_avg[end]
n2 = n_avg[end - 1]
nothing # hide

# After identifying the two modes, we add two cascaded virtual output cavities.
# The coupling of the second cavity must be corrected for the reshaping caused by
# the first output cavity, which is done with [`v_eff`](@ref).

gv1_t = v_to_gv(v1_mode, T)
v2_eff = v_eff([v1_mode, v2_mode], T, 2) 
gv2_t = v_to_gv(v2_eff, T)

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

# We monitor the excited-state population, the cavity population, and the
# populations transferred to the two output modes.

P_e_t = real.(expect(σ_qo(3, 3), ρt_2))
n_cavity_t = real.(expect(a_qo' * a_qo, ρt_2))
n1_t = real.(expect(av1_qo' * av1_qo, ρt_2))
n2_t = real.(expect(av2_qo' * av2_qo, ρt_2))
nothing # hide

# Finally, we inspect the reduced density matrix of the emitter and the two
# captured modes. The paper displays the real part of this matrix in the basis
# |g₁/g₂, n₁, n₂⟩ with n₁, n₂ ∈ {0, 1}.

ρ_atom_modes = ptrace(ρt_2[end], 1)
ρ_plot = real.(ρ_atom_modes.data[1:8, 1:8])
ρ_abs_max = maximum(abs, ρ_plot)

labels = [
    L"|g_1,0,0\rangle",
    L"|g_1,0,1\rangle",
    L"|g_1,1,0\rangle",
    L"|g_1,1,1\rangle",
    L"|g_2,0,0\rangle",
    L"|g_2,0,1\rangle",
    L"|g_2,1,0\rangle",
    L"|g_2,1,1\rangle",
]

function hinton_plot(ρ, labels)
    n = size(ρ, 1)
    p = plot(
        xlims = (0.5, n + 0.5),
        ylims = (n + 0.5, 0.5),
        xticks = (1:n, labels),
        yticks = (1:n, labels),
        xrotation = 90,
        aspect_ratio = 1,
        framestyle = :box,
        legend = false,
        grid = true,
        colorbar = true,
        clims = (-ρ_abs_max, ρ_abs_max),
    )

    for i = 1:n, j = 1:n
        if abs(ρ[i, j]) > 1e-6
            scatter!(
                p,
                [j],
                [i];
                marker = :rect,
                markersize = 30 * sqrt(abs(ρ[i, j]) / ρ_abs_max),
                markerstrokewidth = 0,
                marker_z = [ρ[i, j]],
                c = cgrad([:red, :white, :blue]),
                label = "",
            )
        end
    end

    return p
end

#

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
plot!(
    p_b,
    T,
    real.(v2_mode);
    color = :red,
    lw = 2,
    label = "n₂ = $(round(n2; digits = 2))",
)
plot!(
    p_b;
    xlabel = L"\gamma t",
    ylabel = L"\Re[v_i(t)]",
    legend = :topright,
)

p_c = plot(
    T,
    n_cavity_t;
    color = :green,
    ls = :dot,
    lw = 2,
    label = L"n_\mathrm{cavity}",
)
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

p_d = hinton_plot(ρ_plot, labels) # TODO: adjust aspect ratio, label size, size
plot!(p_d)

plot(p_a, p_b, p_c, p_d; layout = (1, 4), size = (1400, 320))

# ## Package versions

# These results were obtained using the following versions:

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
