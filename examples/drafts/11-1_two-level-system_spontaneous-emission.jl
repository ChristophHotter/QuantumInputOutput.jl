# # Dominant Output Mode of a Decaying Two-Level System

# A two-level system initially in its excited state emits one photon by spontaneous
# emission. We determine the dominant temporal output mode from the first-order field
# correlation function, capture that mode in a virtual output cavity, and calculate the
# overlap of its final state with the one-photon Fock state.

using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using QuantumOpticsBase: dagger
using DataInterpolations
using LinearAlgebra
using Plots
using LaTeXStrings

#
## symbolic Hilbert space
hs = NLevelSpace(:tls, 2)
hv = FockSpace(:v)
h = hs ⊗ hv

## symbolic operators
σ(i, j) = Transition(h, :σ, i, j, 1)
av = Destroy(h, :a_v, 2)

## symbolic parameters
@variables γ::Real g_v::Number

# The emitter couples to a single Markovian output channel. Cascading a virtual
# cavity after it allows a chosen temporal mode to be represented as a discrete
# oscillator mode.

G_s = SLH(1, √(γ) * σ(1, 2), 0)
G_v = SLH(1, g_v' * av, 0)
G = cascade(G_s, G_v)

H = hamiltonian(G)
L = lindblad(G)[1]
nothing # hide

# We use the spontaneous-emission rate as the frequency unit. At
# `T_end = 10 / γ`, the residual excited-state population is `exp(-10)`.

γ_ = 1.0
T_end = 10 / γ_
T = collect(0.0:0.01:1)*T_end
ΔT = T[2] - T[1]

## numerical bases and operators
bs = NLevelBasis(2)
bv = FockBasis(1)
b = bs ⊗ bv

σm_qo = transition(bs, 1, 2) ⊗ one(bv)
Pe_qo = transition(bs, 2, 2) ⊗ one(bv)
av_qo = one(bs) ⊗ destroy(bv)
nv_qo = av_qo' * av_qo

ψ0 = nlevelstate(bs, 2) ⊗ fockstate(bv, 0)
nothing # hide

# ## Extract the dominant temporal output mode

# In the first evolution the virtual cavity is disconnected. The first-order
# correlation matrix of the emitted field is
#
# ```math
# g^{(1)}(t_1,t_2) =
# \langle L_s^\dagger(t_1)L_s(t_2)\rangle,
# \qquad L_s = \sqrt{\gamma}\,\sigma_-.
# ```

dict_p_1 = Dict(γ => γ_, g_v => 0.0)
H_QO_1 = to_numeric(H, b; parameter = dict_p_1)
L_QO_1 = to_numeric(L, b; parameter = dict_p_1)

t_1, ρt_1 = timeevolution.master(T, ψ0, H_QO_1, [L_QO_1])

Ls = √(γ_) * σm_qo
g1 = correlation_matrix(T, ρt_1, H_QO_1, [L_QO_1], Ls)

# We diagonalize the quadrature-weighted integral kernel. Its eigenvalues are
# the mean photon numbers in the corresponding temporal modes. `Hermitian`
# removes insignificant numerical asymmetry before diagonalization.

weights = fill(ΔT, length(T))
weights[[1, end]] ./= 2
√weights = sqrt.(weights)
weighted_g1 = √weights .* ((g1 + g1') / 2) .* transpose(√weights)
mode_decomposition = eigen(Hermitian(weighted_g1))
mode_occupations = real.(mode_decomposition.values)
dominant_index = argmax(mode_occupations)
n_dominant = mode_occupations[dominant_index]
v_mode = mode_decomposition.vectors[:, dominant_index] ./ √weights

# Fix the arbitrary eigenvector phase by aligning the numerical mode with the
# analytical spontaneous-emission envelope
# `v(t) = sqrt(γ) exp(-γ t / 2)`.

v_analytic = √(γ_) .* exp.(-γ_ .* T ./ 2)
v_analytic ./= √(sum(weights .* abs2.(v_analytic)))
phase_overlap = sum(weights .* conj.(v_analytic) .* v_mode)
v_mode .*= exp(-1im * angle(phase_overlap))
mode_overlap = abs(sum(weights .* conj.(v_analytic) .* v_mode))^2
nothing # hide

# ## Capture the mode and calculate its one-photon fidelity

v_mode_interpolation = LinearInterpolation(v_mode, T)
T_2 = [0:0.0001:1;]*T[end]
gv_t = coupling_output(v_mode_interpolation, T_2)

# Without the finer time grid the fidelity is much smaller. 

dict_p_2 = Dict(γ => γ_)
dict_p_t_2 = Dict(g_v => gv_t)
H_QO_2 = to_numeric(H, b; parameter = dict_p_2, time_parameter = dict_p_t_2)
L_QO_2 = to_numeric(L, b; parameter = dict_p_2, time_parameter = dict_p_t_2)

function input_output_2(t, ρ)
    Ht = H_QO_2(t)
    J = [L_QO_2(t)]
    return Ht, J, dagger.(J)
end

t_2, ρt_2 = timeevolution.master_dynamic(T, ψ0, input_output_2)

Pe_t = real.(expect(Pe_qo, ρt_2))
nv_t = real.(expect(nv_qo, ρt_2))

# Trace out the two-level system. Since the target is pure, its state fidelity
# is simply the projector overlap `Tr[ρ_v |1⟩⟨1|]`.

ρv_final = ptrace(ρt_2[end], 1)
one_photon_projector = projector(fockstate(bv, 1))
F_1 = real(expect(one_photon_projector, ρv_final))

@show n_dominant
@show mode_overlap
@show F_1
nothing # hide

#
using Plots
p_mode = plot(
    T,
    real.(v_mode);
    lw = 2,
    color = :blue,
    label = "dominant numerical mode",
    xlabel = L"\gamma t",
    ylabel = L"v(t)/\sqrt{\gamma}",
)
plot!(
    p_mode,
    T,
    v_analytic;
    lw = 2,
    ls = :dash,
    color = :black,
    label = L"\sqrt{\gamma}e^{-\gamma t/2}",
)

p_capture = plot(
    T,
    Pe_t;
    lw = 2,
    color = :red,
    label = L"P_e",
    xlabel = L"\gamma t",
    ylabel = "population",
    ylims = (0, 1.02),
)
plot!(p_capture, T, nv_t; lw = 2, color = :blue, label = L"\langle n_v\rangle")
annotate!(
    p_capture,
    0.58T_end,
    0.45,
    text("F₁ = $(round(F_1; digits = 6))", 11, :left),
)

plot(p_mode, p_capture; layout = (1, 2), size = (900, 350))

