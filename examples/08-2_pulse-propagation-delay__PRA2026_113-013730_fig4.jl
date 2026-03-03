# # TODO

# # Ramsey Interference with Delayed Fock Pulses

# In this example, we reproduces the Ramsey-like interference pattern for a partially delayed input Fock state with `n = 9`, studied in  [V. R. Christiansen and K. Mølmer, Phys. Rev. A 113, 013730 (2026)](https://doi.org/10.1103/PhysRevA.113.013730). 
# We model a single input Fock pulse that is split on a 50/50 beam splitter (with vacuum padding on the second port),
# creating two channels that interact with the two-level system. We then scan the detuning to extract the excited-state
# population at the evaluation time used in the paper and compare with a coherent-state drive of the same mean photon number.

using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using LinearAlgebra
using PyPlot

#

## symbolic Hilbert space
hu1 = FockSpace(:u1)
hu2 = FockSpace(:u2)
hs = NLevelSpace(:atom, 2)
h = hu1 ⊗ hu2 ⊗ hs

## symbolic operators
au1 = Destroy(h, :au_1, 1)
au2 = Destroy(h, :au_2, 2)
σ(i,j) = Transition(h, :σ, i, j, 3)

## symbolic parameters
@rnumbers γ Δ
gu1 = cnumber("gu_1")
nothing # hide

#

# Two-channel coupling to the atom with a beam splitter.
# The second input is a padding (vacuum) channel.

G_u1 = SLH(1, gu1*au1, 0)
G_u2 = SLH(1, 0, 0)
G_u = concatenate(G_u1, G_u2)

S_bs = [1 1; 1 -1] ./ √2
G_bs = SLH(S_bs, [0, 0], 0)

S_atom = Matrix(I, 2, 2)
L_atom = [√(γ/2)*σ(1,2), √(γ/2)*σ(1,2)]
G_atom = SLH(S_atom, L_atom, Δ*σ(2,2))

G_cas = cascade(G_u, G_bs, G_atom)
H = get_hamiltonian(G_cas)
L = get_lindblad(G_cas)

#

## Pulse parameters 
n = 9
γ_ = 1.0
τ = 0.5/γ_
tw = π^(3/2) / (8*γ_*n)
tp = 6*tw

u(t_) = 1/(√(tw)*π^(1/4)) * exp( -(t_ - tp)^2 / (2*tw^2) ) # TODO
u1(t_) = u(t_)

T0 = 0.0
Tend = tp + τ + 6*tw
dt = tw/30
T = [T0:dt:Tend;]
t1 = tp + τ + 2tw

sum(u.(T))*(T[2]-T[1]) # TODO

gu1_ = u_to_gu_Gauss(tp, tw)
dict_p_t = Dict([gu1, conj(gu1)] .=> [gu1_, t -> conj(gu1_(t))])

#

## numerical bases
bu1 = FockBasis(n)
bu2 = FockBasis(n)
bs = NLevelBasis(2)
b = bu1 ⊗ bu2 ⊗ bs

σ22_QO = translate(σ(2,2), b)

dict_p_Δ(Δn) = Dict([γ, Δ] .=> [γ_, Δn])
H_QO_Δ(Δn) = translate(H, b; parameter=dict_p_Δ(Δn), time_parameter=dict_p_t)
L_QO_Δ(Δn) = [translate(L[i], b; parameter=dict_p_Δ(Δn), time_parameter=dict_p_t) for i=1:length(L)]
nothing # hide
#

ψ0_fock = fockstate(bu1, n) ⊗ fockstate(bu2, 0) ⊗ nlevelstate(bs, 1)

α = sqrt(n/2)
ψ0_coh = coherentstate(bu1, α) ⊗ fockstate(bu2, 0) ⊗ nlevelstate(bs, 1)

#

# Detuning scan and excited-state population at t1.

Δ_ls = range(-6*γ_, 6*γ_, length=41)
pop_fock = zeros(length(Δ_ls))
pop_coh = zeros(length(Δ_ls))
idx_t1 = findmin(abs.(T .- t1))[2]

for (it, Δn) in enumerate(Δ_ls)
    H_QO = H_QO_Δ(Δn)
    L_QO = L_QO_Δ(Δn)
    input_output(t, ρ) = (H_QO(t), [L_QO[i](t) for i in 1:length(L_QO)], [dagger(L_QO[i](t)) for i in 1:length(L_QO)])

    t_, ρt = timeevolution.master_dynamic(T, ψ0_fock, input_output)
    pop_fock[it] = real(expect(σ22_QO, ρt[idx_t1]))

    # t_, ρt = timeevolution.master_dynamic(T, ψ0_coh, input_output)
    # pop_coh[it] = real(expect(σ22_QO, ρt[idx_t1]))
end
nothing # hide

#

close("ramsey-population") # hide
figure("ramsey-population")
plot(Δ_ls, pop_fock, label="Fock, n=9")
# plot(Δ_ls, pop_coh, ls="--", label="coherent, ⟨n⟩=9")
xlabel("detuning Δ/γ")
ylabel(L"\langle \sigma^{22} \rangle (t_1)")
grid(true)
legend()
tight_layout()

#

# ## Package versions

using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    ["QuantumInputOutput", "SecondQuantizedAlgebra", "QuantumOptics", "PyPlot"],
    mode = PKGMODE_MANIFEST,
)
