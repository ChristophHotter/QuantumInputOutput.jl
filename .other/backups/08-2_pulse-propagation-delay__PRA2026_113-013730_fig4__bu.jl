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
hu = FockSpace(:u)
hd1 = FockSpace(:d1)
hd2 = FockSpace(:d2)
hs = NLevelSpace(:atom, 2)
h = hu ⊗ hd1 ⊗ hd2 ⊗ hs

## symbolic operators
au = Destroy(h, :au, 1)
ad1 = Destroy(h, :ad_1, 2)
ad2 = Destroy(h, :ad_2, 3)
σ(i, j) = Transition(h, :σ, i, j, 4)

## symbolic parameters
@rnumbers γ Δ r t
gu, gd1in, gd1out, gd2in, gd2out = cnumbers("gu gd1_{in} gd1_{out} gd2_{in} gd2_{out}")
nothing # hide

#

# Two-channel coupling to the atom with a beam splitter.
# The second input is a padding (vacuum) channel.

G_u = concatenate(SLH(1, gu*au, 0), SLH(1, 0, 0))

S_bs = [r t; t -r]
G_bs = SLH(S_bs, [0, 0], 0)

I2 = Matrix(I, 2, 2)
G_u_bs = concatenate(G_u ▷ G_bs, SLH(I2, [0, 0], 0))

I4 = Matrix(I, 4, 4)
G_d1_d2 = SLH(I4, [gd1in*ad1, gd2in*ad2, gd1out*ad1, gd2out*ad2], 0) # input at port 1 & 2, output at port 3 & 4

G_u_bs_d1_d2 = G_u_bs ▷ G_d1_d2
nothing #hide 

#

G_u_bs_d1_d2.hamiltonian

L_atom = [0, 0, √(γ/2)*σ(1, 2), √(γ/2)*σ(1, 2)]
G_atom = SLH(I4, L_atom, Δ*σ(2, 2))

G_u_bs_d1_d2_atom = G_u_bs_d1_d2 ▷ G_atom
H = get_hamiltonian(G_u_bs_d1_d2_atom)
L = get_lindblad(G_u_bs_d1_d2_atom)

#

## Pulse parameters 
n = 5#9
γ_ = 1.0
τ = 0.5/γ_
tw = π^(3/2) / (8*γ_*n)
tp = 4*tw

u(t) = 1/(√(tw)*π^(1/4)) * exp(-(t - tp)^2 / (2*tw^2))
ud1(t) = u(t-tp)
ud2(t) = u(t-tp-τ)

Tend = 2tp + τ + 2tw
dt = tw/30
T = [0:dt:Tend;]
t1 = 2tp + τ + 2tw

rn = 1/√(2)
tn = 1/√(2)

gu_ = u_to_gu(u, T)
gd1in_ = v_to_gv(u, T)
gd1out_ = u_to_gu(ud1, T)
gd2in_ = v_to_gv(u, T)
gd2out_ = u_to_gu(ud2, T)
dict_p_t =
    Dict([gu, gd1in, gd1out, gd2in, gd2out] .=> [gu_, gd1in_, gd1out_, gd2in_, gd2out_])

#

## numerical bases
bu = FockBasis(n)
bd1 = FockBasis(n)
bd2 = FockBasis(n)
bs = NLevelBasis(2)
b = bu ⊗ bd1 ⊗ bd2 ⊗ bs

σ22_QO = translate_qo(σ(2, 2), b)

dict_p_Δ(Δn) = Dict([γ, Δ, r, t] .=> [γ_, Δn, rn, tn])
H_QO_Δ(Δn) = translate_qo(H, b; parameter = dict_p_Δ(Δn), time_parameter = dict_p_t)
L_QO_Δ(Δn) = [
    translate_qo(L[i], b; parameter = dict_p_Δ(Δn), time_parameter = dict_p_t) for
    i = 1:length(L)
]
nothing # hide
#

ψ0_fock = fockstate(bu, n) ⊗ fockstate(bd1, 0) ⊗ fockstate(bd2, 0) ⊗ nlevelstate(bs, 1)

#

# Detuning scan and excited-state population at t1.

Δ_ls = range(-20γ_, 20γ_, length = 80)
pop_fock = zeros(length(Δ_ls))
pop_coh = zeros(length(Δ_ls))
idx_t1 = findmin(abs.(T .- t1))[2]
using ProgressMeter
prog = Progress(length(Δ_ls))

for (it, Δn) in enumerate(Δ_ls)
    Δn = Δ_ls[it]

    H_QO = H_QO_Δ(Δn)
    L_QO = L_QO_Δ(Δn)

    input_output =
        (t, ρ) -> (
            H_QO(t),
            [L_QO[i](t) for i = 1:length(L_QO)],
            [dagger(L_QO[i](t)) for i = 1:length(L_QO)],
        )

    t_, ρt = timeevolution.master_dynamic([T[1], T[end]], ψ0_fock, input_output)
    # pop_fock[it] = real(expect(σ22_QO, ρt[idx_t1]))
    pop_fock[it] = real(expect(σ22_QO, ρt[end]))

    # t_, ρt = timeevolution.master_dynamic(T, ψ0_coh, input_output)
    # pop_coh[it] = real(expect(σ22_QO, ρt[idx_t1]))
    next!(prog)
end
nothing # hide

#

pygui(true)
close("ramsey-population") # hide
figure("ramsey-population")
plot(Δ_ls, pop_fock, label = "Fock, n=$(n)")
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
