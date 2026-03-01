# # Interaction Picture: Rabi Oscillations with a Quantum Pulse (Fig. 2, PRA 107, 013706)
#
# Reproduces the interaction-picture example from
# Christiansen et al., Phys. Rev. A 107, 013706 (2023), Fig. 2.

using QuantumInputOutput
using QuantumOptics
using SecondQuantizedAlgebra
using SymbolicUtils
using PyPlot

# Parameters
γ_ = 1.0
n_photons = 20

# Gaussian pulse u(t), Eq. (21)
τ = 1 / γ_
t_p = 4 / γ_
u(t) = 1 / (sqrt(τ) * π^(1/4)) * exp(-0.5 * ((t - t_p) / τ)^2)

T_end = 12.0
T = [0:0.002:1;] * T_end

# Virtual-cavity couplings for input/output modes
gu_t = u_to_gu(u, T)
gv_t = v_to_gv(u, T) # identical output mode v(t) = u(t)

# Interaction-picture coefficient matrix M(t) for (u, v)
A_uv = interaction_picture_A_2modes(gu_t, gv_t)
M_t = interaction_picture_M(A_uv, T)

# Symbolic SLH setup for u -> s -> v and u -> v
hu = FockSpace(:u)
hs = NLevelSpace(:s, 2)
hv = FockSpace(:v)
h = hu ⊗ hs ⊗ hv

au_sym = Destroy(h, :a_u, 1)
av_sym = Destroy(h, :a_v, 3)
σ12_sym = Transition(h, :σ, 1, 2, 2)

gu, γ, gv = rnumbers("gu γ gv")

G_u = SLH(1, gu * au_sym, 0)
G_s = SLH(1, sqrt(γ) * σ12_sym, 0)
G_v = SLH(1, gv * av_sym, 0)

H = get_hamiltonian(▷(G_u, G_s, G_v))
H0 = get_hamiltonian(▷(G_u, G_v))
H_int_sym_ = simplify(H - H0)
L_sym_ = get_lindblad(▷(G_u, G_s, G_v))[1]

M(i,j) = Int(i!=j)*cnumber("M$(i)$(j)") # M_ii = 0
a0_ls = [au_sym, av_sym]
la = length(a0_ls)
a_int_ls = [sum(M(i,j)*a0_ls[j] for j=1:la) for i=1:la]

int_dict = Dict(a0_ls .=> a_int_ls)
H_int_sym = simplify(substitute(H_int_sym_, int_dict))
L_sym = simplify(substitute(L_sym_, int_dict))

# Bases: u-mode (truncated), two-level atom, v-mode (small truncation)
bu = FockBasis(n_photons)#, n_photons-4)
ba = NLevelBasis(2)
bv = FockBasis(4)
b = bu ⊗ ba ⊗ bv

# # Operators
# au = destroy(bu) ⊗ one(ba) ⊗ one(bv)
# av = one(bu) ⊗ one(ba) ⊗ destroy(bv)
# σm = one(bu) ⊗ transition(ba, 1, 2) ⊗ one(bv)
# σe = one(bu) ⊗ transition(ba, 2, 2) ⊗ one(bv)

dict_p = Dict(γ => γ_)
p_t_sym = [gu, gv, M(1,2), M(2,1), conj(M(1,2)), conj(M(2,1))]
p_t_num = [gu_t, gv_t, t -> M_t(t)[1,2], t -> M_t(t)[2,1], t -> conj(M_t(t)[1,2]), t -> conj(M_t(t)[2,1])]
dict_p_t = Dict(p_t_sym .=> p_t_num)

#

H_int_QO = translate(H_int_sym, b; parameter=dict_p, time_parameter=dict_p_t)
L_QO = translate(L_sym, b; parameter=dict_p, time_parameter=dict_p_t)


# continue here! # TODO


function au_I(t)
    M = M_t(t)
    return M[1, 1] * au + M[1, 2] * av
end

function av_I(t)
    M = M_t(t)
    return M[2, 1] * au + M[2, 2] * av
end

# Interaction-picture Hamiltonian and Lindblad operator
function H_I(t)
    ops = Dict(au_sym => au_I(t), av_sym => av_I(t), σ12_sym => σm)
    return translate(H_int_sym, b; parameter=dict_p, time_parameter=dict_t, operators=ops)(t)
end

function L_I(t)
    ops = Dict(au_sym => au_I(t), av_sym => av_I(t), σ12_sym => σm)
    return translate(L_sym, b; parameter=dict_p, time_parameter=dict_t, operators=ops)(t)
end

function input_output(t, ρ)
    Ht = H_I(t)
    J = [L_I(t)]
    return Ht, J, dagger.(J)
end

# Initial state: |n⟩ in u mode, atom in ground state, vacuum in v mode
ψ0 = fockstate(bu, n_photons) ⊗ nlevelstate(ba, 1) ⊗ fockstate(bv, 0)

# Time evolution
_, ρt = timeevolution.master_dynamic(T, ψ0, input_output)

# Observables
n_u = real.(expect(au' * au, ρt))
n_v = real.(expect(av' * av, ρt))
P_e = real.(expect(σe, ρt))

# Plotting
close("interaction picture fig2")
figure("interaction picture fig2", figsize = (5.6, 4.2))

subplot(2, 1, 1)
plot(T, n_u, label = L"\langle n_u \rangle")
plot(T, n_v, label = L"\langle n_v \rangle", ls = "--")
xlabel("time")
ylabel("mean photons")
legend()

subplot(2, 1, 2)
plot(T, P_e, label = L"\langle \sigma_{ee} \rangle")
plot(T, n_u, label = L"\langle n_u \rangle", ls = "--")
xlabel("time")
ylabel("population")
legend()

tight_layout()
