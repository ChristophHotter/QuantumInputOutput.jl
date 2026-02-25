# # Deterministic cat-state generation

# [H. Ritsch, et al., Phys. Rev. A 70, 033804 (2004)](https://doi.org/10.1103%2FPhysRevA.70.033804)
# As usual, we start by loading the packages and specifying the system. 

using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using QuantumCumulants
using PyPlot

#

@rnumbers κ g γ Δa Δc
@cnumbers η gv
Natoms = 2

hc = FockSpace(:cavity)
ha = NLevelSpace(:atom,2)
hv = FockSpace(:output)
h = hc ⊗ ha ⊗ hv

a = Destroy(h,:a,1) # cavity 
σ(i,j) = Transition(h,"σ",i,j,2) # two-level atom
av = Destroy(h,:a_v,3) # output cavity
nothing # hide 

# We couple a classical drive into the cavity through the left mirror $(\kappa_L)$. The decay through the right mirror can be added in several ways: with concatenation, including it already in the initial cavity SLH triple or by simply including the decay term to the Lindblad by hand. In this example, we use the first option. 

H_ac = -Δc*a'a - Δa*σ(2,2) - g*( σ(2,1)*a + σ(1,2)*a' ) + η*σ(2,1) + η'*σ(1,2)
J = [√(κ)*a, √(γ)*σ(1,2)]

# needed later # TODO: move down
G_ac = SLH(1,[√(κ)*a], H_ac) 
G_v = SLH(1,[gv*av], 0)

## numerical parameters
κ_ = γ_ = 1.0
g_ = 100.0
Δa_ = 0.0
Δc_ = 0.0

p_sym = [κ , g , γ , Δa , Δc ]
p_num = [κ_, g_, γ_, Δa_, Δc_]
dict_p1 = Dict(p_sym .=> p_num)

# Gaussian pulse (definiton?)
# η0 = 9g_
# τ = 5/g_
# η_t(t) = η0*exp( -(t-2τ)^2 / (0.5*τ^2) )
# Tend = 8τ

# rectengular pulse
η0 = √(8)*g_
T_pulse = 2π/g_
η_t(t) = (t<T_pulse)*1.0
dict_t = Dict( [η, conj(η)] .=> [η_t, ηc_t] )

Tend = 4/κ_
T = [0:0.001:1;]*Tend

# pygui(true)
# close("pulse")
# figure("pulse")
# plot(T*g_, η_t.(T)/η0)
# grid(true)
# xlabel("gt")


## numerical basis
bc = FockBasis(10)
ba = NLevelBasis(2)
bv = FockBasis(10)
b = bc ⊗ ba ⊗ bv
#
a_QO = destroy(bc) ⊗ one(ba) 
ad_QO = dagger(a_QO)
σ_QO(i,j) = one(bc) ⊗ transition(ba,i,j)

ops_sym = [a, a', σ(2,2), σ(1,2), σ(2,1)]
ops_QO = [a_QO, dagger(a_QO), σ_QO(2,2), σ_QO(1,2), σ_QO(2,1)]
ops_dict = Dict(ops_sym .=> ops_QO)

H_QO = translate(H_ac, b; parameter=dict_p1, time_dep_param=dict_t, operators=ops_dict)
L_QO = translate(√(κ)*a, b; parameter=dict_p1, time_dep_param=dict_t, operators=ops_dict)
J_add_QO = []

#

## time evolution
T = [0:0.01:1;]*20
ψ0 = fockstate(bc1,0) 
t_, ρt = timeevolution.master(T, ψ0, H1_QO, J1_QO)
nothing # hide

#

n_cavity = real(expect(a_QO'a_QO, ρt))
n_ref = real(expect(dagger(J1_QO[1])*J1_QO[1], ρt))
n_trans = real(expect(dagger(J1_QO[2])*J1_QO[2], ρt))
nothing # hide

#

close("time evolution") # hide
figure("time evolution")
subplot(2,1,1)
title("n(t)")
plot(t_, n_cavity)
xlabel("t")
ylabel("cavity photons")
grid(true)

subplot(2,1,2)
title("transmission - reflection")
plot(t_, n_ref; label="reflection")
plot(t_, n_trans; label="transmission", ls="--")
xlabel("t")
ylabel("intensity rate")
grid(true)
legend()
tight_layout()
gcf()

# Now we scan the laser-cavity detuning $\Delta$ to plot the transmission and reflection spectrum. 

dict_p_Δ(Δn) = Dict(p_sym .=> [En, κ_Rn, κ_Ln, Δn])
H1_QO_Δ_(Δn) = translate(H1, bc1; parameter=dict_p_Δ(Δn), operators=ops_dict)

n_ref_Δ = zeros(lΔ)
n_trans_Δ = zeros(lΔ)

for it=1:lΔ
    Δn_ = Δn_ls[it]
    t_it, ρt_it = timeevolution.master(T, ψ0, H1_QO_Δ_(Δn_), J1_QO)

    n_ref_Δ[it] = real(expect(dagger(J1_QO[1])*J1_QO[1], ρt_it[end]))
    n_trans_Δ[it] = real(expect(dagger(J1_QO[2])*J1_QO[2], ρt_it[end]))
end
nothing # hide

#

close("spectrum") # hide
figure("spectrum")
plot(Δn_ls, n_ref_Δ; label="reflection")
plot(Δn_ls, n_trans_Δ; label="transmission", ls="--")
xlabel("Δ")
grid(true)
legend()
gcf()

# ## Two-sided cavity with atoms

# In the following, we include $N=2$ two-level atoms in the cavity and simulate the transmission and reflection of coherent Gaussian pulse with with a mean photon number of $|\alpha|^2 = 1/10$. We assume that the atoms are on resonance with the cavity, i.e. $\Delta = \Delta_c = \Delta_a$.

H_ac = -Δ*(a'a + ∑σ(2,2)) + g*(a'∑σ(1,2) + a*∑σ(2,1))
G_ac = SLH(1, √κ_L*a, H_ac)
G_ac_drive = (G_d ▷ G_ac) ⊞ SLH(1, √κ_R*a, 0)
nothing # hide

#

H2 = G_ac_drive.hamiltonian

#

L2_L = G_ac_drive.lindblad[1]

#

L2_R = G_ac_drive.lindblad[2]

#

## numerical parameter
κ_Rn2 = κ_Ln2 = 2π*1
gn2 = 2π*0.4
Δn2 = 0.0
γn = 2π*0.1

p_sym2 = [κ_R, κ_L, Δ, g]
p_num2 = [κ_Rn2, κ_Ln2, Δn2, gn2]

σp = 10/κ_Ln2 # pulse width 
Tp = 4σp # pulse peak 
Tend = 3Tp
α0 = √(0.1) # √ of total photon number
Ω0 = α0*2*√(κ_Ln2)/(π^(1/4)*√(σp)) 
Ω1(t) = Ω0/2*exp( -(t-Tp)^2 / (2*σp^2) )
E_t(t) = Ω1(t)/√(κ_Ln2)

T = [0:0.001:1;]*Tend
ΔT = T[2] - T[1]
n_pulse = round(sum(abs2.(E_t.(T)))*ΔT, digits=7)
@show n_pulse 

dict_p2 = Dict(p_sym2 .=> p_num2)
dict_p_t2 = Dict( [E, conj(E)] .=> [E_t, E_t] )
nothing # hide

#

## numeric bases
bc1 = FockBasis(4)
ba = NLevelBasis(2)
b = bc1 ⊗ tensor([ba for i=1:Natoms]...)

a_QO2 = to_numeric(a,b)
σ_QO(α,i,j) = to_numeric(σ(α,i,j),b)

## translate to numeric Hamiltonian and Lindblad
H_QO = translate(H2, b; parameter=dict_p2, time_dep_param=dict_p_t2)
L2_L_QO = translate(L2_L, b; parameter=dict_p2, time_dep_param=dict_p_t2)
L2_R_QO = translate(L2_R, b; parameter=dict_p2)

## additional atomic decay into free space
J_add = [√(γn)*σ_QO(α,1,2) for α=1:Natoms]

function input_output(t,ρ)
    H = H_QO(t)
    J = [L2_L_QO(t), L2_R_QO, J_add...]
    return H, J, dagger.(J)
end
nothing # hide

#

## time evolution
ψ0 = fockstate(bc1,0) ⊗ tensor([nlevelstate(ba,1) for i=1:Natoms ]...)
t2_, ρt2 = timeevolution.master_dynamic(T, ψ0, input_output)
nothing # hide

#

L2_L_QO_dag(t) = dagger(L2_L_QO(t))
l_t = length(t2_)
n_trans2 = zeros(l_t)
n_ref2 = zeros(l_t)
for it=1:l_t
    n_trans2[it] = abs(expect(dagger(L2_R_QO)*L2_R_QO, ρt2[it]))
    n_ref2[it] = abs(expect( L2_L_QO_dag(t2_[it])*L2_L_QO(t2_[it]), ρt2[it]))
end
nothing # hide

close("time evolution") # hide
figure("time evolution")
plot(t2_, n_trans2, label="transmission = $(round(sum(n_trans2)*ΔT/n_pulse*100))%")
plot(t2_, n_ref2, ls="--", label="reflection = $(round(sum(n_ref2)*ΔT/n_pulse*100))%")
xlabel("t")
legend()
grid(true)
tight_layout()
gcf()

# We can see that only about 5% is transmitted and 71% are reflected. The rest is scattered into free space by the atoms.