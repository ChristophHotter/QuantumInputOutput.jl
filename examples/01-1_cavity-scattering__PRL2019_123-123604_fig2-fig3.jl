# # Cavity Scattering of a Single Photon

# In this example, we simulate the scattering of a resonant single photon on an empty one-sided cavity. The temporal mode of the light pulse is a Gaussian with width $\sigma$ and the cavity has a decay rate $\gamma$. This system has been studied in [A. Kiilerich, et. al., Phys. Rev. Lett. 123, 123604 (2019)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.123604). We start by loading the needed packages and specifying the model. 

using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using PyPlot
using LinearAlgebra

#

## symbolic Hilbert space
hu1 = FockSpace(:u1) 
hc1 = FockSpace(:c1)
hv1 = FockSpace(:v1)
h = hu1 ⊗ hc1 ⊗ hv1

## symbolic operators
au = Destroy(h,:a_u,1) 
c = Destroy(h,:c,2)
av = Destroy(h,:a_v,3)

## symbolic parameters
gu, Δ, γ = rnumbers("g_u Δ γ") 
gv = cnumber("g_v")
nothing # hide

# We use the symbolic operators and paramters to define the SLH triples and cascade them to obtain the Hamiltonian and Lindblad for the system. 

G_u = SLH(1, gu*au, 0) # input cavity 
G_c = SLH(1, √(γ)*c, Δ*c'c) # system cavity
G_v = SLH(1, gv*av, 0) # output cavity

G_cas = ▷(G_u, G_c, G_v)
nothing # hide

#

H = get_hamiltonian(G_cas)

#

L = get_lindblad(G_cas)[1] # only one Lindblad term in this example

# To solve the dynamics of the system we translate the symbolic expressions into numeric operators (matrices) of QuantumOpitcs.jl. To do so, we define the numerical parameters and operator basis. 

## numerical parameters 
γ_ = 1.0
Δ_ = 0.0

p_sym = [γ , Δ , gv]
p_num = [γ_, Δ_, 0 ] # gv=0
dict_p = Dict(p_sym .=> p_num);

## Gaussian input mode
T_p = 1/γ_
T_end = 12T_p
σ = sqrt(0.5)*T_p
u1(t) = sqrt(1/(σ*√(2π))*exp( -0.5*(t - 4T_p)^2/σ^2 ))
T = [0:0.002:1;]*T_end
ΔT = T[2] - T[1]

gu_int = u_to_gu(u1, T) # interpolation
gu_t(t) = gu_int(t)
dict_p_t = Dict(gu => gu_t)

## numeric bases 
bu1 = FockBasis(1) 
bc1 = FockBasis(1)
bv1 = FockBasis(1)
b = bu1 ⊗ bc1 ⊗ bv1
nothing # hide

# We now use the function `translate()` to create the numeric operators. If the kwarg `time_dep_param` is provided the created operator is a time-dependent function. 

H_QO = translate(H, b; parameter=dict_p, time_dep_param=dict_p_t)
L_QO = translate(L, b; parameter=dict_p, time_dep_param=dict_p_t)
nothing # hide


# To solve the dynamics we use the QuantumOptics.jl function `timeevolution.master_dynamic()`. 

## time-depedent function for timeevolution.master_dynamic() that returns H(t), J(t) and Jd(t)
function input_output_1(t,ρ)
    H = H_QO(t)
    J = [L_QO(t)]
    return H, J, dagger.(J)
end;

## initial state
ψ0 = fockstate(bu1,1) ⊗ fockstate(bc1,0) ⊗ fockstate(bv1,0)

## time evolution
t_, ρt = timeevolution.master_dynamic(T, ψ0, input_output_1) 
nothing # hide

# To calculate expectation values we create the desired numerical operators.

au_qo = translate(au, b)
c_qo = translate(c, b)
av_qo = translate(av, b)

n_c_t = real.(expect(c_qo'*c_qo, ρt))
n_u1_t = real.(expect(au_qo'*au_qo, ρt))
nothing # hide

# In order to determine suitable temporal output modes we calculate the two-time autocorrelation function $g^{(1)}(t_1,t_2) = \langle L_s^\dagger(t_1) L_s(t_2) \rangle$ and diagonalize the matrix to obtain the eigenvalues with the corresponding eigenvectors. The eigenvalues correspond to the mean photon number $n_i$ in the corresponding temporal eigenvector mode $v_i$.  

Ls(t) = gu_t(t)*au_qo + √(γ_)*c_qo
g1_m = two_time_corr_matrix(T, ρt, input_output_1, Ls)
nothing # hide

#

close("all") # hide
figure("g1(t1,t2) matrix", figsize=(4,3.5))
pcolormesh(T, T, real.(g1_m), cmap="inferno")
xlabel(L"\gamma t_2")
ylabel(L"\gamma t_1")
tight_layout()
colorbar(label=L"g^{(1)}(t_1,t_2)")
gcf() 

# The eigenvalues and corresponding eigenvectors are sorted in ascending order, which means the last eigenvalue corresponds to the highest populated temporal mode. 

F = eigen(g1_m)
n_avg = round.(real.(F.values)*ΔT; digits=3)
modes = F.vectors
v_mode = (modes[:,end]) / sqrt(ΔT)

@show n_avg[end-1:end]
nothing # hide

# We use now the mode with the highest mean photon number as our out-mode to determine its quantum state. Due to a problem with the conjugate of a function we also need to provide the $g_v^*(t)$ in the dictionary for the time-dependent paramters. 

p_sym_2 = [γ , Δ ]
p_num_2 = [γ_, Δ_]
dict_p_2 = Dict(p_sym_2 .=> p_num_2);

## time-depedent coupling for the output mode $v(t)$
gv_t_ = v_to_gv(v_mode, T)
gv_t(t) = gv_t_(t)
gvc_t(t) = conj(gv_t_(t))

dict_p_t_2 = Dict([gu, gv, conj(gv)] .=> [gu_t, gv_t, gvc_t]);

H_QO_2 = translate(H, b; parameter=dict_p_2, time_dep_param=dict_p_t_2)
L_QO_2 = translate(L, b; parameter=dict_p_2, time_dep_param=dict_p_t_2)
function input_output_2(t,ρ)
    H = H_QO_2(t)
    J = [L_QO_2(t)]
    return H, J, dagger.(J)
end;

## time evolution for the system including the output cavity
t_2, ρt_2 = timeevolution.master_dynamic(T, ψ0, input_output_2)

n_v1_t = real.(expect(av_qo'*av_qo, ρt_2))
nothing # hide

#

figure("modes")
subplot(2,1,1)
plot(T, u1.(T), ls="--", label="u", color="red")
fill_between(T, 0., u1.(T), alpha=0.5, color="red")

plot(T, real.(v_mode), color="blue", label="v")
fill_between(T, 0., real.(v_mode), alpha=0.5, color="blue")
xlim(0,12)
ylim(-0.8,0.8)
yticks([-0.8,0,0.8])
ylabel("modes")
legend()

twinx()
xlim(0,12)
ylim(0,8)
plot(T, abs2.(gu_t.(T)), color="red")
plot(T[3:end], abs2.(gv_t.(T))[3:end], color="blue", ls="--")
ylabel("Rates (γ)")

subplot(2,1,2)
plot(T, n_u1_t, label=L"\langle a^\dagger a \rangle_u", color="red")
plot(T, n_c_t, label=L"\langle c^\dagger c \rangle", color="green")
plot(T, n_v1_t, label=L"\langle a^\dagger a \rangle_v", color="blue")
xlim(0,12)
ylim(0,1)
xlabel("time (1/γ)")
ylabel("Exciations")
legend()
gcf() 

# ## Cavity with phase noise

# We slightly adapt the above example by assuming the intial pulse to be in a coherent state and adding phase noise to the cavity. This results in scattering into multiple modes. 

## new basis of the system
bu1_3 = FockBasis(12)
bc1_3 = FockBasis(6)
bv1_3 = FockBasis(6)
b_3 = tensor(bu1_3,bc1_3,bv1_3)

## new operators of the system
au_3 = embed(b_3, 1, destroy(bu1_3))
c_3 = embed(b_3, 2, destroy(bc1_3))
av_3 = embed(b_3, 3, destroy(bv1_3))
cdc_3 = c_3'c_3

## we use the same Hamiltonian as before but add a depasing term to the dissipation 
H_QO_3 = translate(H, b_3; parameter=dict_p_2, time_dep_param=dict_p_t_2)
L_QO_3 = translate(L, b_3; parameter=dict_p_2, time_dep_param=dict_p_t_2)
function input_output_3(t,ρ)
    H = H_QO_3(t)
    J = [L_QO_3(t), √(γ_)*cdc_3]
    return H, J, dagger.(J)
end;

#

# Due to the larger Hilbert space the time evolution takes a few seconds. 

ψ0_3 = coherentstate(bu1_3,2) ⊗ fockstate(bc1_3,0) ⊗ fockstate(bv1_3,0)
t_3, ρt_3 = timeevolution.master_dynamic(T, ψ0_3, input_output_3)
nothing # hide

#

L0(t) = √(γ_)*c_3 + gu_t(t)*au_3 + gv_t(t)*av_3
I_out = [expect(dagger(L0(t_3[i]))*L0(t_3[i]), ρt_3[i]) for i=1:length(t_3)]

n_u1_t_3 = real.(expect(au'*au, ρt_3))
n_v1_t_3 = real.(expect(av'*av, ρt_3))
nothing # hide

#

figure("dephasing", figsize=(6,3))
plot(t_3, n_u1_t_3, label=L"\langle a^\dagger a \rangle_u", color="red")
plot(t_3, n_v1_t_3, label=L"\langle a^\dagger a \rangle_v", color="blue", ls="--")
plot(t_3, I_out, label=L"I_{out}", color="black", ls="dotted")
xlim(0,12)
ylim(0,4)
xlabel("time (1/γ)")
ylabel("expectation values")
legend()
gcf()

# ## Package versions

# These results were obtained using the following versions:

using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    ["QuantumInputOutput", "SecondQuantizedAlgebra", "QuantumOpitcs", "PyPlot"],
    mode = PKGMODE_MANIFEST,
)
