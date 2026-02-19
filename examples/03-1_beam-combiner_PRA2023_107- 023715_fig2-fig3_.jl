# # Perfect Splitting of a Two-photon Pulse

# TODO: Intro

using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using PyPlot
using LinearAlgebra

#

## symbolic Hilbert space
hu2 = FockSpace(:u2)
hu1 = FockSpace(:u1)
hs1 = NLevelSpace(:atom, 2)
hv1 = FockSpace(:v1)
h = hu2 ⊗ hu1 ⊗ hs1 ⊗ hv1

## symbolic operators
au2 = Destroy(h,:au_2,1)
au1 = Destroy(h,:au_1,2)
σ(i,j) = Transition(h,:σ,i,j,3)
av1 = Destroy(h,:av_2,4)

## symbolic parameters
@rnumbers γ Δ
gu1, gu2, gv1 = cnumbers("gu_1 gu_2 gv_1"); 

# We use the symbolic operators and paramters to define the SLH triples and cascade them to obtain the Hamiltonian and Lindblad for the system. 

G_u2 = SLH(1, gu2*au2, 0) # input cavity 2
G_u1 = SLH(1, gu1*au1, 0) # input cavity 1
G_a = SLH(1, √(γ)*σ(1,2), Δ*σ(2,2)) # scattering atom
G_v1 = SLH(1, gv1*av1, 0) # output cavity 1

G_cas = cascade(G_u2, G_u1, G_a, G_v1)
nothing # hide

#

H = get_hamiltonian(G_cas)

#

L = get_lindblad(G_cas)[1] # only one Lindblad in this example

# Next, the numerical parameters and functions of the system are defined.

γ_ = 1.0
Δ_ = 0.0

p_sym = [γ , Δ , gu2, gv1]
p_num = [γ_, Δ_, 0  , 0  ]
dict_p = Dict(p_sym .=> p_num);

## GAussian input pulse
τ = 0.38; tp = 4/γ_
u(t) = 1/(sqrt(τ)*π^(1/4)) * exp( -(t - tp)^2 / (2*τ^2) )
T = [0:0.002:1;]*20
ΔT = T[2] - T[1]

gu_int = u_to_gu(u, T) # interpolation
gu_(t) = gu_int(t)
gu_c(t) = conj.(gu_int(t))

dict_p_t = Dict([gu1, conj(gu1)] .=> [gu_, gu_c]);

# We translate the symbolic expressions to numerical operators and solve the time-dependent master equation with QuantumOptics.jl. 

## numeric bases 
bu2 = FockBasis(2)
bu1 = FockBasis(2)
bs1 = NLevelBasis(2)
bv1 = FockBasis(2)
b = bu2 ⊗ bu1 ⊗ bs1 ⊗ bv1;

H_QO = translate(H, b; parameter=dict_p, time_dep_param=dict_p_t)
L_QO = translate(L, b; parameter=dict_p, time_dep_param=dict_p_t)
function input_output(t,ρ)
    H = H_QO(t)
    J = [L_QO(t)]
    return H, J, dagger.(J)
end
nothing # hide

#

## time evolution
ψ0 = fockstate(bu2,0) ⊗ fockstate(bu1,2) ⊗ nlevelstate(bs1,1) ⊗ fockstate(bv1,0) 
t_, ρt = timeevolution.master_dynamic(T, ψ0, input_output; abstol, reltol)
nothing # hide


###########################################
################ Hier weiter ! ############


# autocorrelation function
au1_qo = to_numeric(au1,b)
σ_qo(i,j) = to_numeric(σ(i,j),b)

Ls(t) = (gu_(t))'*au1_qo + √(γ_)*σ_qo(1,2)
g1_m = two_time_corr_matrix(T, ρt, input_output, Ls);

close("g1(t1,t2) matrix")
figure("g1(t1,t2) matrix", figsize=(4,3.5))
pcolormesh(T, T, real.(g1_m), cmap="inferno")
xlabel("γ t2")
ylabel("γ t1")
colorbar();
tight_layout()
display(gcf())
savefig("figures/05-1_PRA2023_g1.png")

F = eigen(g1_m)
n_avg = real.(F.values)*ΔT
modes = F.vectors
v1_mode = (modes[:,end]) / √(ΔT)
v2_mode = (modes[:,end-1]) / √(ΔT)

v1_p = 1/√(2) * ( v1_mode - v2_mode ) 
v2_p = 1/√(2) * ( v1_mode + v2_mode );

@show sum(abs2.(v1_p))*ΔT;
@show sum(abs2.(v2_p))*ΔT;

close("modes")
figure("modes")
plot(t_, -real.(v1_mode), color="black")
plot(t_, -real.(v2_mode), color="red", ls="--")
xlabel("time (1/γ)")
ylabel("output mode")
display(gcf())
savefig("figures/05-1_PRA2023_output-modes_v1-v2.png")

# output mode
v1_new(t) = (u(T[end]-t))'

# input modes
u1_new = conj.(reverse(v1_p))
u2_new = conj.(reverse(v2_p));

# order of the input functions needs to be in ascending order [1, 2, 3, ...]
# [cascaded from right to left]
u_new_data = [u1_new, u2_new]
u_new_fct = [LinearInterpolation(u, T) for u in u_new_data];

# the coupling of the u1 cavity needs no adaptation 
gu1_int = u_to_gu(u1_new, T)
gu1_(t) = gu1_int(t)
gu1_c(t) = conj(gu1_(t))

# gu2 needs to take into account to go scatter at the u1 cavity
# the last argument in ui_to_u_i_im1 (=2) describes the number of the input cavity
u2_for_gu2 =  ui_to_u_i_im1(u_new_fct, T, 2)
gu2_int = u_to_gu(u2_for_gu2,T)
gu2_(t) = gu2_int(t)
gu2_c(t) = conj(gu2_(t))
#
gv1_int = v_to_gv(v1_new,T)
gv1_(t) = gv1_int(t)
gv1_c(t) = conj(gv1_(t))

g_sym = [gu1, gu2, gv1, conj(gu1), conj(gu2), conj(gv1)]
g_num = [gu1_, gu2_, gv1_, gu1_c, gu2_c, gv1_c]

dict_p_t_out = Dict(g_sym .=> g_num);
;

p_sym_out = [γ , Δ ]
p_num_out = [γ_, Δ_]
dict_p_out = Dict(p_sym_out .=> p_num_out);

H_QO_2 = translate(H, b; parameter=dict_p_out, time_dep_param=dict_p_t_out)
L_QO_2 = translate(L, b; parameter=dict_p_out, time_dep_param=dict_p_t_out)
function input_output_2(t,ρ)
    H = H_QO_2(t)
    J = [L_QO_2(t)]
    return H, J, dagger.(J)
end;

# time evolution
ψ0_out = fockstate(bu2,1) ⊗ fockstate(bu1,1) ⊗ nlevelstate(bs1,1) ⊗ fockstate(bv1,0) 
t_2, ρt_2 = timeevolution.master_dynamic(T, ψ0_out, input_output_2; abstol, reltol);

nu1_t_comb = real(expect(au1'*au1, ρt_2))
nu2_t_comb = real(expect(au2'*au2, ρt_2))
nv1_t_comb = real(expect(av1'*av1, ρt_2))
s22_t_comb = real(expect(σ(2,2), ρt_2));

close("beam combiner")
figure("beam combiner")
plot(T, nu2_t_comb, color="red", ls="--", label="⟨a⁺a⟩ u2")
plot(T, nu1_t_comb, color="blue", ls="-.", label="⟨a⁺a⟩ u1")
plot(T, s22_t_comb, color="black", ls="dotted", label="⟨σ²²⟩")
plot(T, nv1_t_comb, color="green", ls="-", label="⟨a⁺a⟩ v1")
ylim(0,2)
xlim(10,18)
xlabel("time (1/γ)")
ylabel("Mean Excitation")
legend()
tight_layout()
display(gcf())
savefig("figures/05-1__PRA2023_beam-combiner.png")
