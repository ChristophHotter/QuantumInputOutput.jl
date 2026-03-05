```@meta
EditURL = "../../../examples/08-1_pulse-delay__simple.jl"
```

# Simple Pulse Delay with a Virtual Cavity

In this example, a single-photon pulse is emitted from an input cavity, delayed by a virtual delay cavity,
and finally captured by an output cavity. The delay cavity is driven by an incoming pulse `u(t)` and
simultaneously emits a delayed pulse `u(t-τ)` using the pulse-shaping couplings introduced in
[V. R. Christiansen and K. Mølmer, Phys. Rev. A 113, 013730 (2026)](https://doi.org/10.1103/PhysRevA.113.013730).

````@example 08-1_pulse-delay__simple
using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using SymbolicUtils
using LinearAlgebra
using PyPlot
````

````@example 08-1_pulse-delay__simple
# symbolic Hilbert space
hu = FockSpace(:u)
hd = FockSpace(:d)
hv = FockSpace(:v)
h = hu ⊗ hd ⊗ hv

# symbolic operators
au = Destroy(h, :a_u, 1)
ad = Destroy(h, :a_d, 2)
av = Destroy(h, :a_v, 3)

# symbolic parameters
gu, gin, gout, gv = cnumbers("g_u g_in g_out g_v")
nothing # hide
````

The input cavity couples ($g_u(t)*a_u$) into the input port of the delay cavity ($g_{in}(t)*a_d$) and the delay cavity couples
the photons via the output port ($g_{out}(t)*a_d$) into the the output cavity ($g_v(t)*a_v$).
This leads to the following cascade of SLH elements.

````@example 08-1_pulse-delay__simple
G_u = SLH(1, gu*au, 0)
G_u2 = concatenate(G_u, SLH(1, 0, 0))

S2 = Matrix(I, 2, 2)
G_d = SLH(S2, [gin*ad, gout*ad], 0)

G_v = SLH(1, gv*av, 0)
G_v2 = concatenate(SLH(1, 0, 0), G_v)

G_cas = cascade(G_u2, G_d, G_v2)
H = get_hamiltonian(G_cas)
L = get_lindblad(G_cas)
nothing # hide
````

````@example 08-1_pulse-delay__simple
# short pulse delay
σ = 1.0
tp = 6*σ
τ = 0.5σ # pulse delay
u(t) = 1/(√(σ)*π^(1/4)) * exp( -(t - tp)^2 / (2*σ^2) )
u_del(t_) = u(t_ - τ)

Tend = 2tp + τ
dt = Tend/5e2
T = [0:dt:Tend;]

gu_ = u_to_gu(u, T)
gout_ = uv_to_gout(u_del, u, T)
gin_ = uv_to_gin(u_del, u, T)
gv_ = v_to_gv(u_del, T)

dict_p_t = Dict([gu, gout, gin, gv] .=> [gu_, gout_, gin_, gv_])
nothing # hide
````

````@example 08-1_pulse-delay__simple
# numeric bases
n = 3
bu = FockBasis(n)
bd = FockBasis(n)
bv = FockBasis(n)
b = bu ⊗ bd ⊗ bv

H_QO = translate(H, b; time_parameter=dict_p_t)
L_QO = [translate(L[i], b; time_parameter=dict_p_t) for i=1:length(L)]

function input_output(t, ρ)
    Ht = H_QO(t)
    Jt = [L_QO[i](t) for i=1:length(L_QO)]
    return Ht, Jt, dagger.(Jt)
end
nothing # hide
````

````@example 08-1_pulse-delay__simple
# time evolution
ψ0 = fockstate(bu, n) ⊗ fockstate(bd, 0) ⊗ fockstate(bv, 0)
t_, ρt = timeevolution.master_dynamic(T, ψ0, input_output)

au_qo = translate(au, b)
ad_qo = translate(ad, b)
av_qo = translate(av, b)

nu = real.(expect(dagger(au_qo)*au_qo, ρt))
nd = real.(expect(dagger(ad_qo)*ad_qo, ρt))
nv = real.(expect(dagger(av_qo)*av_qo, ρt))

nothing # hide
````

````@example 08-1_pulse-delay__simple
close("delay short pulse") # hide
figure("delay short pulse", figsize=(5,3))
plot(T, nu, label=L"\langle a_u^\dagger a_u \rangle")
plot(T, nd, label=L"\langle a_d^\dagger a_d \rangle")
plot(T, nv, label=L"\langle a_v^\dagger a_v \rangle")
xlabel("time")
ylabel("mean photon number")
grid(true)
legend()
tight_layout()
gcf()
````

## Interaction picture for the input and delay cavities

Introducing a separate cavity to delay the pulse is a big disadvantage if one has multiple modes.
To eliminate the delay cavity, we can transform into the interaction picture of the output and delay cavity coupling.
This is, however, only possible if the delay is larger than the pulse, because only then we have $g_{in}(t) \approx g_{v=u}(t)$. Nevertheless, since in most cases only the relative
delay between different modes is crucial, e.g. for two arms of an interferometer, we can simply add a constant delay `Tc >> σ` to all modes.

````@example 08-1_pulse-delay__simple
G_d_in = SLH(S2, [gin*ad, 0], 0)
H_ud = get_hamiltonian(cascade(G_u2, G_d_in))
H_int_sym_ = simplify(H - H_ud)

M(i,j) = cnumber("M_{$(i)$(j)}") # TODO: real? analytic expression
a0_ls = [au, ad]
la = length(a0_ls)
a_int_ls = [sum(M(i,j)*a0_ls[j] for j=1:la) for i=1:la]

int_dict = Dict([a0_ls; adjoint.(a0_ls)] .=> [a_int_ls; adjoint.(a_int_ls)])
nothing # hide
````

````@example 08-1_pulse-delay__simple
H_int_sym = simplify(substitute_operators(H_int_sym_, int_dict))
````

````@example 08-1_pulse-delay__simple
L_int_sym = simplify.(substitute_operators.(L, Ref(int_dict)))
````

````@example 08-1_pulse-delay__simple
# long pulse delay
σ = 1.0
tp = 6*σ
τ = 6σ # pulse delay
u(t) = 1/(√(σ)*π^(1/4)) * exp( -(t - tp)^2 / (2*σ^2) )
u_del(t_) = u(t_ - τ)

Tend = 2tp + τ
dt = Tend/5e2
T = [0:dt:Tend;]

gu_ = u_to_gu(u, T)
gout_ = uv_to_gout(u_del, u, T)
gin_ = uv_to_gin(u_del, u, T)
gv_ = v_to_gv(u_del, T)
nothing # hide
````

````@example 08-1_pulse-delay__simple
# interaction-picture coefficient matrix M(t) for u ↔ d
A_ud = interaction_picture_A_2modes(gu_, gin_)
M_t = interaction_picture_M(A_ud, T) # TODO: analytic M

M_ls = [M(i,j) for i=1:la for j=1:la]
M_t_ls = [t -> M_t(t)[i,j] for i=1:la for j=1:la]
nothing # hide

p_t_sym = [gu, gin, gout, gv, M_ls...]
p_t_num = [gu_, gin_, gout_, gv_, M_t_ls...]
dict_p_t_int = Dict(p_t_sym .=> p_t_num)

au_int = destroy(bu) ⊗ one(bv)
ad_int = one(bu ⊗ bv) # this is amazing
av_int = one(bu) ⊗ destroy(bv)
operators = Dict( [au, au', ad, ad', av, av'] .=>  [au_int, au_int', ad_int, ad_int', av_int, av_int'])

H_int_QO = translate(H_int_sym, b; time_parameter=dict_p_t_int, operators)
L_int_QO = [translate(L_int_sym[i], b; time_parameter=dict_p_t_int, operators) for i=1:length(L_int_sym)]
# H_int_QO = translate(H_int_sym, b; time_parameter=dict_p_t_int)
# L_int_QO = [translate(L_int_sym[i], b; time_parameter=dict_p_t_int) for i=1:length(L_int_sym)]
nothing # hide

function input_output_int(t, ρ)
    Ht = H_int_QO(t)
    Jt = [L_int_QO[i](t) for i=1:length(L_int_QO)]
    return Ht, Jt, dagger.(Jt)
end

ψ0_int = fockstate(bu,3) ⊗ fockstate(bv,0)
t_int, ρt_int = timeevolution.master_dynamic(T, ψ0_int, input_output_int)

nu_int = real.(expect(dagger(au_int)*au_int, ρt_int))
# nd_int = real.(expect(dagger(ad_int)*ad_int, ρt_int)) # hide
nv_int = real.(expect(dagger(av_int)*av_int, ρt_int))
# @show maximum(nd_int) # hide
nothing # hide
````

````@example 08-1_pulse-delay__simple
close("delay int. picture") # hide
figure("delay int. picture", figsize=(5,3))
plot(T, nu_int, label=L"\langle a_u^\dagger a_u \rangle_{IP}")
# plot(T, nd_int, label=L"\langle a_d^\dagger a_d \rangle_{IP}") # hide
plot(T, nv_int, label=L"\langle a_v^\dagger a_v \rangle_{IP}")
xlabel("time")
ylabel("mean photon number")
grid(true)
legend()
tight_layout()
gcf()
````

## Package versions

````@example 08-1_pulse-delay__simple
using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    ["QuantumInputOutput", "SecondQuantizedAlgebra", "QuantumOptics", "PyPlot"],
    mode = PKGMODE_MANIFEST,
)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

