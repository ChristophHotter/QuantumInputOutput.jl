```@meta
EditURL = "../../../examples/06-1_interaction-picture__PRA2023_107-013706_fig2.jl"
```

# Interaction Picture Scattering with a Quantum Pulse

In this example, we study the scattering of a Fock state $| n = 20 \rangle$ on a two-level system, using the interaction picture introduced in
[Christiansen et al., Phys. Rev. A 107, 013706 (2023)](https://doi.org/10.1103/PhysRevA.107.013706).
We start by loading the packages and specifying the model.

````@example 06-1_interaction-picture__PRA2023_107-013706_fig2
using QuantumInputOutput
using QuantumOptics
using SecondQuantizedAlgebra
using SymbolicUtils
using PyPlot
````

````@example 06-1_interaction-picture__PRA2023_107-013706_fig2
# symbolic Hilbert space
hu = FockSpace(:u)
hs = NLevelSpace(:s, 2)
hv = FockSpace(:v)
h = hu ⊗ hs ⊗ hv

# symbolic operators
au_sym = Destroy(h, :a_u, 1)
av_sym = Destroy(h, :a_v, 3)
σ12_sym = Transition(h, :σ, 1, 2, 2)

# symbolic parameters
gu, γ, gv = rnumbers("gu γ gv")

# cascade the SLH elements
G_u = SLH(1, gu * au_sym, 0)
G_s = SLH(1, sqrt(γ) * σ12_sym, 0)
G_v = SLH(1, gv * av_sym, 0)
G_cas = ▷(G_u, G_s, G_v)
nothing # hide
````

````@example 06-1_interaction-picture__PRA2023_107-013706_fig2
H = get_hamiltonian(G_cas)
````

````@example 06-1_interaction-picture__PRA2023_107-013706_fig2
L = get_lindblad(G_cas)[1]
````

Usually we deal with the above derived Hamiltonian and Lindblad. In this example, however, we transform the system into the rotating frame of the virtual cavity-cavity interaction Hamiltonian $H_{uv}$.

````@example 06-1_interaction-picture__PRA2023_107-013706_fig2
H_uv = get_hamiltonian(▷(G_u, G_v))
````

To do so, we first subtract $H_{uv}$ from $H$ and then replace the virtual cavity operators $a_v$ and $a_u$ as described in the [Theory](@ref) section.

````@example 06-1_interaction-picture__PRA2023_107-013706_fig2
H_int_sym_ = simplify(H - H_uv)

# symbolic coefficient matrix $M(t)$
M(i,j) = cnumber("M_{$(i)$(j)}")
a0_ls = [au_sym, av_sym]
la = length(a0_ls)
a_int_ls = [sum(M(i,j)*a0_ls[j] for j=1:la) for i=1:la]

# substitute interaction picture operators
int_dict = Dict([a0_ls; adjoint.(a0_ls)] .=> [a_int_ls; adjoint.(a_int_ls)])
nothing # hide
````

````@example 06-1_interaction-picture__PRA2023_107-013706_fig2
H_int_sym = simplify(substitute_operators(H_int_sym_, int_dict))
````

````@example 06-1_interaction-picture__PRA2023_107-013706_fig2
L_int_sym = simplify(substitute_operators(L, int_dict))
````

The above Hamiltonian and Lindblad operator are the ones in the interaction picture of the virtual cavity interaction.
We define the numerical parameters of the system, calculate the solution for the coefficient matrix $M(t)$ and solve the time evolution of the system.

````@example 06-1_interaction-picture__PRA2023_107-013706_fig2
# numerical parameters
γ_ = 1.0
n_photons = 20

τ = 1 / γ_
t_p = 4 / γ_
u(t) = 1 / (sqrt(τ) * π^(1/4)) * exp(-0.5 * ((t - t_p) / τ)^2)

T_end = 12.0
T = [0:0.005:1;] * T_end

# virtual-cavity couplings for input/output modes
gu_t = u_to_gu(u, T)
gv_t = v_to_gv(u, T) # identical output mode v(t) = u(t)

# interaction-picture coefficient matrix M(t) for u ↔ v
A_uv = interaction_picture_A_2modes(gu_t, gv_t)
M_t = interaction_picture_M(A_uv, T)

# numerical basis
bu = FockBasis(n_photons, n_photons-5)
ba = NLevelBasis(2)
bv = FockBasis(5)
b = bu ⊗ ba ⊗ bv

# constant and time-dependent parameters
dict_p = Dict(γ => γ_)
M_ls = [M(i,j) for i=1:la for j=1:la]
M_t_ls = [t -> M_t(t)[i,j] for i=1:la for j=1:la]
p_t_sym = [gu, gv, M_ls...]
p_t_num = [gu_t, gv_t, M_t_ls...]
dict_p_t = Dict(p_t_sym .=> p_t_num)
nothing # hide
````

````@example 06-1_interaction-picture__PRA2023_107-013706_fig2
H_int_QO = translate(H_int_sym, b; parameter=dict_p, time_parameter=dict_p_t)
L_QO = translate(L_int_sym, b; parameter=dict_p, time_parameter=dict_p_t)

function input_output(t, ρ)
    Ht = H_int_QO(t)
    J = [L_QO(t)]
    return Ht, J, dagger.(J)
end
nothing # hide
````

````@example 06-1_interaction-picture__PRA2023_107-013706_fig2
# time evolution
ψ0 = fockstate(bu, n_photons) ⊗ nlevelstate(ba, 1) ⊗ fockstate(bv, 0)
t, ρt = timeevolution.master_dynamic(T, ψ0, input_output)
nothing # hide

# numerical operators
au = destroy(bu) ⊗ one(ba) ⊗ one(bv)
av = one(bu) ⊗ one(ba) ⊗ destroy(bv)
σee = one(bu) ⊗ transition(ba, 2, 2) ⊗ one(bv)

n_u = real.(expect(au' * au, ρt))
n_v = real.(expect(av' * av, ρt))
P_e = real.(expect(σee, ρt))
nothing # hide
````

We can see that the mean photon number of the initial temporal mode $u$ is reduced by less than two photons.

````@example 06-1_interaction-picture__PRA2023_107-013706_fig2
close("interaction picture fig2") # hide
figure("interaction picture fig2", figsize = (5.6, 4.2))
subplot(2, 1, 1)
plot(T, n_u, label = L"\langle n_u \rangle")
xlabel("tγ")
ylabel("excitations")
ylim(17, 20.2)
grid(true)
legend()

subplot(2, 1, 2)
plot(T, P_e, label = L"\langle \sigma_{ee} \rangle")
plot(T, n_v, label = L"\langle n_v \rangle", ls = "--")
xlabel("tγ")
ylabel("excitations")
ylim(0, 1)
grid(true)
legend()
tight_layout()
gcf()
````

## Package versions

These results were obtained using the following versions:

````@example 06-1_interaction-picture__PRA2023_107-013706_fig2
using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    ["QuantumInputOutput", "QuantumOptics", "SecondQuantizedAlgebra", "SymbolicUtils", "PyPlot"],
    mode = PKGMODE_MANIFEST,
)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

