```@meta
EditURL = "../../../examples/03-1_beam-combiner__PRA2023_107-023715_fig2-fig3.jl"
```

# Perfect Splitting and Combining of a Two-photon Pulse

In this example, we simulate the perfect splitting of a two-photon pulse into two orthogonal temporal modes with one photon each.
We then show the reverse process to combine the two orthogonal photons into a single temporal mode with two photons.
The system is described in [M. Lund , et al., Phys. Rev. A 107, 023715 (2023)](https://doi.org/10.1103/PhysRevA.107.023715).

As usual, we start by loading the packages and defining the symbolic operators and parameters.

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using PyPlot
using LinearAlgebra
using DataInterpolations
````

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
# symbolic Hilbert space
hu2 = FockSpace(:u2)
hu1 = FockSpace(:u1)
hs1 = NLevelSpace(:atom, 2)
hv1 = FockSpace(:v1)
h = hu2 ⊗ hu1 ⊗ hs1 ⊗ hv1

# symbolic operators
au2 = Destroy(h, :au_2, 1)
au1 = Destroy(h, :au_1, 2)
σ(i, j) = Transition(h, :σ, i, j, 3)
av1 = Destroy(h, :av_2, 4)

# symbolic parameters
@rnumbers γ Δ
gu1, gu2, gv1 = cnumbers("gu_1 gu_2 gv_1");
nothing #hide
````

We use the symbolic operators and parameters to define the SLH triples and cascade them to obtain the Hamiltonian and Lindblad for the system.

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
G_u2 = SLH(1, gu2*au2, 0) # input cavity 2
G_u1 = SLH(1, gu1*au1, 0) # input cavity 1
G_a = SLH(1, √(γ)*σ(1, 2), Δ*σ(2, 2)) # scattering atom
G_v1 = SLH(1, gv1*av1, 0) # output cavity 1

G_cas = cascade(G_u2, G_u1, G_a, G_v1)
nothing # hide
````

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
H = get_hamiltonian(G_cas)
````

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
L = get_lindblad(G_cas)[1] # only one Lindblad in this example
````

Next, the numerical parameters and functions of the system are defined.

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
γ_ = 1.0
Δ_ = 0.0

p_sym = [γ, Δ, gu2, gv1]
p_num = [γ_, Δ_, 0, 0]
dict_p = Dict(p_sym .=> p_num)

# Gaussian input pulse
τ = 0.38;
tp = 4/γ_
u(t) = 1/(sqrt(τ)*π^(1/4)) * exp(-(t - tp)^2 / (2*τ^2))
T = [0:0.002:1;]*20
ΔT = T[2] - T[1]

gu_ = u_to_gu(u, T)
dict_p_t = Dict(gu1 => gu_)
nothing # hide
````

We translate the symbolic expressions to numerical operators and solve the time-dependent master equation with QuantumOptics.jl.

To obtain the output modes we do not use the second input mode and the output mode cavity.
However, to keep the example short we include them already from the beginning since they are needed later.
To perform time consuming parameter scans one should merely use the necessary Hilbert spaces. In this case, this would correspond to one input cavity and the two-level system.
The kwarg `operators` of the function [translate](@ref) provides a convenient way to use predefined numerical operators, see the example `Two-sided Cavity with Atom`.

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
# numeric bases
bu2 = FockBasis(2)
bu1 = FockBasis(2)
bs1 = NLevelBasis(2)
bv1 = FockBasis(2)
b = bu2 ⊗ bu1 ⊗ bs1 ⊗ bv1;

H_QO = translate(H, b; parameter = dict_p, time_parameter = dict_p_t)
L_QO = translate(L, b; parameter = dict_p, time_parameter = dict_p_t)
function input_output(t, ρ)
    H = H_QO(t)
    J = [L_QO(t)]
    return H, J, dagger.(J)
end
nothing # hide
````

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
# time evolution
ψ0 = fockstate(bu2, 0) ⊗ fockstate(bu1, 2) ⊗ nlevelstate(bs1, 1) ⊗ fockstate(bv1, 0)
t_, ρt = timeevolution.master_dynamic(T, ψ0, input_output)
nothing # hide
````

Now we analyze the output modes with the two-time autocorrelation function $g^{(1)}(t_1,t_2) = \langle L_s^\dagger(t_1) L_s(t_2) \rangle$.

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
au1_qo = translate(au1, b)
σ_qo(i, j) = translate(σ(i, j), b)

Ls(t) = (gu_(t))'*au1_qo + √(γ_)*σ_qo(1, 2)
g1_m = two_time_corr_matrix(T, ρt, input_output, Ls);

close("g1(t1,t2) matrix") # hide
figure("g1(t1,t2) matrix", figsize = (4.5, 3.5))
pcolormesh(T, T, real.(g1_m), cmap = "inferno")
xlabel(L"\gamma t_2")
ylabel(L"\gamma t_1")
colorbar(label = L"g^{(1)}(t_1,t_2)")
tight_layout()
gcf()
````

The eigenvalues correspond to the mean photon number $n_i$ in the corresponding temporal eigenvector mode $v_i$. We find two modes with a mean photon number of one.

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
F = eigen(g1_m)
n_avg = real.(F.values)*ΔT
modes = F.vectors
v1_mode = (modes[:, end]) / √(ΔT)
v2_mode = (modes[:, end-1]) / √(ΔT)

@show n_avg[(end-1):end]
nothing # hide
````

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
close("modes") # hide
figure("modes")
plot(t_, -real.(v1_mode), color = "black")
plot(t_, real.(v2_mode), color = "red", ls = "--")
xlabel("time (1/γ)")
ylabel("output mode")
gcf()
````

As described in the paper, we can define a rotated basis in which the two modes are not entangled and equally populated by a single photon Fock state.
In the following, we define these rotated modes and use them to combine two single photons into a two photon Fock state. The temporal output mode of this two photon Fock state is the same as the previous input mode which separated the two single photons before.

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
v1_p = 1/√(2) * (v1_mode - v2_mode)
v2_p = 1/√(2) * (v1_mode + v2_mode)

# new output mode = old input mode
v1_new(t) = (u(T[end]-t))'

# new input modes
u1_new = conj.(reverse(v1_p))
u2_new = conj.(reverse(v2_p))
nothing # hide
````

The pulse from input cavity $u_2$ is scattered on the cavity $u_1$. This distortion needs to be taken into account for the coupling of $u_2$,
which is done with the function [`u_eff`](@ref).
The coupling of the $u_1$ cavity needs no adaptation, since it directly couples to the two-level system.

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
gu1_ = u_to_gu(u1_new, T)

u_new_data = [u1_new, u2_new]
u_new_fct = [LinearInterpolation(u, T) for u in u_new_data]
````

effective u2 mode and corresponding coupling

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
u2_for_gu2 = u_eff(u_new_fct, T, 2)
gu2_ = u_to_gu(u2_for_gu2, T)

# coupling of the output mode
gv1_ = v_to_gv(v1_new, T)

# dictionary for the time-dependent functions
g_sym = [gu1, gu2, gv1]
g_num = [gu1_, gu2_, gv1_]
dict_p_t_out = Dict(g_sym .=> g_num)

# dictionary for the constant parameters
p_sym_out = [γ, Δ]
p_num_out = [γ_, Δ_]
dict_p_out = Dict(p_sym_out .=> p_num_out)
nothing # hide
````

The time-dependent couplings are used to define the numeric Hamiltonian and Lindblad term, and then solve the dynamics of the system.

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
H_QO_2 = translate(H, b; parameter = dict_p_out, time_parameter = dict_p_t_out)
L_QO_2 = translate(L, b; parameter = dict_p_out, time_parameter = dict_p_t_out)
function input_output_2(t, ρ)
    H = H_QO_2(t)
    J = [L_QO_2(t)]
    return H, J, dagger.(J)
end

# time evolution
ψ0_out = fockstate(bu2, 1) ⊗ fockstate(bu1, 1) ⊗ nlevelstate(bs1, 1) ⊗ fockstate(bv1, 0)
t_2, ρt_2 = timeevolution.master_dynamic(T, ψ0_out, input_output_2)
nothing # hide
````

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
nu1_t_comb = real(expect(au1'*au1, ρt_2))
nu2_t_comb = real(expect(au2'*au2, ρt_2))
nv1_t_comb = real(expect(av1'*av1, ρt_2))
s22_t_comb = real(expect(σ(2, 2), ρt_2))
nothing # hide
````

We can see that the two single photons combine to a two photon Fock-state in one temporal mode.

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
close("beam combiner") # hide
figure("beam combiner")
plot(T, nu2_t_comb, color = "red", ls = "--", label = L"\langle a^\dagger a \rangle_{u_2}")
plot(T, nu1_t_comb, color = "blue", ls = "-.", label = L"\langle a^\dagger a \rangle_{u_1}")
plot(T, s22_t_comb, color = "black", ls = "dotted", label = L"\langle \sigma^{22} \rangle")
plot(T, nv1_t_comb, color = "green", ls = "-", label = L"\langle a^\dagger a \rangle_{v_1}")
ylim(0, 2)
xlim(10, 18)
xlabel("time (1/γ)")
ylabel("Mean Excitation")
legend()
tight_layout()
gcf()
````

## Package versions

These results were obtained using the following versions:

````@example 03-1_beam-combiner__PRA2023_107-023715_fig2-fig3
using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    [
        "QuantumInputOutput",
        "SecondQuantizedAlgebra",
        "QuantumOptics",
        "PyPlot",
        "DataInterpolations",
    ],
    mode = PKGMODE_MANIFEST,
)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

