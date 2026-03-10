```@meta
EditURL = "../../../examples/01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3.jl"
```

# Cavity Scattering of a Single Photon

In this example, we simulate the scattering of a resonant single photon on an empty one-sided cavity. The temporal mode of the light pulse is a Gaussian with width $\sigma$ and the cavity has a decay rate $\gamma$. This system has been studied in [A. Kiilerich, et al., Phys. Rev. Lett. 123, 123604 (2019)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.123604).

We start by loading the needed packages and specifying the model.

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using Plots
using LaTeXStrings
using LinearAlgebra
````

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
# symbolic Hilbert space
hu1 = FockSpace(:u1)
hc1 = FockSpace(:c1)
hv1 = FockSpace(:v1)
h = hu1 ⊗ hc1 ⊗ hv1

# symbolic operators
au = Destroy(h, :a_u, 1)
c = Destroy(h, :c, 2)
av = Destroy(h, :a_v, 3)

# symbolic parameters
gu, Δ, γ = rnumbers("g_u Δ γ")
gv = cnumber("g_v")
nothing # hide
````

We use the symbolic operators and parameters to define the SLH triples and cascade them to obtain the Hamiltonian and Lindblad for the system.

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
G_u = SLH(1, gu*au, 0) # input cavity
G_c = SLH(1, √(γ)*c, Δ*c'c) # system cavity
G_v = SLH(1, gv*av, 0) # output cavity

G_cas = ▷(G_u, G_c, G_v)
nothing # hide
````

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
H = get_hamiltonian(G_cas)
````

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
L = get_lindblad(G_cas)[1] # only one Lindblad term in this example
````

To solve the dynamics of the system we translate the symbolic expressions into numeric operators (matrices) of QuantumOptics.jl. To do so, we define the numerical parameters and operator basis.

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
# numerical parameters
γ_ = 1.0
Δ_ = 0.0

p_sym = [γ, Δ, gv]
p_num = [γ_, Δ_, 0] # gv=0
dict_p = Dict(p_sym .=> p_num);

# Gaussian input mode
σ = 1/γ_
T_end = 12σ
# u1(t) = sqrt(1/(σ*√(2π))*exp( -0.5*(t - 4T_p)^2/σ^2 )) # hide
u1(t) = 1/(sqrt(σ)*π^(1/4)) * exp(-(t - 4σ)^2 / (2*σ^2))
T = [0:0.002:1;]*T_end
ΔT = T[2] - T[1]

gu_t = u_to_gu(u1, T)
dict_p_t = Dict(gu => gu_t)

# numeric bases
bu1 = FockBasis(1)
bc1 = FockBasis(1)
bv1 = FockBasis(1)
b = bu1 ⊗ bc1 ⊗ bv1
nothing # hide
````

We now use the function [`translate_qo`](@ref) to create the numeric operators. If the kwarg `time_parameter` is provided the created operator is a time-dependent function.

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
H_QO = translate_qo(H, b; parameter = dict_p, time_parameter = dict_p_t)
L_QO = translate_qo(L, b; parameter = dict_p, time_parameter = dict_p_t)
nothing # hide
````

To solve the dynamics we use the QuantumOptics.jl function `timeevolution.master_dynamic`.

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
# time-dependent function for timeevolution.master_dynamic that returns H(t), J(t) and Jd(t)
function input_output_1(t, ρ)
    H = H_QO(t)
    J = [L_QO(t)]
    return H, J, dagger.(J)
end;

# initial state
ψ0 = fockstate(bu1, 1) ⊗ fockstate(bc1, 0) ⊗ fockstate(bv1, 0)

# time evolution
t_, ρt = timeevolution.master_dynamic(T, ψ0, input_output_1)
nothing # hide
````

We create the desired numerical operators to calculate expectation values.

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
au_qo = translate_qo(au, b)
c_qo = translate_qo(c, b)
av_qo = translate_qo(av, b)

n_c_t = real.(expect(c_qo'*c_qo, ρt))
n_u1_t = real.(expect(au_qo'*au_qo, ρt))
nothing # hide
````

In order to determine suitable temporal output modes we calculate the two-time autocorrelation function $g^{(1)}(t_1,t_2) = \langle L_s^\dagger(t_1) L_s(t_2) \rangle$ and diagonalize the matrix to obtain the eigenvalues with the corresponding eigenvectors. The eigenvalues correspond to the mean photon number $n_i$ in the corresponding temporal eigenvector mode $v_i$.

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
Ls(t) = gu_t(t)*au_qo + √(γ_)*c_qo
g1_m = two_time_corr_matrix(T, ρt, input_output_1, Ls)
nothing # hide
````

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
p = heatmap(
    T,
    T,
    real.(g1_m);
    c = :inferno,
    xlabel = L"\gamma t_2",
    ylabel = L"\gamma t_1",
    colorbar_title = L"g^{(1)}(t_1,t_2)",
    size = (400, 350),
)
p
````

The eigenvalues and corresponding eigenvectors are sorted in ascending order, which means the last eigenvalue corresponds to the highest populated temporal mode.

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
F = eigen(g1_m)
n_avg = round.(real.(F.values)*ΔT; digits = 3)
modes = F.vectors
v_mode = (modes[:, end]) / sqrt(ΔT)

@show n_avg[(end-1):end]
nothing # hide
````

We use now the mode with the highest mean photon number as our out-mode to determine its quantum state.

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
p_sym_2 = [γ, Δ]
p_num_2 = [γ_, Δ_]
dict_p_2 = Dict(p_sym_2 .=> p_num_2)

# time-dependent coupling for the output mode $v(t)$
gv_t = v_to_gv(v_mode, T)

dict_p_t_2 = Dict([gu, gv] .=> [gu_t, gv_t])

H_QO_2 = translate_qo(H, b; parameter = dict_p_2, time_parameter = dict_p_t_2)
L_QO_2 = translate_qo(L, b; parameter = dict_p_2, time_parameter = dict_p_t_2)
function input_output_2(t, ρ)
    H = H_QO_2(t)
    J = [L_QO_2(t)]
    return H, J, dagger.(J)
end;

# time evolution for the system including the output cavity
t_2, ρt_2 = timeevolution.master_dynamic(T, ψ0, input_output_2)

# mean photon number in the output mode
n_v1_t = real.(expect(av_qo'*av_qo, ρt_2))
nothing # hide
````

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
p1 = plot(T, u1.(T); ls = :dash, label = "u", color = :red)
plot!(p1, T, u1.(T); fillrange = 0, fillalpha = 0.5, color = :red, label = "")
plot!(p1, T, real.(v_mode); color = :blue, label = "v")
plot!(p1, T, real.(v_mode); fillrange = 0, fillalpha = 0.5, color = :blue, label = "")
plot!(
    p1;
    xlims = (0, 12),
    ylims = (-0.8, 0.8),
    yticks = [-0.8, 0, 0.8],
    ylabel = "modes",
    legend = :best,
)

p1r = twinx(p1)
plot!(p1r, T, abs2.(gu_t.(T)); color = :red, label = "")
plot!(p1r, T[3:end], abs2.(gv_t.(T))[3:end]; color = :blue, ls = :dash, label = "")
plot!(p1r; xlims = (0, 12), ylims = (0, 8), ylabel = "Rates (γ)")

p2 = plot(T, n_u1_t; label = L"\langle a^\dagger a \rangle_u", color = :red)
plot!(p2, T, n_c_t; label = L"\langle c^\dagger c \rangle", color = :green)
plot!(p2, T, n_v1_t; label = L"\langle a^\dagger a \rangle_v", color = :blue)
plot!(
    p2;
    xlims = (0, 12),
    ylims = (0, 1),
    xlabel = "time (1/γ)",
    ylabel = "Exciations",
    legend = :best,
)

plot(p1, p2; layout = (2, 1), size = (600, 550))
````

## Cavity with phase noise

We slightly adapt the above example by assuming the initial pulse to be in a coherent state and adding phase noise to the cavity. This results in scattering into multiple modes.

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
# new basis of the system
bu1_3 = FockBasis(12)
bc1_3 = FockBasis(6)
bv1_3 = FockBasis(6)
b_3 = tensor(bu1_3, bc1_3, bv1_3)

# new operators of the system
au_3 = embed(b_3, 1, destroy(bu1_3))
c_3 = embed(b_3, 2, destroy(bc1_3))
av_3 = embed(b_3, 3, destroy(bv1_3))
cdc_3 = c_3'c_3

# we use the same Hamiltonian as before but add a depasing term to the dissipation
H_QO_3 = translate_qo(H, b_3; parameter = dict_p_2, time_parameter = dict_p_t_2)
L_QO_3 = translate_qo(L, b_3; parameter = dict_p_2, time_parameter = dict_p_t_2)
function input_output_3(t, ρ)
    H = H_QO_3(t)
    J = [L_QO_3(t), √(γ_)*cdc_3]
    return H, J, dagger.(J)
end;
nothing #hide
````

Due to the larger Hilbert space the time evolution takes a few seconds.

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
ψ0_3 = coherentstate(bu1_3, 2) ⊗ fockstate(bc1_3, 0) ⊗ fockstate(bv1_3, 0)
t_3, ρt_3 = timeevolution.master_dynamic(T, ψ0_3, input_output_3)
nothing # hide
````

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
L0(t) = √(γ_)*c_3 + gu_t(t)*au_3 + gv_t(t)*av_3
I_out = [real(expect(dagger(L0(t_3[i]))*L0(t_3[i]), ρt_3[i])) for i = 1:length(t_3)]

n_u1_t_3 = real.(expect(au'*au, ρt_3))
n_v1_t_3 = real.(expect(av'*av, ρt_3))
nothing # hide
````

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
p = plot(t_3, n_u1_t_3; label = L"\langle a^\dagger a \rangle_u", color = :red)
plot!(p, t_3, n_v1_t_3; label = L"\langle a^\dagger a \rangle_v", color = :blue, ls = :dash)
plot!(p, t_3, I_out; label = L"I_{out}", color = :black, ls = :dot)
plot!(
    p;
    xlims = (0, 12),
    ylims = (0, 4),
    xlabel = "time (1/γ)",
    ylabel = "expectation values",
    legend = :best,
    size = (600, 300),
)
p
````

## Package versions

These results were obtained using the following versions:

````@example 01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3
using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    ["QuantumInputOutput", "SecondQuantizedAlgebra", "QuantumOptics", "Plots", "LaTeXStrings"],
    mode = PKGMODE_MANIFEST,
)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

