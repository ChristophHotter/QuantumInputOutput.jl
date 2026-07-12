# # Bi-Directional Waveguide
#
# This example constructs an SLH model for `N=2` quantum dots coupled to a
# bi-directional waveguide. A coherent input pulse enters from the left (right-moving
# mode), and we compute the time evolution of the transmitted and reflected intensities. 

using QuantumInputOutput
using SecondQuantizedAlgebra
using Symbolics: Symbolics
using QuantumOptics
using QuantumOpticsBase: dagger
using Plots
using LaTeXStrings

# 

N = 2 # number of quantum dots

## symbolic Hilbert space
ha(i) = NLevelSpace(Symbol("a$(i)"), 2)
h = tensor([ha(i) for i = 1:N]...)

## symbolic operators
σ(α, i, j) = Transition(h, Symbol("σ_$(α)"), i, j, α)

## symbolic parameters
γR(i) = Symbolics.variable(Symbol("γ^{($(i))}_R"); T = Real) # right-moving decay rate
γL(i) = Symbolics.variable(Symbol("γ^{($(i))}_L"); T = Real) # left-moving decay rate
Δ(i) = Symbolics.variable(Symbol("Δ_{$(i)}"); T = Real) # detuning
ϕ(i, j) = Symbolics.variable(Symbol("ϕ_{$(i)$(j)}"); T = Real) # phase between QD-i and QD-j
Ein = Symbolics.variable(Symbol("E_{in}"); T = Real) # coherent drive in the right-moving input
nothing # hide

# We use the symbolic operators and parameters to define the SLH triples, cascade the left and right moving channels, and concatenate them to obtain the Hamiltonian and Lindblad for the system. 

G_d = SLH(1, Ein, 0) # coherent drive in the right-moving input
G_ϕ(i, j) = SLH(exp(1im * ϕ(i, j)), 0, 0) # phase shift
G_R(i) = SLH(1, √(γR(i)) * σ(i, 1, 2), -Δ(i) * σ(i, 2, 2)) # right-moving decay
G_L(i) = SLH(1, √(γL(i)) * σ(i, 1, 2), 0) # left-moving decay


## Cascade right-moving channel
## G_R_t = G_d ▷ cascade([G_R(i) ▷ G_ϕ(i, i + 1) for i=1:N-1]...) ▷ G_R(N) # for N > 1 # hide
G_R_t = G_d ▷ G_R(1) ▷ G_ϕ(1, 2) ▷ G_R(2)

## Cascade left-moving channel (reverse order)
## cascade([G_R(i) ▷ G_ϕ(i-1, i) for i=N:-1:2]...) ▷ G_R(1) # hide
G_L_t = G_L(2) ▷ G_ϕ(1, 2) ▷ G_L(1)

## Concatenate both channels
G_t = G_R_t ⊞ G_L_t
nothing # hide

H = hamiltonian(G_t)

#

L = lindblad(G_t)
L_R = L[1]

#

L_L = L[2]

# Note that this Hamiltonian and Lindblad terms (without the drive) describe the collective decay of the quantum dots. 

# Next, the numerical parameters and functions of the system are defined, and we translate the symbolic expression to [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) operators to numerically solve the dynamics. 

γ_ = 1.0
β = 0.9 # waveguide coupling fraction
γRn = fill(γ_ * β / 2, N)
γLn = fill(γ_ * β / 2, N)
γ_add = fill(γ_ * (1-β), N) # free space decay
Δn = fill(0.0, N)
ϕn = fill(π/10, max(N - 1, 0))

σt = 0.8 # pulse with
α0 = √(0.1) # √ of total photon number
t0 = 4σt # pulse peak
Tend = 3t0
## u1(t) = sqrt(1 / (σt * √(2π)) * exp(-0.5 * (t - t0)^2 / σt^2)) # hide
u1(t) = 1/(sqrt(σt)*π^(1/4)) * exp(-(t - t0)^2 / (2*σt^2))
Ein_t(t) = α0*u1(t)

p_sym = [
    [γR(i) for i = 1:N];
    [γL(i) for i = 1:N];
    [Δ(i) for i = 1:N];
    [ϕ(i, i + 1) for i = 1:(N-1)]
]
p_num = [γRn; γLn; Δn; ϕn]
dict_p = Dict(p_sym .=> p_num)
dict_p_t = Dict(Ein => Ein_t)
nothing # hide

# 

## numeric bases
ba = NLevelBasis(2)
b = tensor([ba for i = 1:N]...)

H_QO = to_numeric(H, b; parameter = dict_p, time_parameter = dict_p_t)
L_R_QO = to_numeric(L_R, b; parameter = dict_p, time_parameter = dict_p_t)
L_L_QO = to_numeric(L_L, b; parameter = dict_p, time_parameter = dict_p_t)

σ_qo(α, i, j) = to_numeric(σ(α, i, j), b)
J_add = [√(γ_add[i])*σ_qo(i, 1, 2) for i = 1:N]

function input_output(t, ρ)
    Ht = H_QO(t)
    J = [L_R_QO(t), L_L_QO(t), J_add...]
    return Ht, J, dagger.(J)
end
nothing # hide

#

## time evolution
T = [0:0.005:1;]*Tend
ψ0 = tensor([nlevelstate(ba, 1) for _ = 1:N]...)
t, ρt = timeevolution.master_dynamic(T, ψ0, input_output)
nothing # hide

#

## transmitted and reflected intensity
I_R = zeros(length(t))
I_L = zeros(length(t))

for (i, ti) in enumerate(t)
    LR = L_R_QO(ti)
    LL = L_L_QO(ti)
    I_R[i] = real(expect(LR'LR, ρt[i]))
    I_L[i] = real(expect(LL'LL, ρt[i]))
end
nothing # hide

#

p = plot(t, I_R; label = "Transmission")
plot!(p, t, I_L; label = "Reflection")
plot!(p, t, abs2.(Ein_t.(t)); color = :grey, ls = :dash, label = "Input")
plot!(
    p;
    xlabel = "time",
    ylabel = "intensity",
    legend = :best,
    grid = true,
    size = (500, 320),
)
p

# ## Quantum regression theorem

# In the following, we calculate the two-time correlation function $G^{(2)}(t_1,t_2)$ for the transmitted and reflected pulse via the quantum regression theorem.

## two-time correlation function G2(t1, t2)
lT = length(T)
G2 = zeros(lT, lT) # transmission
G2_ref = zeros(lT, lT) # reflection

L0(t) = L_R_QO(t)
L0_dag(t) = dagger(L0(t))
L0_ref(t) = L_L_QO(t)
L0_ref_dag(t) = dagger(L0_ref(t))

for it1 = 1:(lT-1)
    ρ_t1 = ρt[it1]

    t_2, ρ_2 = timeevolution.master_dynamic(
        T[it1:end],
        L0(T[it1]) * ρ_t1 * L0_dag(T[it1]),
        input_output,
    )

    ## transmission
    G2_ls = real.([expect(L0_dag(t_2[j]) * L0(t_2[j]), ρ_2[j]) for j = 1:length(t_2)])
    G2[it1, it1:end] = G2_ls
    G2[it1:end, it1] = G2_ls

    t_2_r, ρ_2_r = timeevolution.master_dynamic(
        T[it1:end],
        L0_ref(T[it1]) * ρ_t1 * L0_ref_dag(T[it1]),
        input_output,
    )

    ## reflection
    G2_ls_r = real.([
        expect(L0_ref_dag(t_2_r[j]) * L0_ref(t_2_r[j]), ρ_2_r[j]) for j = 1:length(t_2_r)
    ])
    G2_ref[it1, it1:end] = G2_ls_r
    G2_ref[it1:end, it1] = G2_ls_r
end
nothing # hide

#

p_ref = heatmap(
    T,
    T,
    G2_ref' / maximum(G2_ref);
    c = :inferno,
    title = "reflection",
    xlabel = L"t_1",
    ylabel = L"t_2",
    colorbar_title = L"G^{(2)}(t_1, t_2)[a.u.]",
)
p_trans = heatmap(
    T,
    T,
    G2' / maximum(G2);
    c = :inferno,
    title = "transmission",
    xlabel = L"t_1",
    ylabel = L"t_2",
    colorbar_title = L"G^{(2)}(t_1, t_2)[a.u.]",
)
plot(p_ref, p_trans; layout = (1, 2), size = (700, 300))

# ## Package versions

# These results were obtained using the following versions:

using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    [
        "QuantumInputOutput",
        "SecondQuantizedAlgebra",
        "QuantumOptics",
        "Plots",
        "LaTeXStrings",
    ],
    mode = PKGMODE_MANIFEST,
)
