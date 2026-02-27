```@meta
EditURL = "../../../examples/05-1_N-QDs_bidirectional-waveguide_coherent-pulse.jl"
```

# Bi-Directional Waveguide

This example constructs an SLH model for `N=2` quantum dots coupled to a
bi-directional waveguide. A coherent input pulse enters from the left (right-moving
mode), and we compute the time evolution of the transmitted and reflected intensities.

````@example 05-1_N-QDs_bidirectional-waveguide_coherent-pulse
using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using PyPlot
````

````@example 05-1_N-QDs_bidirectional-waveguide_coherent-pulse
N = 2 # number of quantum dots

# symbolic Hilbert space
ha(i) = NLevelSpace("a$(i)", 2)
h = tensor([ha(i) for i = 1:N]...)

# symbolic operators
σ(α, i, j) = Transition(h, "σ_$(α)", i, j, α)

# symbolic parameters
γR(i) = rnumber("γ^{($(i))}_R") # right-moving decay rate
γL(i) = rnumber("γ^{($(i))}_L") # left-moving decay rate
Δ(i) = rnumber("Δ_{$(i)}") # detuning
ϕ(i, j) = rnumber("ϕ_{$(i)$(j)}") # phase between QD-i and QD-j
Ein = rnumber("E_{in}") # coherent drive in the right-moving input
nothing # hide
````

We use the symbolic operators and paramters to define the SLH triples, cascade the left and right moving channels, and concatenate them to obtain the Hamiltonian and Lindblad for the system.

````@example 05-1_N-QDs_bidirectional-waveguide_coherent-pulse
G_d = SLH(1, Ein, 0) # coherent drive in the right-moving input
G_ϕ(i, j) = SLH(exp(1im * ϕ(i, j)), 0, 0) # phase shift
G_R(i) = SLH(1, √(γR(i)) * σ(i, 1, 2), -Δ(i) * σ(i, 2, 2)) # right-moving decay
G_L(i) = SLH(1, √(γL(i)) * σ(i, 1, 2), 0) # left-moving decay

# Cascade right-moving channel
G_R_t = cascade(G_d, [G_R(i) for i=1:N]...)

# Cascade left-moving channel (reverse order)
G_L_t = cascade([G_L(i) for i=N:-1:1]...)

# Concatenate both channels
G_t = G_R_t ⊞ G_L_t
nothing # hide

H = get_hamiltonian(G_t)
````

````@example 05-1_N-QDs_bidirectional-waveguide_coherent-pulse
L = get_lindblad(G_t)
L_R = L[1]
````

````@example 05-1_N-QDs_bidirectional-waveguide_coherent-pulse
L_L = L[2]
````

Next, the numerical parameters and functions of the system are defined, and we translate the symbolic expression to `QuantumOptics` operators to numerically solve the dynamics.

````@example 05-1_N-QDs_bidirectional-waveguide_coherent-pulse
γ_ = 1.0
β = 0.9 # waveguide coupling fraction
γRn = fill(γ_ * β / 2, N)
γLn = fill(γ_ * β / 2, N)
γ_add = fill(γ_ * (1-β), N) # free space decay
Δn = fill(0.0, N)
ϕn = fill(0.25 * π, max(N - 1, 0))

σt = 1.5 # pulse with
α0 = √(0.1) # √ of total photon number
Ω0 = α0/(π^(1/4)*√(σt))
t0 = 3σt
Tend = 3t0
Ein_t(t) = Ω0 * exp(-0.5 * (t - t0)^2 / σt^2)

p_sym = [ [γR(i) for i = 1:N];
          [γL(i) for i = 1:N];
          [Δ(i) for i = 1:N];
          [ϕ(i, i + 1) for i = 1:N-1] ]
p_num = [ γRn; γLn; Δn; ϕn ]
dict_p = Dict(p_sym .=> p_num)
dict_p_t = Dict(Ein => Ein_t)
nothing # hide
````

````@example 05-1_N-QDs_bidirectional-waveguide_coherent-pulse
# numeric bases
ba = NLevelBasis(2)
b = tensor([ba for i = 1:N]...)

H_QO = translate(H, b; parameter=dict_p, time_parameter=dict_p_t)
L_R_QO = translate(L_R, b; parameter=dict_p, time_parameter=dict_p_t)
L_L_QO = translate(L_L, b; parameter=dict_p, time_parameter=dict_p_t)

σ_qo(α,i,j) = translate(σ(α,i,j), b)
J_add = [√(γ_add[i])*σ_qo(i,1,2) for i=1:N]

function input_output(t, ρ)
    Ht = H_QO(t)
    J = [L_R_QO(t), L_L_QO(t), J_add...]
    return Ht, J, dagger.(J)
end
nothing # hide
````

````@example 05-1_N-QDs_bidirectional-waveguide_coherent-pulse
# time evolution
T = collect(0.0:0.05:10.0)
ψ0 = tensor([nlevelstate(ba, 1) for _ = 1:N]...)
t, ρt = timeevolution.master_dynamic(T, ψ0, input_output)
nothing # hide
````

````@example 05-1_N-QDs_bidirectional-waveguide_coherent-pulse
# transmitted and reflected intensity
I_R = zeros(length(t))
I_L = zeros(length(t))

for (i, ti) in enumerate(t)
    LR = L_R_QO(ti)
    LL = L_L_QO(ti)
    I_R[i] = real(expect(LR'LR, ρt[i]))
    I_L[i] = real(expect(LL'LL, ρt[i]))
end
nothing # hide
````

````@example 05-1_N-QDs_bidirectional-waveguide_coherent-pulse
close("time evolution")
figure("time evolution", figsize = (5, 3.2))
plot(t, I_R, label = "Transmission")
plot(t, I_L, label = "Reflection")
plot(t, abs2.(Ein_t.(t)), color="grey", ls="--", label = "Input")
xlabel("time")
ylabel("intensity")
legend()
grid(true)
tight_layout()
gcf()
````

TODO: simulate g2

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

