# # Quantum Pulse Bi-Directional Waveguide
#
# This example mirrors the example `Bi-Directional Waveguide` but uses the numeric SLH struct [`SLHqo`](@ref), to circumvent the symbolic derivation part. 
# Furthermore it drives the system with a *quantum* single-photon pulse via a virtual cavity. 

# As usual, we start by loading the packages and defining the operators and parameters of the system.

using QuantumInputOutput
using QuantumOptics
using Plots

# 

N = 2 # number of quantum dots

## numeric Hilbert space
bu = FockBasis(6) # virtual input cavity 
ba = NLevelBasis(2)
b_qds = tensor([ba for _ = 1:N]...)
b = bu ⊗ b_qds

## QD operators
σ(i, j, k) = embed(b, i+1, transition(ba, j, k))

## input mode operator
a_u = destroy(bu) ⊗ one(b_qds)

## parameters
γ = 1.0
β = 0.9
γRn = fill(γ * β / 2, N)
γLn = fill(γ * β / 2, N)
γ_add = fill(γ * (1 - β), N)
Δn = fill(0.0, N)
ϕn = fill(π/10, max(N - 1, 0))

## quantum pulse u(t) and virtual cavity coupling g_u(t)
σt = 0.8 # pulse with
t0 = 4σt  # pulse peak
Tend = 3t0
T = [0:0.005:1;]*Tend
ΔT = T[2]-T[1]
## u1(t) = sqrt(1 / (σt * √(2π)) * exp(-0.5 * (t - t0)^2 / σt^2)) # hide
u1(t) = 1/(sqrt(σt)*π^(1/4)) * exp(-(t - t0)^2 / (2*σt^2))

gu_t = u_to_gu(u1, T)
nothing # hide

# We use the [`SLHqo`](@ref) to directly use [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) operators and functions to the model the system.

G_u = SLHqo(1, t -> gu_t(t) * a_u, 0*one(b))
G_ϕ(i, j) = SLHqo(exp(1im * ϕn[i]), 0*one(b), 0*one(b))
G_R(i) = SLHqo(1, √(γRn[i]) * σ(i, 1, 2), -Δn[i] * σ(i, 2, 2))
G_L(i) = SLHqo(1, √(γLn[i]) * σ(i, 1, 2), 0*one(b))

## Cascade right-moving channel
## G_R_t = G_d ▷ cascade([G_R(i) ▷ G_ϕ(i, i + 1) for i=1:N-1]...) ▷ G_R(N) # for N > 1 # hide
G_R_t = G_u ▷ G_R(1) ▷ G_ϕ(1, 2) ▷ G_R(2)

## Cascade left-moving channel (reverse order)
## cascade([G_R(i) ▷ G_ϕ(i-1, i) for i=N:-1:2]...) ▷ G_R(1) # hide
G_L_t = G_L(2) ▷ G_ϕ(1, 2) ▷ G_L(1)

G_t = G_R_t ⊞ G_L_t
nothing # hide

# The full Hamiltonian and Lindblad terms are extracted from the final SLH element. Note that as soon as one time-dependent function is involved in a cascade or concatenate, the returned $H$ and $L$ will also be time-dependent. 

H = hamiltonian(G_t)
L = lindblad(G_t)
L_R = L[1]
L_L = L[2]

J_add = [√(γ_add[i]) * σ(i, 1, 2) for i = 1:N]

function input_output(t, ρ)
    Ht = H(t)
    J = [L_R(t), L_L(t), J_add...]
    return Ht, J, dagger.(J)
end
nothing # hide

#

## time evolution
α0 = √(0.1) # √ of total photon number 
ψ0 = coherentstate(bu, α0) ⊗ tensor([nlevelstate(ba, 1) for _ = 1:N]...)
t, ρt = timeevolution.master_dynamic(T, ψ0, input_output)
nothing # hide

#

## transmitted and reflected intensity
I_R = zeros(length(t))
I_L = zeros(length(t))
for (i, ti) in enumerate(t)
    LR = L_R(ti)
    LL = L_L(ti)
    I_R[i] = real(expect(LR' * LR, ρt[i]))
    I_L[i] = real(expect(LL' * LL, ρt[i]))
end
nothing # hide

#

p = plot(t, I_R; label = "Transmission")
plot!(p, t, I_L; label = "Reflection")
plot!(p, t, abs2.(α0*u1.(t)); color = :grey, ls = :dash, label = "Input")
plot!(
    p;
    xlabel = "time",
    ylabel = "intensity",
    legend = :best,
    grid = true,
    size = (500, 320),
)
p

# ## Package versions

# These results were obtained using the following versions:

using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(["QuantumInputOutput", "QuantumOptics", "Plots"], mode = PKGMODE_MANIFEST)
