# # Mean-field Two-sided Cavity

# Here we show how to solve the dynamics of the example `Two-sided Cavity with Atoms` in the Heisenberg picture with a higher-order mean-field approach (cumulant expansion), which is done with the package [QuantumCumulants.jl](https://github.com/qojulia/QuantumCumulants.jl).    

using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumCumulants
using ModelingToolkitBase
using OrdinaryDiffEq
using QuantumOpticsBase
using QuantumOpticsBase: dagger
using Plots

#

@variables E::Real κ_L::Real κ_R::Real Δ::Real g::Real γ::Real
Natoms = 2

hc = FockSpace(:cavity)
ha(i) = NLevelSpace(Symbol("a_$i"), 2)
h = hc ⊗ tensor([ha(i) for i = 1:Natoms]...);

a = Destroy(h, :a, 1) # cavity 
σ(α, i, j) = Transition(h, Symbol("σ_$(α)"), i, j, 1+α) # two-level atom α
∑σ(i, j) = sum(σ(α, i, j) for α = 1:Natoms) # collective atomic operator
nothing # hide 

# ## Empty two-sided cavity

# We couple a classical drive into the cavity through the left mirror $(\kappa_L)$. The decay through the right mirror can be added in several ways: with the `concatenate` rule, including it already in the initial cavity SLH triple or by simply including the decay term to the Lindblad by hand. In this example, we use the first option. 

G_d = SLH(1, E, 0) # classical drive
H_cavity = -Δ*a'a
G_c_L = SLH(1, [√(κ_L)*a], H_cavity)

G_cav_L_drive = G_d ▷ G_c_L
G_c_R = SLH(1, [√(κ_R)*a], 0)

G_cav_L_R_drive = G_cav_L_drive ⊞ G_c_R
nothing # hide 

# Note that one needs to be careful to not double-count the Hamiltonian terms with the concatenation rule. 

H1 = hamiltonian(G_cav_L_R_drive)

# 

J1_L = jump_operator(G_cav_L_R_drive)[1]

# 

J1_R = jump_operator(G_cav_L_R_drive)[2]

# The typical cavity drive-term $\sqrt{\kappa_L} E (a^\dagger + a)$ is a combination of Hamiltonian term and Lindblad. 
# We use the function `meanfield` to obtain the equation for the intra-cavity field, which leads to a closed set of equations in this particular case.  

eqs_a = meanfield([a], H1, [J1_L, J1_R])
complete!(eqs_a)
nothing # hide

# ```math
# \begin{align} \frac{d}{dt} \langle a\rangle &= 1 i \Delta \langle a\rangle -1.0 \sqrt{\kappa{L}} E -0.5 \left( \left( \sqrt{\kappa{L}} \right)^{2} + \left( \sqrt{\kappa_{R}} \right)^{2} \right) \langle a\rangle \end{align}
# ```

# We defined the numerical parameters and the initial state of the system, create the ODE problem and solve the dynamics.

## numerical parameters
En = 0.5
κ_Rn = 1.0
κ_Ln = 1.0
Δn = 0.0
Δn_ls = [-5.0:0.1:5.0;];
lΔ=length(Δn_ls)

p_sym = [E, κ_R, κ_L, Δ]
p_num = [En, κ_Rn, κ_Ln, Δn]
dict_p1 = Dict(p_sym .=> p_num)

## initial state
u0 = [0.0im]
nothing # hide

#

## numerical system
T = [0:0.01:1;]*20
sys_a = mtkcompile(System(eqs_a; name = :sys_a))
dict = merge(Dict(p_sym .=> p_num), initial_values(eqs_a, u0))
prob_a = ODEProblem(sys_a, dict, (0.0, T[end]))
nothing # hide

#

sol_a = solve(prob_a, Tsit5(); saveat = T)
nothing # hide

#

n_cavity = abs2.(get_solution(sol_a, a, eqs_a).(T))
n_ref = abs2.(√(κ_Ln) .* get_solution(sol_a, a, eqs_a).(T) .+ En)
n_trans = abs2.(√(κ_Rn) .* get_solution(sol_a, a, eqs_a).(T))
nothing # hide

#

p1 = plot(T, n_cavity; label = "", xlabel = "t", ylabel = "cavity photons", grid = true)
p2 = plot(T, n_ref; label = "reflection")
plot!(p2, T, n_trans; label = "transmission", ls = :dash)
plot!(p2; xlabel = "t", ylabel = "intensity rate", grid = true, legend = :best)
plot(p1, p2; layout = (2, 1), size = (600, 500))

# Now we scan the laser-cavity detuning $\Delta$ to plot the transmission and reflection spectrum. 

dict_p_Δ(Δn) = merge(Dict(p_sym .=> [En, κ_Rn, κ_Ln, Δn]), initial_values(eqs_a, u0))
n_ref_Δ = zeros(lΔ)
n_trans_Δ = zeros(lΔ)

for it = 1:lΔ
    Δn_ = Δn_ls[it]
    prob_ = ODEProblem(sys_a, dict_p_Δ(Δn_), (0.0, T[end]))
    sol_ = solve(prob_, Tsit5(); saveat = T)

    n_ref_Δ[it] = abs2(√(κ_Ln) * get_solution(sol_, a, eqs_a)(T[end]) + En)
    n_trans_Δ[it] = abs2(√(κ_Ln) * get_solution(sol_, a, eqs_a)(T[end]))
end
nothing # hide

#

p = plot(Δn_ls, n_ref_Δ; label = "reflection")
plot!(p, Δn_ls, n_trans_Δ; label = "transmission", ls = :dash)
plot!(p; xlabel = "Δ", grid = true, legend = :best, size = (600, 350))
p

# ## Two-sided cavity with atoms

# In the following, we include $N=2$ two-level atoms in the cavity and simulate the transmission and reflection of a coherent Gaussian pulse with a mean photon number of $|\alpha|^2 = 1/10$. We assume that the atoms are on resonance with the cavity, i.e. $\Delta = \Delta_c = \Delta_a$.

# Define the ModelingToolkit independent variable and register the classical drive as a
# function of it.
@independent_variables t
@register_symbolic Et(tt)

G_d_t = SLH(1, Et(t), 0)
H_ac = -Δ*(a'a + ∑σ(2, 2)) + g*(a'∑σ(1, 2) + a*∑σ(2, 1))
G_ac = SLH(1, √κ_L*a, H_ac)
G_ac_drive = (G_d_t ▷ G_ac) ⊞ SLH(1, √κ_R*a, 0)
nothing # hide

#

H2 = G_ac_drive.hamiltonian

#

J2_L = G_ac_drive.jump_operator[1]

#

J2_R = G_ac_drive.jump_operator[2]

# We derive the equations of motion for system with a second-order mean-field approximation. 

J_add = [√(γ)*σ(α, 1, 2) for α = 1:Natoms]
eqs2 = meanfield([a'a, σ(1, 2, 2)], H2, [J2_L, J2_R, J_add...]; order = 2)
nothing # hide

# ```math
# \begin{align} \frac{d}{dt} \langle a^\dagger a\rangle &= -1 i g \langle a^\dagger {\sigma1}^{{12}}\rangle + \langle a^\dagger a\rangle \left( -1.0 \left( \sqrt{\kappa{L}} \right)^{2} -1.0 \left( \sqrt{\kappa{R}} \right)^{2} \right) -1 i g \langle a^\dagger {\sigma2}^{{12}}\rangle + 1 i \langle a {\sigma2}^{{21}}\rangle g -1.0 \sqrt{\kappa{L}} \langle a\rangle \mathrm{Et}\left( t \right) -1.0 \sqrt{\kappa{L}} \langle a^\dagger\rangle \mathrm{Et}\left( t \right) + 1 i g \langle a {\sigma1}^{{21}}\rangle \\ \frac{d}{dt} \langle {\sigma1}^{{22}}\rangle &= 1 i g \langle a^\dagger {\sigma1}^{{12}}\rangle -1.0 \left( \sqrt{\gamma} \right)^{2} \langle {\sigma1}^{{22}}\rangle -1 i g \langle a {\sigma1}^{{21}}\rangle \end{align}
# ```

#

eqs2_c = complete(eqs2)
length(eqs2_c)

# Again, we defined the numerical parameters and the initial state of the system, create the ODE problem and solve the dynamics.

## numerical parameter
κ_Rn2 = κ_Ln2 = 2π*1
gn2 = 2π*0.4
Δn2 = 0.0
γn = 2π*0.1

p_sym2 = [κ_R, κ_L, Δ, g, γ]
p_num2 = [κ_Rn2, κ_Ln2, Δn2, gn2, γn]

σp = 10/κ_Ln2 # pulse width 
Tp = 4σp # pulse peak 
Tend = 3Tp
α0 = √(0.1) # √ of total photon number
Ω0 = α0*2*√(κ_Ln2)/(π^(1/4)*√(σp))
Ω1(t) = Ω0/2*exp(-(t-Tp)^2 / (2*σp^2))
Et(t) = Ω1(t)/√(κ_Ln2)

T2 = [0:0.001:1;]*Tend
ΔT = T2[2] - T2[1]
n_pulse = round(sum(abs2.(Et.(T2)))*ΔT, digits = 7)
@show n_pulse

dict_p2 = Dict(p_sym2 .=> p_num2)
nothing # hide

#

## initial state
bc1 = FockBasis(4)
ba = NLevelBasis(2)
b = bc1 ⊗ tensor([ba for i = 1:Natoms]...)
ψ0 = LazyKet(b, (fockstate(bc1, 0), [nlevelstate(ba, 1) for i = 1:Natoms]...))
u0_2 = initial_values(eqs2_c, initial_values(eqs2_c, ψ0)) # state -> vector -> u0 dict
nothing # hide

#

sys2 = mtkcompile(System(eqs2_c; name = :sys2)) # initial state
dict2 = merge(Dict(p_sym2 .=> p_num2), u0_2)
prob2 = ODEProblem(sys2, dict2, (0.0, T2[end]))
nothing # hide

#

sol2 = solve(prob2, Tsit5(); saveat = T2)
nothing # hide

#

n_ref2 = abs2.(√(κ_Ln2) .* get_solution(sol2, a, eqs2_c).(T2) .+ Et.(T2))
n_trans2 = abs2.(√(κ_Rn2) .* get_solution(sol2, a, eqs2_c).(T2))
nothing # hide

p = plot(
    T2,
    n_trans2;
    label = "transmission = $(round(sum(n_trans2) * ΔT / n_pulse * 100))%",
)
plot!(
    p,
    T2,
    n_ref2;
    ls = :dash,
    label = "reflection = $(round(sum(n_ref2) * ΔT / n_pulse * 100))%",
)
plot!(p; xlabel = "t", legend = :best, grid = true, size = (600, 350))
p

# We can see that all results agree with the full quantum dynamics of the example `Two-sided Cavity with Atoms`, which means the second order cumulant expansion is good approximation here. 

# ## Package versions

# These results were obtained using the following versions:

using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    [
        "QuantumInputOutput",
        "SecondQuantizedAlgebra",
        "QuantumCumulants",
        "ModelingToolkitBase",
        "OrdinaryDiffEq",
        "QuantumOpticsBase",
        "Plots",
    ],
    mode = PKGMODE_MANIFEST,
)
