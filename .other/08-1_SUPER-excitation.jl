# # SUPER Excitation with Two Off-Resonant Photon Modes

# TODO: rewrite
# In this example, we model SUPER (swing-up of a quantum emitter) using two off-resonant photon modes that jointly invert a two-level emitter. This is the few-photon, fully quantum analog of the detuned-pulse swing-up scheme introduced in [T. K. Bracht et al., PRX Quantum 2, 040354 (2021)](https://doi.org/10.1103/PRXQuantum.2.040354) and extended to two quantized modes in [Q. W. Richter et al., Phys. Rev. Research 7, 013079 (2025)](https://doi.org/10.1103/PhysRevResearch.7.013079). We use a cascaded input-output model with two input virtual cavities, a two-level scatterer, and two output cavities, and we evolve the equations of motion with QuantumCumulants.jl.

using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumCumulants
using ModelingToolkit
using QuantumOptics
using OrdinaryDiffEq
using PyPlot
using SpecialFunctions: loggamma # TODO: move to QOB.jl

#

## symbolic Hilbert space
hu1 = FockSpace(:u1)
hu2 = FockSpace(:u2)
hs1 = NLevelSpace(:atom, 2)
hv1 = FockSpace(:v1)
hv2 = FockSpace(:v2)
h = hu1 ⊗ hu2 ⊗ hs1 ⊗ hv1 ⊗ hv2

## symbolic operators
au1 = Destroy(h, :au_1, 1)
au2 = Destroy(h, :au_2, 2)
σ(i,j) = Transition(h, :σ, i, j, 3)
av1 = Destroy(h, :av_1, 4)
av2 = Destroy(h, :av_2, 5)

## symbolic parameters
@rnumbers γ Δ
gu1, gu2, gv1, gv2 = cnumbers("gu_1 gu_2 gv_1 gv_2")
@syms t::Real # time parameter

#

# We build the cascaded SLH model and extract the Hamiltonian and the Lindblad operator.

G_u2 = SLH(1, gu2*au2, 0)            # input cavity 2
G_u1 = SLH(1, gu1*au1, 0)            # input cavity 1
G_tl = SLH(1, √(γ)*σ(1,2), Δ*σ(2,2)) # two-level scatterer
G_v1 = SLH(1, gv1*av1, 0)            # output cavity 1
G_v2 = SLH(1, gv2*av2, 0)            # output cavity 2

G_cas = cascade(G_u2, G_u1, G_tl, G_v1, G_v2)
H = get_hamiltonian(G_cas)
L = get_lindblad(G_cas)[1]
Ld = adjoint(L)

#

# We insert time-dependent couplings for the input/output virtual cavities.
# TODO: explain @register_symbolic

@register_symbolic gu1_t(t)
@register_symbolic gu2_t(t)
@register_symbolic gv1_t(t)
@register_symbolic gv2_t(t)
@register_symbolic gu1_c_t(t)
@register_symbolic gu2_c_t(t)
@register_symbolic gv1_c_t(t)
@register_symbolic gv2_c_t(t)

g_syms = [gu1, gu2, gv1, gv2]
g_time = [gu1_t(t), gu2_t(t), gv1_t(t), gv2_t(t)]
g_time_c = [gu1_c_t(t), gu2_c_t(t), gv1_c_t(t), gv2_c_t(t)]
dict_p_t = Dict([g_syms; conj.(g_syms)] .=> [g_time; g_time_c])

H_t = substitute(H, dict_p_t)
L_t = substitute(L, dict_p_t)
Ld_t = substitute(Ld, dict_p_t)

#

# Two detuned Gaussian input pulses define the SUPER driving modes.

dt = 1e-4
T0 = dt
Tend = 24.0
T = [T0:dt:Tend;]

γ_ = 1e-2
Δ_ = 0.0
Δ1_ = -2π*1.934
Δ2_ = -2π*4.634
α1_ = 22.65π
α2_ = 19.29π
σ1_ = 2.4
σ2_ = 3.04
τ1_ = 12.0
τ2_ = 11.27

u1(t_) = 1/(√(σ1_)*π^(1/4)) * exp( -(t_ - τ1_)^2 / (2*σ1_^2) ) * exp(-1im*Δ1_*t_)
u2(t_) = 1/(√(σ2_)*π^(1/4)) * exp( -(t_ - τ2_)^2 / (2*σ2_^2) ) * exp(-1im*Δ2_*t_)

gu1_t_ = u_to_gu_Gauss(τ1_, σ1_)
gu1_t(t) = gu1_t_(t) # TODO: needed?
gv1_t_ = v_to_gv_Gauss(τ1_, σ1_)
gv1_t(t) = gv1_t_(t)
gu1_c_t(t) = conj(gu1_t(t))
gv1_c_t(t) = conj(gv1_t(t))

u_fcts = [u1, u2]
gu2_ = u_to_gu(ui_to_u_i_im1(u_fcts, T, 2), T)
gu2_t(t) = gu2_(t)
gu2_c_t(t) = conj(gu2_t(t))

v_fcts = [u1, u2]
gv2_ = v_to_gv(vi_to_v_i_im1(v_fcts, T, 2), T)
gv2_t(t) = gv2_(t)
gv2_c_t(t) = conj(gv2_t(t))

#

# Mean-field equations (first order) and time evolution with QuantumCumulants.jl.

order = 1
ops = [au1, au2, σ(2,2), σ(2,1), av1, av2]
eqs = meanfield(ops, H_t, [L_t]; Jdagger=[Ld_t], order=order)
complete!(eqs)

#

# Initial state: coherent input modes, atom in the ground state, vacuum outputs.

α1_io = α1_ / (2*√(2)*π^(1/4)*√(σ1_*γ_))
α2_io = α2_ / (2*√(2)*π^(1/4)*√(σ2_*γ_))

poisson(k, λ) = exp(k*log(λ) - λ - loggamma(k+1))
function QuantumOpticsBase.coherentstate!(ket::Ket, b::FockBasis, alpha::Number)
    C = eltype(ket)
    alpha = C(alpha)
    data = ket.data
    λ = abs2(alpha)
    ϕ = angle(alpha)
    offset = b.offset
    @inbounds for n=offset:b.N
        i = n - offset + 1
        data[i] = sqrt(poisson(n, λ)) * exp(1im*ϕ*n)
    end
    return ket
end

rd(x) = round(Int, x)
function ψ0_u1u2s1v1v2(α1, α2, eqs)
    bu1 = FockBasis(rd(α1^2 + 20α1), rd(α1^2 - 20α1))
    bu2 = FockBasis(rd(α2^2 + 20α2), rd(α2^2 - 20α2))
    bs1 = NLevelBasis(2)
    bv1 = FockBasis(1)
    bv2 = FockBasis(1)
    b = tensor([bu1, bu2, bs1, bv1, bv2]...)

    ψu1 = coherentstate(bu1, α1)
    ψu2 = coherentstate(bu2, α2)
    ψs1 = nlevelstate(bs1, 1)
    ψv1 = fockstate(bv1, 0)
    ψv2 = fockstate(bv2, 0)
    ψ0 = LazyKet(b, (ψu1, ψu2, ψs1, ψv1, ψv2))
    return initial_values(eqs, ψ0)
end

u0 = ψ0_u1u2s1v1v2(α1_io, α2_io, eqs)

@named sys = ODESystem(eqs)
ps = [γ, Δ]
p0 = [γ_, Δ_]
dict_u0 = merge(Dict(unknowns(sys) .=> u0), Dict(ps .=> p0))
prob = ODEProblem(sys, dict_u0, (T0, Tend))
nothing # hide

#

abstol = 1e-10; reltol = 1e-10
sol = solve(prob, Tsit5(); abstol, reltol)
nothing # hide

#

# Plot population and photon number transfer between input and output modes.

t_ = sol.t
s22 = real.(sol[σ(2,2)])
nu1 = abs2.(sol[au1])
nu2 = abs2.(sol[au2])
nv1 = abs2.(sol[av1])
nv2 = abs2.(sol[av2])

close("population") # hide
figure("population")
plot(t_, s22, label=L"\langle \sigma^{22} \rangle")
xlabel("time")
ylabel("population")
grid(true)
legend()
tight_layout()
gcf()

#

# pygui(true)
close("modes") # hide
figure("modes", figsize=(5,8))
subplot(311)
plot(t_, nu1, ls="--", label=L"\langle a^\dagger a \rangle_{u_1}")
plot(t_, nv1, ls="-", label=L"\langle a^\dagger a \rangle_{v_1}")
grid(true)
ylabel("photons")
legend(loc="center right")

subplot(312)
plot(t_, nu2, ls="--", label=L"\langle a^\dagger a \rangle_{u_2}")
plot(t_, nv2, ls="-", label=L"\langle a^\dagger a \rangle_{v_2}")
grid(true)
xlabel("time")
ylabel("photons")
legend(loc="center right")

subplot(313)
plot(t_, nu1 .+ nv1 .- nu1[1], label="mode 1")
plot(t_, nu2 .+ nv2 .- nu2[1], label="mode 2")
grid(true)
xlabel("time")
ylabel("photons")
legend(loc="lower right")
tight_layout()
display(gcf())

# TODO: final time of the pulse??

#

# ## Package versions

using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    ["QuantumInputOutput", "SecondQuantizedAlgebra", "QuantumCumulants", "OrdinaryDiffEq", "PyPlot"],
    mode = PKGMODE_MANIFEST,
)
