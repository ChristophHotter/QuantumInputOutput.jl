```@meta
EditURL = "../../../examples/07-1_beamsplitter_loss__quantum-pulse.jl"
```

# Beam Splitter Loss

This example models loss of pulse in a Fock-state by mixing the pulse with a vacuum
port on a beam splitter. We analyze the output mode $v(t)$ which is the same as the
input mode $u(t)$.

````@example 07-1_beamsplitter_loss__quantum-pulse
using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using PyPlot
using LinearAlgebra
````

````@example 07-1_beamsplitter_loss__quantum-pulse
# symbolic Hilbert space and operators (virtual input and output modes)
hu = FockSpace(:u)
hv = FockSpace(:v)
h = hu ⊗ hv

au = Destroy(h, :a_u, 1)
av = Destroy(h, :a_v, 2)

# symbolic parameters
@rnumbers gu gv r t
nothing # hide
````

In our example, we only have one intput mode and one output mode, however, the beam splitter has two input and two output ports.
Since the number of input and output ports needs to match to cascade a system we need to create a padding element and add concatenate it to the corresponding SLH elements. The padding element `(1,0,0)` represents vacuum input and a non-tracked output, respectively.

````@example 07-1_beamsplitter_loss__quantum-pulse
# padding element for the unused port
G_p = SLH(1, 0, 0)

# beam splitter scattering matrix
S_bs = [r t; t -r]

# input cavity, beam splitter, and output cavity
G_u = SLH(1, gu * au, 0)
G_in = G_u ⊞ G_p
G_bs = SLH(S_bs, [0, 0], 0)
G_v = SLH(1, gv * av, 0)
````

G_out = G_v ⊞ G_p # reflection is tracked

````@example 07-1_beamsplitter_loss__quantum-pulse
G_out = G_p ⊞ G_v # transmission is tracked

G = G_in ▷ G_bs ▷ G_out
nothing # hide
````

````@example 07-1_beamsplitter_loss__quantum-pulse
H = get_hamiltonian(G)
````

````@example 07-1_beamsplitter_loss__quantum-pulse
L = get_lindblad(G)
````

````@example 07-1_beamsplitter_loss__quantum-pulse
# Gaussian input mode
γ_ = 1.0
σ = 1 / γ_
T_end = 12σ
u(t) = 1/(sqrt(σ)*π^(1/4)) * exp( -(t - 4σ)^2 / (2*σ^2) )

T = [0:0.004:1;] * T_end
ΔT = T[2] - T[1]

# time-dependent coupling for the virtual cavities
gu_t = u_to_gu(u, T)
gv_t = v_to_gv(u, T) # v(t) = u(t)
dict_p_t = Dict(gu => gu_t, gv => gv_t)

# beam splitter parameters
η = 0.2 # loss
r_ = sqrt(η) # reflection
t_ = sqrt(1 - η) # transmission
dict_p = Dict([t,r] .=> [t_, r_])
nothing # hide
````

As usual, we translate the symbolic system into numeric expressions and solve the dynamics with `QuantumOptics.jl`.

````@example 07-1_beamsplitter_loss__quantum-pulse
# numeric basis
n_ph = 4
bu = FockBasis(n_ph)
bv = FockBasis(n_ph)
b = bu ⊗ bv
au_qo = destroy(bu) ⊗ one(bv)
av_qo = one(bu) ⊗ destroy(bv)

# translate to numeric operators
H_QO = translate(H, b; parameter=dict_p, time_parameter=dict_p_t)
L_QO = [translate(Li, b; parameter=dict_p, time_parameter=dict_p_t) for Li in L]

function input_output(t, ρ)
    Ht = H_QO(t)
    J = [L_QO[1](t), L_QO[2](t)]
    return Ht, J, dagger.(J)
end
nothing # hide
````

````@example 07-1_beamsplitter_loss__quantum-pulse
# time evolution
ψ0 = fockstate(bu, n_ph) ⊗ fockstate(bv, 0)
time, ρt = timeevolution.master_dynamic(T, ψ0, input_output)
nothing # hide
````

````@example 07-1_beamsplitter_loss__quantum-pulse
n_u_t = real(expect(au_qo'au_qo, ρt))
n_v_t = real(expect(av_qo'av_qo, ρt))

ρv_end = ptrace(ρt[end], 1)
pop_n_ls = [ρv_end.data[i,i] for i=1:n_ph+1]
nothing # hide
````

We plot the mean photon number and the distribution of the Fock state components after the beam splitter interaction.
We can see that the mean photon number is reduced by $\eta = 20%$.

````@example 07-1_beamsplitter_loss__quantum-pulse
close("beam splitter loss") # hide
figure("beam splitter loss", figsize=(5.4, 4.2))
subplot(2,1,1)
plot(T, n_u_t .+ n_v_t)
xlabel("time")
ylabel("photon number")
grid(true)

subplot(2,1,2)
bar([0:n_ph;], pop_n_ls)
xlabel("Fock state component n")
ylabel("population")
tight_layout()
gcf()
````

Note that the calculation can also be performed in the interaction picture, which would be numerically beneficial and the loss of the input mode $u(t)$ can be directly observed.

## Package versions

These results were obtained using the following versions:

````@example 07-1_beamsplitter_loss__quantum-pulse
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

