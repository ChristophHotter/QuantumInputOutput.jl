# # Coherent-Feedback Squeezing with a Beam-Splitter Loop
#
# This example implements the feedback loop to enhance field squeezing using coherent feedback described in  
# [J. Gough and S. Wildfeuer, PRA 80, 042107 (2009)](http://dx.doi.org/10.1103/PhysRevA.80.042107), 
# see also Example VI.1 of [J. Combes, et al. Advances in Physics: X, 2:3, 784-888 (2017)](https://doi.org/10.1080/23746149.2017.1343097). 
# A degenerate parametric oscillator is concatenated with a beam splitter,
# and the internal wires are eliminated with the SLH feedback reduction rule.

using QuantumInputOutput
using SecondQuantizedAlgebra
using Symbolics: Symbolics
using SymbolicUtils
using Plots

#

## symbolic Hilbert space
hc = FockSpace(:c)

## symbolic operator
a = Destroy(hc, :a)

## symbolic parameters
κ = Symbolics.variable(Symbol("κ"); T = Real)
ϵ = Symbolics.variable(Symbol("ϵ"); T = Real)
η = Symbolics.variable(Symbol("η"); T = Real)
nothing # hide

# The open-loop OPO has one port, while the beam splitter has two ports. We first
# form the unconnected network and then apply the two feedback reductions. 

r = √(1 - η^2)
G_opo = SLH(1, √(κ) * a, 1im * ϵ * (a'^2 - a^2))
G_bs = SLH([-r η; η r], [0, 0], 0)
G_unconnected = G_opo ⊞ G_bs
G_loop = feedback(G_unconnected, 1 => 2, 2 => 1)
nothing # hide

#

S_loop = scattering(G_loop)
L_loop = lindblad(G_loop)[1]

#

H_loop = hamiltonian(G_loop)

# The feedback loop leaves the OPO Hamiltonian unchanged but rescales the
# coupling operator by $l = \eta / (1 + \sqrt{1-\eta^2})$.
# For a numerical illustration, we evaluate the effective damping factor $l^2$ as
# a function of the beam-splitter transmission coefficient.

η_grid = collect(0.0:0.005:0.999)
l2_grid = @. (η_grid / (1 + sqrt(1 - η_grid^2)))^2

p = plot(
    η_grid,
    l2_grid;
    lw = 2,
    label = "",
    xlabel = "η",
    ylabel = "effective damping fraction",
    grid = true,
    size = (520, 320),
)
hline!(p, [1.0]; color = :grey, ls = :dash, label = "open loop")
p

# ## Package versions

# These results were obtained using the following versions:

using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    ["QuantumInputOutput", "SecondQuantizedAlgebra", "Plots"],
    mode = PKGMODE_MANIFEST,
)
