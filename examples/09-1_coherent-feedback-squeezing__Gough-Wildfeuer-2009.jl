# # Coherent-Feedback Squeezing with a Beam-Splitter Loop
#
# This example implements the feedback loop from
# Gough and Wildfeuer, "Enhancement of field squeezing using coherent feedback"
# (2009). A degenerate parametric oscillator is concatenated with a beam splitter,
# and the internal wires are eliminated with the SLH feedback reduction rule.

# TODO: article reference, plot wigner function of squeezed state(?), read carefully  

using QuantumInputOutput
using SecondQuantizedAlgebra
using SymbolicUtils
using Plots

#

## symbolic Hilbert space
hc = FockSpace(:c)

## symbolic operator
a = Destroy(hc, :a, 1)

## symbolic parameters
κ = rnumber("κ")
ϵ = rnumber("ϵ")
η = rnumber("η")
nothing # hide

# The open-loop OPO has one port, while the beam splitter has two ports. We first
# form the unconnected network and then apply the two feedback reductions
# described in Example VI.1 of the SLH review by Combes et al.

r = √(1 - η^2)
G_opo = SLH(1, √(κ) * a, 1im * ϵ * (a'^2 - a^2))
G_bs = SLH([-r η; η r], [0, 0], 0)
G_unconnected = G_opo ⊞ G_bs
G_loop = feedback(G_unconnected, 1 => 2, 2 => 1)
nothing # hide

#

S_loop = get_scattering(G_loop)
L_loop = get_lindblad(G_loop)[1]
H_loop = get_hamiltonian(G_loop)

# The feedback loop leaves the OPO Hamiltonian unchanged but rescales the
# coupling operator by
# $l = \eta / (1 + \sqrt{1-\eta^2})$.

l = η / (1 + r)
L_expected = l * √(κ) * a
nothing # hide

#

@show S_loop
@show simplify(L_loop - L_expected)
@show H_loop

# For a numerical illustration, we evaluate the effective damping factor $l^2$ as
# a function of the beam-splitter transmission coefficient.

η_grid = collect(0.0:0.005:0.999)
l2_grid = @. (η_grid / (1 + sqrt(1 - η_grid^2)))^2

p = plot(
    η_grid,
    l2_grid;
    lw = 2,
    label = "l²",
    xlabel = "η",
    ylabel = "effective damping fraction",
    title = "Feedback-reduced OPO coupling",
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
