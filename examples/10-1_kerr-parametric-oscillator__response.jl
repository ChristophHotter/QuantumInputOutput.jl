# # Driven response of a Kerr parametric oscillator

# A Kerr parametric oscillator (KPO) is a single bosonic mode with a Kerr nonlinearity
# that is pumped by a two-photon (parametric) drive. Written in the frame rotating at half
# the pump frequency, its Hamiltonian is time independent,
#
# ```math
# H = \Delta\, a^\dagger a + K\, a^\dagger a^\dagger a a + \frac{G}{2}\left(a^\dagger a^\dagger + a a\right),
# ```
#
# with detuning ``\Delta``, Kerr coefficient ``K`` and parametric pump amplitude ``G``. We
# couple the mode to a transmission line through two ports with rates ``\kappa_1`` (input)
# and ``\kappa_2`` (output), so the device is a **two-port** scatterer.
#
# This example shows how to measure such a device the way an experimentalist would, with a
# weak probe tone, using three tools from `QuantumInputOutput`:
#
# * [`scattering_parameter`](@ref) — the transmission ``S_{21}(\omega)`` and reflection ``S_{11}(\omega)``,
# * [`power_spectrum`](@ref) — the output emission and squeezing spectra,
# * [`susceptibility`](@ref) — the underlying Kubo linear-response engine.
#
# The probe is treated only to linear order (the weak-tone regime), but the system itself is
# never linearised: the full nonlinear master equation enters through the Liouvillian. No
# steady-state Jacobian or fluctuation expansion is needed.

using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using Plots
using LaTeXStrings

# ## Building the two-port KPO as an SLH component

# The mode lives in a single Fock space. A two-port cavity has a two-component Lindblad
# vector ``L = (\sqrt{\kappa_1}\,a,\ \sqrt{\kappa_2}\,a)`` and an identity scattering matrix.

h = FockSpace(:kpo)
a = Destroy(h, :a)

@variables Δ::Real K::Real G::Real κ1::Real κ2::Real

H = Δ*a'*a + K*a'*a'*a*a + (G/2)*(a'*a' + a*a)
Gkpo = SLH([1 0; 0 1], [√(κ1)*a, √(κ2)*a], H)
nothing # hide

# The numeric basis and a working point below the parametric oscillation threshold
# ``G_\mathrm{th} = (\kappa_1+\kappa_2)/2``.

b = FockBasis(20)
κ1_, κ2_ = 0.5, 0.5
params(G_; Δ_ = 0.0, K_ = 0.1) = Dict(Δ => Δ_, K => K_, G => G_, κ1 => κ1_, κ2 => κ2_)
nothing # hide

# ## (A) Steady state and output field

# The pumped steady state is the working point we linearise the probe around, and it is the
# one expensive, method-sensitive ingredient shared by all three tools below: you compute it
# once with whichever solver suits your system and pass it in. Here we translate the network to
# numeric operators with [`translate_qo`](@ref) and take the exact null vector of the
# Liouvillian with `QuantumOptics.steadystate`. The output field operator of a port follows the
# input–output relation ``b_{\mathrm{out},k} = (S\,b_{\mathrm{in}})_k - L_k``.

function steady(p)
    Hn, Jn = translate_qo(Gkpo, b; parameter = p)
    return steadystate.eigenvector(Hn, collect(Jn))
end
ρ_passive = steady(params(0.0))
ρ_mid     = steady(params(0.3))
ρ_near    = steady(params(0.45))
b_out = output_field(Gkpo, 2)             # symbolic output operator of port 2
nothing # hide

# ## (B) Transmission ``S_{21}(\omega)`` and parametric gain

# We sweep the probe detuning ``\omega``. With the pump off (``G=0``) the device is a passive
# two-port cavity and ``|S_{21}| \le 1``. As ``G`` approaches threshold the parametric drive
# adds gain and ``|S_{21}|`` exceeds one: the KPO acts as a phase-insensitive amplifier.

ω = collect(range(-1.5, 1.5; length = 401))
S21_passive = scattering_parameter(Gkpo, b, ρ_passive; in_port = 1, out_port = 2, omega = ω, parameter = params(0.0))
S21_mid     = scattering_parameter(Gkpo, b, ρ_mid;     in_port = 1, out_port = 2, omega = ω, parameter = params(0.3))
S21_near    = scattering_parameter(Gkpo, b, ρ_near;    in_port = 1, out_port = 2, omega = ω, parameter = params(0.45))

p1 = plot(ω, abs.(S21_passive); label = "G = 0 (passive)", color = :black)
plot!(p1, ω, abs.(S21_mid);  label = "G = 0.3", color = :blue)
plot!(p1, ω, abs.(S21_near); label = "G = 0.45 (near threshold 0.5)", color = :red)
hline!(p1, [1.0]; ls = :dot, color = :gray, label = "")
plot!(p1; xlabel = L"probe detuning $\omega$", ylabel = L"|S_{21}(\omega)|",
      title = "KPO transmission: parametric gain", legend = :topright, size = (650, 350))
p1

# For a lossless passive two-port, transmission and reflection conserve energy,
# ``|S_{11}|^2 + |S_{21}|^2 = 1``. This is a useful sanity check on the working point.

S11_passive = scattering_parameter(Gkpo, b, ρ_passive; in_port = 1, out_port = 1, omega = ω, parameter = params(0.0))
@info "max energy deviation (passive): " maximum(abs.(abs2.(S11_passive) .+ abs2.(S21_passive) .- 1))

# ## (C) Output emission and squeezing spectra

# With the pump on, the KPO emits parametric fluorescence even with no probe. The emission
# spectrum of the output port is the Fourier transform of ``\langle b_\mathrm{out}^\dagger(\tau)\,b_\mathrm{out}(0)\rangle``
# and has a zero vacuum floor. Passing a `quadrature` angle instead returns the vacuum-normalised
# homodyne (squeezing) spectrum, whose shot-noise floor is one: the squeezed quadrature drops
# below it while the orthogonal one is anti-squeezed above it.

S_emit = power_spectrum(Gkpo, b, ρ_mid; port = 2, omega = ω, parameter = params(0.3))
S_sqz  = power_spectrum(Gkpo, b, ρ_mid; port = 2, omega = ω, parameter = params(0.3), quadrature = π/4)
S_anti = power_spectrum(Gkpo, b, ρ_mid; port = 2, omega = ω, parameter = params(0.3), quadrature = 3π/4)

p2 = plot(ω, S_emit; label = "emission", color = :black)
plot!(p2, ω, S_sqz;  label = L"squeezed quadrature $\theta=\pi/4$", color = :blue)
plot!(p2, ω, S_anti; label = L"anti-squeezed $\theta=3\pi/4$", color = :red)
hline!(p2, [1.0]; ls = :dot, color = :gray, label = "vacuum / shot-noise floor")
plot!(p2; xlabel = L"frequency $\omega$", ylabel = "output spectrum",
      title = "KPO output spectra", legend = :topright, size = (650, 380))
p2

# ## (B, general) The susceptibility engine

# [`scattering_parameter`](@ref) is a thin wrapper over [`susceptibility`](@ref), the Kubo response
# ``\chi_{AB}(\omega) = i\,\mathrm{Tr}[A\,(\mathcal{L}+i\omega)^{-1}[B,\rho_{ss}]]``. The
# anomalous (idler) response of a parametric system, the response of ``a`` to a drive that
# couples through ``a`` rather than ``a^\dagger``, is nonzero and is what makes the gain
# phase sensitive. It is one call away:

χ_normal   = susceptibility(Gkpo, b, ρ_mid, a, a',  ω; parameter = params(0.3))   # ⟨a⟩ response to a† drive
χ_anomalous = susceptibility(Gkpo, b, ρ_mid, a, a,  ω; parameter = params(0.3))   # ⟨a⟩ response to a  drive (idler)

p3 = plot(ω, abs.(χ_normal);    label = L"\chi_{a,a^\dagger} (normal)", color = :blue)
plot!(p3, ω, abs.(χ_anomalous); label = L"\chi_{a,a} (anomalous/idler)", color = :red)
plot!(p3; xlabel = L"\omega", ylabel = L"|\chi(\omega)|",
      title = "Normal vs anomalous response", legend = :topright, size = (650, 350))
p3

# ## Package versions

using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(["QuantumInputOutput", "SecondQuantizedAlgebra", "QuantumOptics", "Plots"];
           mode = PKGMODE_MANIFEST)
