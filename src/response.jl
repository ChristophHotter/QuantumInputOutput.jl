################################################################
### linear & nonlinear response: steady state, S-parameters, ###
### susceptibility and output spectra                        ###
################################################################
#
# These tools compute the *driven response* of an SLH network around a
# (pumped) steady state, the way a network/spectrum analyser measures it:
# a weak probe tone is sent into a port and the coherent field leaving
# another port is read out. The probe response is obtained exactly from the
# Liouvillian, without linearising the system Hamiltonian, so it stays valid
# for genuinely nonlinear systems such as a Kerr parametric oscillator. Only
# the probe is treated to linear order (the experimental weak-tone regime).
#
# All three observables reduce to the same primitive: the resolvent of the
# Liouvillian, `(ℒ + μ)⁻¹`, applied to a traceless operator. They share one
# Hessenberg factorisation (`LiouvillianResolvent`) so a frequency sweep costs
# one cheap shifted solve per point.

using QuantumOpticsBase: AbstractOperator, Operator, identityoperator, expect, liouvillian
using LinearAlgebra: hessenberg, tr, I

# `dagger` is imported from QuantumOpticsBase by the top-level module.

# ──────────────────────────────────────────────
# A. Output field
# ──────────────────────────────────────────────

"""
    output_field(G::SLH, k; input=0)

Return the **output-field operator** of port `k` of an SLH network, i.e. the
operator whose expectation value is the travelling field leaving that port.

Following the input–output boundary relation ``b_{\\mathrm{out},k} =
(S\\,b_{\\mathrm{in}})_k - L_k``, the operator part is `-lindblad(G)[k]`. With no
field injected (`input=0`) this is just `-L_k`; passing a vector `input` of
coherent input amplitudes per port adds the directly-scattered classical part
``(S\\,\\alpha)_k``.

This is the canonical observable behind `scattering_parameter` and `power_spectrum`.
"""
function output_field(G::SLH, k::Integer; input = 0)
    L = lindblad(G)[k]
    if input == 0 || (input isa AbstractVector && all(iszero, input))
        return -L
    end
    α = input isa AbstractVector ? input : error("`input` must be a vector of port amplitudes")
    S = scattering(G)
    offset = sum(S[k, j] * α[j] for j in axes(S, 2))
    return offset - L
end

# ──────────────────────────────────────────────
# Liouvillian resolvent (shared engine)
# ──────────────────────────────────────────────

"""
    LiouvillianResolvent

Precomputed factorisation for evaluating ``(\\mathcal{L} + \\mu)^{-1} X`` at many
shifts `μ`, where ``\\mathcal{L}`` is the Liouvillian and `X` is a traceless
operator. The zero mode of ``\\mathcal{L}`` (the steady state) is deflated so the
resolvent is regular at every real frequency, including ``\\omega = 0``; this is
exact because all response right-hand sides used here are traceless.

Build with [`liouvillian_resolvent`](@ref); apply with
[`susceptibility`](@ref), [`scattering_parameter`](@ref) or [`power_spectrum`](@ref).
"""
struct LiouvillianResolvent{BT,OT,FT}
    basis::BT
    ρ_ss::OT
    F::FT
    d::Int
end

"""
    liouvillian_resolvent(H, J, ρ_ss; deflation=1.0)

Construct a [`LiouvillianResolvent`](@ref) from a numeric Hamiltonian `H`,
collapse operators `J` and a precomputed steady state `ρ_ss`. Builds a Hessenberg
factorisation of the deflated Liouvillian for fast shifted solves; the zero mode
is deflated using `ρ_ss` as the right null vector.

Compute `ρ_ss` yourself with whichever solver suits your system, e.g.
`QuantumOptics.steadystate.eigenvector(H, collect(J))`, and reuse the resulting
resolvent across [`susceptibility`](@ref), [`scattering_parameter`](@ref) and
[`power_spectrum`](@ref).
"""
function liouvillian_resolvent(
    H::AbstractOperator,
    J,
    ρ_ss::AbstractOperator;
    deflation = 1.0,
)
    b = H.basis_l
    d = length(b)
    Jvec = collect(J)                                        # liouvillian needs a plain Vector
    ℒ = liouvillian(H, Jvec)
    ρss = ρ_ss
    vecI = reshape(Matrix(identityoperator(b).data), d * d)   # ⟨𝟙| : the trace functional
    vecρ = reshape(Matrix(ρss.data), d * d)                   # |ρ_ss⟩ : the right null vector
    F = hessenberg(Matrix(ℒ.data) + deflation * (vecρ * vecI'))
    return LiouvillianResolvent(b, ρss, F, d)
end

# Solve (ℒ + μ) X = rhs for the operator X. `rhs` must be traceless.
function _resolvent_solve(R::LiouvillianResolvent, μ::Number, rhs::AbstractOperator)
    x = (R.F + μ * I) \ reshape(Matrix(rhs.data), R.d * R.d)
    return Operator(R.basis, R.basis, reshape(x, R.d, R.d))
end

_commutator(A, B) = A * B - B * A

# ──────────────────────────────────────────────
# B. Susceptibility (Kubo) and S-parameters
# ──────────────────────────────────────────────

"""
    susceptibility(R::LiouvillianResolvent, A, B, ω)
    susceptibility(G::SLH, b, ρ_ss, A, B, ω; parameter=Dict(), operators=Dict(), kwargs...)

Kubo linear-response susceptibility

```math
\\chi_{AB}(\\omega) = i\\,\\mathrm{Tr}\\!\\big[A\\,(\\mathcal{L} + i\\omega)^{-1}\\,[B, \\rho_{ss}]\\big],
```

the response of ``\\langle A \\rangle`` at frequency `ω` to a weak perturbation
``H_1(t) = f(t)\\,B + \\text{h.c.}`` with ``f(t) = e^{-i\\omega t}``. The system is
**not** linearised: ``\\mathcal{L}`` is the full (nonlinear) Liouvillian and only
the perturbation is taken to first order. `ω` may be a scalar or a vector.

Pass numeric operators `A`, `B` and a prebuilt resolvent `R`, or symbolic `A`,
`B` together with an `SLH` network `G`, basis `b` and a precomputed steady state
`ρ_ss` (from `QuantumOptics.steadystate` or any solver) to translate and build the
resolvent in one call. This is the general engine; for the standard transmission
or reflection coefficient use [`scattering_parameter`](@ref).
"""
function susceptibility(R::LiouvillianResolvent, A::AbstractOperator, B::AbstractOperator, ω::Real)
    X = _resolvent_solve(R, im * ω, _commutator(B, R.ρ_ss))
    return im * tr(A.data * X.data)
end

susceptibility(R::LiouvillianResolvent, A, B, ω::AbstractVector) =
    [susceptibility(R, A, B, ωk) for ωk in ω]

function susceptibility(
    G::SLH,
    b,
    ρ_ss::AbstractOperator,
    A,
    B,
    ω;
    parameter = Dict(),
    operators = Dict(),
    kwargs...,
)
    H, J = translate_qo(G, b; parameter, operators)
    A_ = translate_qo(A, b; parameter, operators)
    B_ = translate_qo(B, b; parameter, operators)
    R = liouvillian_resolvent(H, J, ρ_ss; kwargs...)
    return susceptibility(R, A_, B_, ω)
end

"""
    scattering_parameter(G::SLH, b, ρ_ss; in_port=1, out_port=2, omega, parameter=Dict(), operators=Dict(), kwargs...)

Scattering parameter ``S_{\\text{out\\_port},\\text{in\\_port}}(\\omega)`` of an SLH
network: the coherent field leaving `out_port` per unit weak probe sent into
`in_port`, as a function of probe detuning `omega` (in the rotating frame of the
pump). With `out_port ≠ in_port` this is a transmission coefficient (``S_{21}``);
with `out_port == in_port`, a reflection coefficient (``S_{11}``).

`ρ_ss` is the precomputed pumped steady state (from `QuantumOptics.steadystate` or any
solver) at the same working point; compute it once and reuse it across
`scattering_parameter`/`power_spectrum`/`susceptibility`.

```math
S_{\\text{out},\\text{in}}(\\omega) = S_{\\text{out},\\text{in}}
    + \\mathrm{Tr}\\!\\big[L_{\\text{out}}\\,(\\mathcal{L} + i\\omega)^{-1}\\,[L_{\\text{in}}^\\dagger, \\rho_{ss}]\\big],
```

where ``L_k`` is the port-`k` coupling operator (`lindblad(G)[k]`) and the first
term is the direct scattering matrix element `scattering(G)[out_port, in_port]`.
The Liouvillian is built from **all** ports of `G`, so internal loss ports
correctly reduce the transmission. Because the probe enters only to linear order,
the result is the exact weak-tone response of the full nonlinear network and
satisfies ``|S_{11}|^2 + |S_{21}|^2 = 1`` for a lossless two-port.

`kwargs` (e.g. `method`, `rho_ss`, `deflation`) are forwarded to
[`liouvillian_resolvent`](@ref).
"""
function scattering_parameter(
    G::SLH,
    b,
    ρ_ss::AbstractOperator;
    in_port::Integer = 1,
    out_port::Integer = 2,
    omega,
    parameter = Dict(),
    operators = Dict(),
    kwargs...,
)
    H, J = translate_qo(G, b; parameter, operators)
    L_in = J[in_port]
    L_out = J[out_port]
    R = liouvillian_resolvent(H, J, ρ_ss; kwargs...)
    direct = scattering(G)[out_port, in_port]
    rhs = _commutator(dagger(L_in), R.ρ_ss)             # [L_in†, ρ_ss], traceless
    return [direct + tr(L_out.data * _resolvent_solve(R, im * ωk, rhs).data) for ωk in omega]
end

# ──────────────────────────────────────────────
# C. Output power / squeezing spectrum
# ──────────────────────────────────────────────

"""
    power_spectrum(G::SLH, b, ρ_ss; port=1, omega, quadrature=nothing,
                   parameter=Dict(), operators=Dict(), subtract_mean=true, kwargs...)

Steady-state output spectrum at `port` of an SLH network around the precomputed
steady state `ρ_ss` X. Two physically
distinct spectra are available.

**Emission spectrum** (default, `quadrature=nothing`). With ``b_{\\mathrm{out}} =
-L`` the port output field and ``\\bar b_{\\mathrm{out}} = b_{\\mathrm{out}} -
\\langle b_{\\mathrm{out}}\\rangle``,

```math
S(\\omega) = -2\\,\\mathrm{Re}\\,\\mathrm{Tr}\\!\\big[\\bar b_{\\mathrm{out}}^\\dagger\\,(\\mathcal{L} - i\\omega)^{-1}\\,(\\bar b_{\\mathrm{out}}\\,\\rho_{ss})\\big],
```

the Fourier transform of ``\\langle \\bar b_{\\mathrm{out}}^\\dagger(\\tau)\\,\\bar
b_{\\mathrm{out}}(0)\\rangle`` (convention ``\\int e^{-i\\omega\\tau}\\,d\\tau``,
matching `QuantumOptics.timecorrelations.spectrum`). This is the
normally-ordered parametric-fluorescence / emission spectrum; its vacuum level
is zero.

**Squeezing spectrum** (`quadrature=θ`). The vacuum-normalised homodyne spectrum
of the output quadrature ``X_\\theta = (e^{-i\\theta} b_{\\mathrm{out}} +
e^{i\\theta} b_{\\mathrm{out}}^\\dagger)/\\sqrt{2}``,

```math
S_\\theta(\\omega) = 1 + 2\\,\\mathrm{Re}\\Big[
    e^{-2i\\theta}\\big(C_{LL}(\\omega)+C_{LL}(-\\omega)\\big)
    + C_{L^\\dagger L}(\\omega)+C_{L^\\dagger L}(-\\omega)\\Big],
\\quad
C_{PQ}(\\omega) = -\\mathrm{Tr}\\!\\big[\\delta P\\,(\\mathcal{L}+i\\omega)^{-1}(\\delta Q\\,\\rho_{ss})\\big],
```

with ``L`` the port coupling operator and ``\\delta O = O - \\langle O\\rangle``.
The shot-noise floor is ``1``; a squeezed quadrature dips below it. For a
degenerate parametric amplifier this reproduces the Collett–Gardiner result
``S_\\mp(\\omega) = 1 \\mp 2\\kappa G/((\\kappa/2 \\pm G)^2 + \\omega^2)`` exactly.

Set `subtract_mean=false` to keep the coherent part. `kwargs` are forwarded to
[`liouvillian_resolvent`](@ref).
"""
function power_spectrum(
    G::SLH,
    b,
    ρ_ss::AbstractOperator;
    port::Integer = 1,
    omega,
    quadrature = nothing,
    parameter = Dict(),
    operators = Dict(),
    subtract_mean = true,
    kwargs...,
)
    H, J = translate_qo(G, b; parameter, operators)
    R = liouvillian_resolvent(H, J, ρ_ss; kwargs...)
    Lp = J[port]
    if quadrature === nothing
        Abar = subtract_mean ? Lp - expect(Lp, R.ρ_ss) * identityoperator(b) : Lp
        Adag = dagger(Abar)
        rhs = Abar * R.ρ_ss                                  # traceless (⟨Ā⟩ = 0)
        return [-2 * real(tr(Adag.data * _resolvent_solve(R, -im * ωk, rhs).data)) for ωk in omega]
    end
    θ = quadrature
    δL = subtract_mean ? Lp - expect(Lp, R.ρ_ss) * identityoperator(b) : Lp
    δLd = dagger(δL)
    # one-sided spectrum C_PQ(ω) = ∫₀^∞ e^{iωτ}⟨δP(τ)δQ(0)⟩dτ = -Tr[δP (ℒ+iω)⁻¹ (δQ ρ_ss)]
    C(P, Q, ω) = -tr(P.data * _resolvent_solve(R, im * ω, Q * R.ρ_ss).data)
    return [
        1 + 2 * real(
            exp(-2im * θ) * (C(δL, δL, ωk) + C(δL, δL, -ωk)) +
            (C(δLd, δL, ωk) + C(δLd, δL, -ωk)),
        ) for ωk in omega
    ]
end
