######################################
### functions for virtual cavities ###
######################################

const _tol_div = 1e-10
const _extrapolate = ExtrapolationType.Extension
const _ϵ = 1e-10

"""
    PulseCoupling(f)

A single-argument time function `g(t)` used as a `to_numeric` `time_parameter` value.
Wraps any evaluator `f` (a sampled interpolation, an analytic closure, ...) behind exactly
one call method, so it satisfies SecondQuantizedAlgebra's single-arity `time_parameter`
contract by construction (a raw interpolation is rejected for having methods of conflicting
arity). Callable as `g(t)` and broadcastable as `g.(T)`.
"""
struct PulseCoupling{F}
    f::F
end
(g::PulseCoupling)(t::Real) = g.f(t)

_interpolate_mode(mode::AbstractVector, T::AbstractVector) =
    LinearInterpolation(mode, T; extrapolation = _extrapolate)
_interpolate_mode(mode, T::AbstractVector) = mode

# ──────────────────────────────────────────────
# Gaussian pulse type
# ──────────────────────────────────────────────

"""
    Gaussian(τ, σ; δ=0)

Gaussian pulse shape descriptor with center time `τ`, width `σ`, and detuning `δ`.
Use with `coupling_input` and `coupling_output` for analytical coupling formulas.
"""
struct Gaussian{T}
    τ::T
    σ::T
    δ::T
end
Gaussian(τ, σ; δ = zero(τ)) = Gaussian(promote(τ, σ, δ)...)

function _gaussian_mode(g::Gaussian)
    τ, σ, δ = g.τ, g.σ, g.δ
    mode = if δ == 0
        t -> 1 / (√(σ) * π^(1 / 4)) * exp(-0.5 * (t - τ)^2 / σ^2)
    else
        t -> 1 / (√(σ) * π^(1 / 4)) * exp(-0.5 * (t - τ)^2 / σ^2) * exp(1im * δ * t)
    end
    ∫mode2 = t -> 0.5 * (erf((t - τ) / σ) + erf(τ / σ))
    return mode, ∫mode2
end

# ──────────────────────────────────────────────
# Shared coupling core
# ──────────────────────────────────────────────

function _coupling_from_mode(mode::Vector, T::Vector, denominator)
    nt = length(T)
    mode_sq = Vector{Float64}(undef, nt)
    map!(abs2, mode_sq, mode)
    ∫mode2 = cumul_integrate(T, mode_sq)
    coupling = zeros(ComplexF64, nt)
    @inbounds for i = 1:nt
        denom = denominator(∫mode2[i])
        if sqrt(abs(denom)) > _tol_div
            coupling[i] = mode[i]' / sqrt(denom)
        end
    end
    return PulseCoupling(LinearInterpolation(coupling, T; extrapolation = _extrapolate))
end

# ──────────────────────────────────────────────
# coupling_input / coupling_output
# ──────────────────────────────────────────────

"""
    coupling_input(u, T)

Compute the virtual-cavity input coupling ``g_u(t)`` from an input mode `u(t)`
sampled on time grid `T`. Returns a [`PulseCoupling`](@ref) (callable as `g(t)`,
broadcastable as `g.(T)`) that plugs directly into a `to_numeric` `time_parameter`.
"""
coupling_input(u::Vector, T::Vector) = _coupling_from_mode(u, T, x -> abs(1 - x) + _ϵ)
coupling_input(u::Function, T::Vector) = coupling_input(u.(T), T)
coupling_input(u::LinearInterpolation, T::Vector) = coupling_input(u.(T), T)

"""
    coupling_input(g::Gaussian)

Analytical input coupling ``g_u(t)`` for a Gaussian pulse.
"""
function coupling_input(g::Gaussian)
    mode, ∫m2 = _gaussian_mode(g)
    return PulseCoupling(t -> mode(t)' / √(abs(1 - ∫m2(t)) + _ϵ))
end

"""
    coupling_output(v, T)

Compute the virtual-cavity output coupling ``g_v(t)`` from an output mode `v(t)`
sampled on time grid `T`. Returns a [`PulseCoupling`](@ref) (callable as `g(t)`,
broadcastable as `g.(T)`) that plugs directly into a `to_numeric` `time_parameter`.
"""
coupling_output(v::Vector, T::Vector) = _coupling_from_mode(-v, T, x -> x + _ϵ)
coupling_output(v::Function, T::Vector) = coupling_output(v.(T), T)
coupling_output(v::LinearInterpolation, T::Vector) = coupling_output(v.(T), T)

"""
    coupling_output(g::Gaussian)

Analytical output coupling ``g_v(t)`` for a Gaussian pulse.
"""
function coupling_output(g::Gaussian)
    mode, ∫m2 = _gaussian_mode(g)
    return PulseCoupling(t -> -mode(t)' / √(∫m2(t) + _ϵ))
end

# ──────────────────────────────────────────────
# effective_output_mode (was v_eff)
# ──────────────────────────────────────────────

"""
    effective_output_mode(v_fcts, gv_fcts, T, i)
    effective_output_mode(v_fcts, T, i)

Compute the effective output mode ``v_i^{\\mathrm{eff}}(t)`` for a system with multiple
output modes, due to the pulse distortion from the preceding output cavities.
The output modes in `v_fcts` must be sorted starting with the first output cavity after the system.

All kwargs are passed on to the ODE solver.
"""
function effective_output_mode(v_fcts, gv_fcts, T, i; alg = Tsit5(), kwargs...)
    @assert i > 1
    n = i - 1
    # Capture as a tuple for concrete closure types
    gv = ntuple(k -> gv_fcts[k], n)
    v_i = v_fcts[i]
    function multiple_outputs_α!(dα, α, p, t)
        gv_buf = p
        @inbounds for k = 1:n
            gv_buf[k] = gv[k](t)
        end
        vi_t = v_i(t)
        @inbounds for j = 1:n
            coupling_sum = zero(ComplexF64)
            for k = 1:(j-1)
                coupling_sum += gv_buf[k]' * α[k]
            end
            dα[j] = -gv_buf[j] * (vi_t + coupling_sum) - 0.5 * abs2(gv_buf[j]) * α[j]
        end
    end
    u0 = zeros(ComplexF64, n)
    p = Vector{ComplexF64}(undef, n)
    tspan = (T[1], T[end])
    prob = ODEProblem(multiple_outputs_α!, u0, tspan, p)
    sol_α = solve(prob, alg; kwargs...)
    function v_i_eff(t)
        α_t = sol_α(t)
        result = v_i(t)
        @inbounds for k = 1:n
            result += gv[k](t)' * α_t[k]
        end
        return result
    end
    return v_i_eff
end

effective_output_mode(v_fcts, T, i; kwargs...) = effective_output_mode(
    v_fcts,
    [coupling_output(v_, T) for v_ in v_fcts],
    T,
    i;
    kwargs...,
)

function effective_output_mode(
    v_data::AbstractVector{<:AbstractVector},
    gv_data::AbstractVector{<:AbstractVector},
    T::AbstractVector,
    i;
    alg = Tsit5(),
    kwargs...,
)
    v_fcts = [_interpolate_mode(v_, T) for v_ in v_data]
    gv_fcts = [_interpolate_mode(gv_, T) for gv_ in gv_data]
    return effective_output_mode(v_fcts, gv_fcts, T, i; alg, kwargs...)
end

effective_output_mode(
    v_data::AbstractVector{<:AbstractVector},
    T::AbstractVector,
    i;
    alg = Tsit5(),
    kwargs...,
) = effective_output_mode(
    v_data,
    [coupling_output(v_, T).(T) for v_ in v_data],
    T,
    i;
    alg,
    kwargs...,
)

# ──────────────────────────────────────────────
# effective_input_mode (was u_eff)
# ──────────────────────────────────────────────

"""
    effective_input_mode(u_fcts, gu_fcts, T, i)
    effective_input_mode(u_fcts, T, i)

Compute the effective input mode ``u_i^{\\mathrm{eff}}(t)`` for a system with multiple
input modes, due to the pulse distortion from the subsequent input cavities.
The input modes in `u_fcts` must be sorted starting with the first input cavity before the system.

All kwargs are passed on to the ODE solver.
"""
function effective_input_mode(u_fcts, gu_fcts, T, i; alg = Tsit5(), kwargs...)
    @assert i > 1
    n = i - 1
    # Capture as a tuple for concrete closure types
    gu = ntuple(k -> gu_fcts[k], n)
    u_i = u_fcts[i]
    gu_i = gu_fcts[i]
    function multiple_inputs_α!(dα, α, p, t)
        gu_buf = p
        @inbounds for k = 1:n
            gu_buf[k] = gu[k](t)
        end
        ui_t = u_i(t)
        gui_t = gu_i(t)
        @inbounds for j = 1:n
            coupling_sum = zero(ComplexF64)
            for k = 1:(j-1)
                coupling_sum += gui_t' * α[k]
            end
            dα[j] = -gu_buf[j] * (ui_t - coupling_sum) + 0.5 * abs2(gu_buf[j]) * α[j]
        end
    end
    u0 = zeros(ComplexF64, n)
    p = Vector{ComplexF64}(undef, n)
    tspan = (T[1], T[end])
    prob = ODEProblem(multiple_inputs_α!, u0, tspan, p)
    sol_α = solve(prob, alg; kwargs...)
    function u_i_eff(t)
        α_t = sol_α(t)
        result = u_i(t)
        @inbounds for k = 1:n
            result -= gu[k](t)' * α_t[k]
        end
        return result
    end
    return u_i_eff
end

effective_input_mode(u_fcts, T, i; kwargs...) =
    effective_input_mode(u_fcts, [coupling_input(u_, T) for u_ in u_fcts], T, i; kwargs...)

function effective_input_mode(
    u_data::AbstractVector{<:AbstractVector},
    gu_data::AbstractVector{<:AbstractVector},
    T::AbstractVector,
    i;
    alg = Tsit5(),
    kwargs...,
)
    u_fcts = [_interpolate_mode(u_, T) for u_ in u_data]
    gu_fcts = [_interpolate_mode(gu_, T) for gu_ in gu_data]
    return effective_input_mode(u_fcts, gu_fcts, T, i; alg, kwargs...)
end

effective_input_mode(
    u_data::AbstractVector{<:AbstractVector},
    T::AbstractVector,
    i;
    alg = Tsit5(),
    kwargs...,
) = effective_input_mode(
    u_data,
    [coupling_input(u_, T).(T) for u_ in u_data],
    T,
    i;
    alg,
    kwargs...,
)

# ──────────────────────────────────────────────
# coupling_delay_out / coupling_delay_in
# ──────────────────────────────────────────────

function _delay_coupling_from_modes(num_mode::Vector, u::Vector, v::Vector, T::Vector)
    nt = length(T)
    mode_sq = Vector{Float64}(undef, nt)
    map!(abs2, mode_sq, u)
    ∫u2 = cumul_integrate(T, mode_sq)
    map!(abs2, mode_sq, v)
    ∫v2 = cumul_integrate(T, mode_sq)
    coupling = zeros(ComplexF64, nt)
    @inbounds for i = 1:nt
        denom = abs(∫v2[i] - ∫u2[i])
        if sqrt(abs(denom)) > _tol_div
            coupling[i] = num_mode[i]' / sqrt(denom + _ϵ)
        end
    end
    return PulseCoupling(LinearInterpolation(coupling, T; extrapolation = _extrapolate))
end

"""
    coupling_delay_out(u, v, T)

Compute the out-coupling strength for a delay cavity.
Returns a [`PulseCoupling`](@ref) (callable as `g(t)`, broadcastable as `g.(T)`).
"""
coupling_delay_out(u::Vector, v::Vector, T::Vector) = _delay_coupling_from_modes(u, u, v, T)
coupling_delay_out(u::Function, v::Function, T::Vector) =
    coupling_delay_out(u.(T), v.(T), T)
coupling_delay_out(u::LinearInterpolation, v::LinearInterpolation, T::Vector) =
    coupling_delay_out(u.(T), v.(T), T)

"""
    coupling_delay_in(u, v, T)

Compute the in-coupling strength for a delay cavity.
Returns a [`PulseCoupling`](@ref) (callable as `g(t)`, broadcastable as `g.(T)`).
"""
coupling_delay_in(u::Vector, v::Vector, T::Vector) = _delay_coupling_from_modes(-v, u, v, T)
coupling_delay_in(u::Function, v::Function, T::Vector) = coupling_delay_in(u.(T), v.(T), T)
coupling_delay_in(u::LinearInterpolation, v::LinearInterpolation, T::Vector) =
    coupling_delay_in(u.(T), v.(T), T)
