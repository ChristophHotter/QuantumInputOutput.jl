######################################
### functions for virtual cavities ###
######################################

_tol_div = 1e-10 #12
_extrapolate = ExtrapolationType.Extension
_ϵu = 1e-10
_ϵv = 1e-10

"""
    u_to_gu(u, T)

Compute the virtual-cavity coupling ``g_u(t)`` from an input mode `u(t)` sampled on `T` with a linear interpolation.
Returns a callable `t -> g_u(t)`.
"""
function u_to_gu(u::Vector, T::Vector)
    l_T = length(T)
    ∫u2_t = cumul_integrate(T, abs2.(u)) .+ 0im
    gu_t = zeros(ComplexF64, l_T)
    for i = 2:l_T
        if abs(sqrt(1 - ∫u2_t[i])) > _tol_div
            gu_t[i] = u[i]' / sqrt(abs(1 - ∫u2_t[i]) + _ϵu)
        end #else 0
        # gu_t[i] = u[i]' / max(sqrt(1 - ∫u2_t[i]), _tol_div) # TODO: test
    end
    gu_int = LinearInterpolation(gu_t, T; extrapolation = _extrapolate)
    return t -> gu_int(t)
end
u_to_gu(u::Function, T::Vector) = u_to_gu(u.(T), T)
u_to_gu(u::LinearInterpolation, T::Vector) = u_to_gu(u.(T), T)

"""
    v_to_gv(v, T)

Compute the virtual-cavity coupling ``g_v(t)`` from an output mode `v(t)` sampled on `T` with a linear interpolation.
Returns a callable `t -> g_v(t)`.
"""
function v_to_gv(v::Vector, T::Vector)
    l_T = length(T)
    ∫v2_t = cumul_integrate(T, abs2.(v)) .+ 0im
    gv_t = zeros(ComplexF64, l_T)
    for i = 2:l_T
        if abs(sqrt(∫v2_t[i])) > _tol_div
            gv_t[i] = -v[i]' / sqrt(∫v2_t[i] + _ϵv)
        end # else 0
        # gv_t[i] = -v[i]' / max(sqrt( ∫v2_t[i]), _tol_div) # TODO: test
    end
    gv_int = LinearInterpolation(gv_t, T; extrapolation = _extrapolate)
    return t -> gv_int(t)
end
v_to_gv(v::Function, T::Vector) = v_to_gv(v.(T), T)
v_to_gv(v::LinearInterpolation, T::Vector) = v_to_gv(v.(T), T)

"""
    u_to_gu_Gauss(τ, σ; δ=0)

Compute ``g_u(t)`` for a Gaussian input mode `u(t)` with delay `τ`, width `σ`, and detuning `δ`:

`` u(t) = 1/\\sqrt{ \\sigma*\\sqrt{\\pi} }*e^{-(t - \\tau)^2 / 2 \\sigma^2 )} * e^{i*\\delta*t} ``


Returns a callable `t -> g_u(t)`.
"""
function u_to_gu_Gauss(τ, σ; δ=0) # slower than u_to_gu
    if δ==0 # type stable
        u = t -> 1/(√(σ)*π^(1/4)) *exp( -0.5*(t - τ)^2/σ^2 )
    else
        u = t -> 1/(√(σ)*π^(1/4)) *exp( -0.5*(t - τ)^2/σ^2 ) * exp(1im*δ*t)
    end
    ∫u_2_t(t_) = 0.5 * (erf((t_ - τ) / σ) + erf(τ / σ))
    f = t_ -> u(t_)' / √(abs(1 - ∫u_2_t(t_)) + _ϵu) 
    return f
end
"""
    v_to_gv_Gauss(τ, σ; δ=0)

Compute ``g_v(t)` for a Gaussian output mode `v(t)` with delay `τ`, width `σ`, and detuning `δ`:

`` v(t) = 1/\\sqrt{ \\sigma*\\sqrt{\\pi} }*e^{-(t - \\tau)^2 / 2 \\sigma^2 )} * e^{i*\\delta*t} ``

Returns a callable `t -> g_v(t)`.
"""
function v_to_gv_Gauss(τ, σ; δ=0) # slower than v_to_gv
    if δ==0 # type stable
        v = t -> 1/(√(σ)*π^(1/4)) *exp( -0.5*(t - τ)^2/σ^2 )
    else
        v = t -> 1/(√(σ)*π^(1/4)) *exp( -0.5*(t - τ)^2/σ^2 ) * exp(1im*δ*t)
    end
    ∫v_2_t(t_) = 0.5 * (erf((t_ - τ) / σ) + erf(τ / σ)) 
    f = t_ -> -v(t_)' / √(∫v_2_t(t_) + _ϵv) 
    return f
end

"""
    vi_to_v_i_im1(v_fcts, gv_fcts, T_ls, i)
    vi_to_v_i_im1(v_fcts, T_ls, i)

Compute the effective output mode `v_i^{(i-1)}(t)` for multiple output modes.
"""
function vi_to_v_i_im1(v_fcts, gv_fcts, T_ls, i)
    @assert i > 1
    function multiple_outputs_α(dα, α, p, t) # only for i>1
        for j = 1:i-1
            dα[j] = -gv_fcts[j](t) * (v_fcts[i](t) +
                sum((gv_fcts[k](t))' * α[k] for k = 1:j-1; init = 0.0)) -
                0.5 * abs2(gv_fcts[j](t)) * α[j]
        end
    end
    u0 = zeros(ComplexF64, i-1)
    tspan = (T_ls[1], T_ls[end])
    prob = ODEProblem(multiple_outputs_α, u0, tspan)
    sol_α = solve(prob) #; abstol, reltol) # TODO: kwarg
    v_i_im1(t) = v_fcts[i](t) + sum((gv_fcts[k](t))' * sol_α(t)[k] for k = 1:i-1)
    return v_i_im1
end
vi_to_v_i_im1(v_fcts, T_ls, i) = vi_to_v_i_im1(v_fcts, [v_to_gv(v_, T_ls) for v_ in v_fcts], T_ls, i)

"""
    ui_to_u_i_im1(u_fcts, gu_fcts, T_ls, i)
    ui_to_u_i_im1(u_fcts, T_ls, i)

Compute the effective input mode `u_i^{(i-1)}(t)` for multiple input modes.
"""
function ui_to_u_i_im1(u_fcts, gu_fcts, T_ls, i)
    @assert i > 1
    function multiple_inputs_α(dα, α, p, t) # only for i>1
        for j = 1:i-1
            dα[j] = -gu_fcts[j](t) * (u_fcts[i](t) -
                sum((gu_fcts[i](t))' * α[k] for k = 1:j-1; init = 0.0)) +
                0.5 * abs2(gu_fcts[j](t)) * α[j]
        end
        # 2 Typos in PRA2020-Kiilerich Eq.(A15): last term "-" → "+" and |g_ui|² → |g_uj|²
    end
    u0 = zeros(ComplexF64, i-1)
    tspan = (T_ls[1], T_ls[end])
    prob = ODEProblem(multiple_inputs_α, u0, tspan)
    sol_α = solve(prob, Tsit5()) #; abstol, reltol) # TODO: kwarg
    u_i_im1(t) = u_fcts[i](t) - sum((gu_fcts[k](t))' * sol_α(t)[k] for k = 1:i-1)
    return u_i_im1
end
ui_to_u_i_im1(u_fcts, T_ls, i) = ui_to_u_i_im1(u_fcts, [u_to_gu(u_, T_ls) for u_ in u_fcts], T_ls, i)
# TODO: rename, better method (Victor)

"""
    uv_to_gout(u, v, T)

Compute the out-coupling strength ``\\tilde g_{\\mathrm{out},u,v}(t)`` for a delay cavity
that absorbs an incoming pulse `v(t)` while simultaneously emitting the desired pulse `u(t)`.
Returns callable `t -> g_out(t)` based on samples on `T` with a linear interpolation.
"""
function uv_to_gout(u::Vector, v::Vector, T::Vector)
    l_T = length(T)
    ∫u2_t = cumul_integrate(T, abs2.(u)) .+ 0im
    ∫v2_t = cumul_integrate(T, abs2.(v)) .+ 0im
    gout_t = zeros(ComplexF64, l_T)
    for i = 2:l_T
        denom = abs(∫v2_t[i] - ∫u2_t[i])
        if abs(sqrt(denom)) > _tol_div
            gout_t[i] = u[i]' / sqrt(denom + _ϵu)
        end
    end
    gout_int = LinearInterpolation(gout_t, T; extrapolation = _extrapolate)
    return t -> gout_int(t)
end
uv_to_gout(u::Function, v::Function, T::Vector) = uv_to_gout(u.(T), v.(T), T)
uv_to_gout(u::LinearInterpolation, v::LinearInterpolation, T::Vector) = uv_to_gout(u.(T), v.(T), T)

"""
    uv_to_gin(u, v, T)

Compute the in-coupling strength ``\\tilde g_{\\mathrm{in},v,u}(t)`` for a delay cavity
that absorbs an incoming pulse `v(t)` while simultaneously emitting the desired pulse `u(t)`.
Returns callable `t -> g_in(t)` based on samples on `T` with a linear interpolation. 
"""
function uv_to_gin(u::Vector, v::Vector, T::Vector)
    l_T = length(T)
    ∫u2_t = cumul_integrate(T, abs2.(u)) .+ 0im
    ∫v2_t = cumul_integrate(T, abs2.(v)) .+ 0im
    gin_t = zeros(ComplexF64, l_T)
    for i = 2:l_T
        denom = abs(∫v2_t[i] - ∫u2_t[i])
        if abs(sqrt(denom)) > _tol_div
            gin_t[i] = -v[i]' / sqrt(denom + _ϵv)
        end
    end
    gin_int = LinearInterpolation(gin_t, T; extrapolation = _extrapolate)
    return t -> gin_int(t)
end
uv_to_gin(u::Function, v::Function, T::Vector) = uv_to_gin(u.(T), v.(T), T)
uv_to_gin(u::LinearInterpolation, v::LinearInterpolation, T::Vector) = uv_to_gin(u.(T), v.(T), T)
