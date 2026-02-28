######################################
### functions for virtual cavities ###
######################################

_tol_div = 1e-10 #12
_extrapolate = ExtrapolationType.Extension
_ϵu = 1e-10
_ϵv = 1e-10

"""
    u_to_gu(u, T)

Compute the virtual-cavity coupling ``g_u(t)`` from an input mode `u(t)` sampled on `T`.
Returns a `LinearInterpolation` over `T`.
"""
function u_to_gu(u::Vector, T::Vector)
    l_T = length(T)
    ∫u2_t = cumul_integrate(T, abs2.(u)) .+ 0im
    gu_t = zeros(ComplexF64, l_T)
    for i = 2:l_T
        if abs(sqrt(1 - ∫u2_t[i])) > _tol_div
            gu_t[i] = u[i]' / sqrt(1 - ∫u2_t[i] + _ϵu)
        end #else 0
        # gu_t[i] = u[i]' / max(sqrt(1 - ∫u2_t[i]), _tol_div) # TODO: test
    end
    return LinearInterpolation(gu_t, T; extrapolation = _extrapolate)
end
u_to_gu(u::Function, T::Vector) = u_to_gu(u.(T), T)
u_to_gu(u::LinearInterpolation, T::Vector) = u_to_gu(u.(T), T)

"""
    v_to_gv(v, T)

Compute the virtual-cavity coupling ``g_v(t)`` from an output mode `v(t)` sampled on `T`.
Returns a `LinearInterpolation` over `T`.
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
    return LinearInterpolation(gv_t, T; extrapolation = _extrapolate)
end
v_to_gv(v::Function, T::Vector) = v_to_gv(v.(T), T)
v_to_gv(v::LinearInterpolation, T::Vector) = v_to_gv(v.(T), T)

"""
    u_to_gu_Gauss(u, τ, σ)

Compute ``g_u(t)`` for a Gaussian input mode `u(t)` with delay `τ` and width `σ`.
Returns a callable `t -> g_u(t)`.
"""
function u_to_gu_Gauss(u, τ, σ)
    # u = mode function (Gauss or Gauss*exp(i*ω*t))
    # τ = time delay # σ = width
    ∫u_2_t(t_) = 0.5 * (erf((t_ - τ) / σ) + erf(τ / σ))
    f(t_) = u(t_)' / max(√(1 - ∫u_2_t(t_)), _tol_div) # TODO: benchmark
    return f
end
"""
    v_to_gv_Gauss(v, τ, σ)

Compute ``g_v(t)` for a Gaussian output mode `v(t)` with delay `τ` and width `σ`.
Returns a callable `t -> g_v(t)`.
"""
function v_to_gv_Gauss(v, τ, σ)
    # v = mode function (Gauss or Gauss*exp(i*ω*t))
    # τ = time delay # σ = width
    ∫v_2_t(t_) = 0.5 * (erf((t_ - τ) / σ) + erf(τ / σ))
    f(t_) = v(t_)' / max(√(∫v_2_t(t_)), _tol_div) # TODO: benchmark
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
