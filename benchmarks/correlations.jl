"""
Correlation function benchmarks — two-time correlation matrices.
Based on example 01-1: single photon cavity scattering.
"""
function benchmark_correlations!(SUITE)
    SUITE["Correlations"] = BenchmarkGroup()

    ## Setup: single photon scattering through a cavity
    hu1 = FockSpace(:u1)
    hc1 = FockSpace(:c1)
    hv1 = FockSpace(:v1)
    h = hu1 ⊗ hc1 ⊗ hv1

    au = Destroy(h, :a_u, 1)
    c = Destroy(h, :c, 2)
    av = Destroy(h, :a_v, 3)

    @variables g_u::Real Δ::Real γ::Real g_v::Number

    G_u = SLH(1, g_u * au, 0)
    G_c = SLH(1, √(γ) * c, Δ * c'c)
    G_v = SLH(1, g_v * av, 0)
    G_cas = ▷(G_u, G_c, G_v)

    H_sym = hamiltonian(G_cas)
    L_sym = lindblad(G_cas)[1]

    γ_ = 1.0
    σ_pulse = 1 / γ_
    T = [0:0.002:1;] * 12σ_pulse
    u1(t) = 1 / (sqrt(σ_pulse) * π^(1 / 4)) * exp(-(t - 4σ_pulse)^2 / (2 * σ_pulse^2))
    gu_t = coupling_input(u1, T)

    bu1 = FockBasis(1)
    bc1 = FockBasis(1)
    bv1 = FockBasis(1)
    b = bu1 ⊗ bc1 ⊗ bv1

    dict_p = Dict([γ, Δ, g_v] .=> [γ_, 0.0, 0])
    dict_p_t = Dict(g_u => gu_t)

    H_QO = to_numeric(H_sym, b; parameter = dict_p, time_parameter = dict_p_t)
    L_QO = to_numeric(L_sym, b; parameter = dict_p, time_parameter = dict_p_t)

    ψ0 = fockstate(bu1, 1) ⊗ fockstate(bc1, 0) ⊗ fockstate(bv1, 0)
    _, ρt = timeevolution.master_dynamic(T, ψ0, H_QO, [L_QO])

    au_qo = to_numeric(au, b)
    c_qo = to_numeric(c, b)
    Ls(t) = gu_t(t) * au_qo + √(γ_) * c_qo

    ## --- Two-time correlation ---
    # Pass the time-dependent operators (`H_QO`, `[L_QO]`) directly to the solver rather
    # than wrapping them in a `(t, ρ)` closure: since v0.10 `to_numeric` returns a lazy
    # `TimeDependentSum`, the direct path lets the solver build its integrator once and is
    # markedly faster than re-reading the operator from a function at every step.

    SUITE["Correlations"]["two-time"] = BenchmarkGroup()

    SUITE["Correlations"]["two-time"]["single photon cavity"] =
        @benchmarkable correlation_matrix($T, $ρt, $H_QO, [$L_QO], $Ls)
    return nothing
end
