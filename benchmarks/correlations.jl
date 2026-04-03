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

    gu, Δ, γ = rnumbers("g_u Δ γ")
    gv = cnumber("g_v")

    G_u = SLH(1, gu * au, 0)
    G_c = SLH(1, √(γ) * c, Δ * c'c)
    G_v = SLH(1, gv * av, 0)
    G_cas = ▷(G_u, G_c, G_v)

    H_sym = get_hamiltonian(G_cas)
    L_sym = get_lindblad(G_cas)[1]

    γ_ = 1.0
    σ_pulse = 1 / γ_
    T = [0:0.002:1;] * 12σ_pulse
    u1(t) = 1 / (sqrt(σ_pulse) * π^(1 / 4)) * exp(-(t - 4σ_pulse)^2 / (2 * σ_pulse^2))
    gu_t = u_to_gu(u1, T)

    bu1 = FockBasis(1)
    bc1 = FockBasis(1)
    bv1 = FockBasis(1)
    b = bu1 ⊗ bc1 ⊗ bv1

    dict_p = Dict([γ, Δ, gv] .=> [γ_, 0.0, 0])
    dict_p_t = Dict(gu => gu_t)

    H_QO = translate_qo(H_sym, b; parameter = dict_p, time_parameter = dict_p_t)
    L_QO = translate_qo(L_sym, b; parameter = dict_p, time_parameter = dict_p_t)

    function input_output(t, ρ)
        Ht = H_QO(t)
        J = [L_QO(t)]
        return Ht, J, dagger.(J)
    end

    ψ0 = fockstate(bu1, 1) ⊗ fockstate(bc1, 0) ⊗ fockstate(bv1, 0)
    _, ρt = timeevolution.master_dynamic(T, ψ0, input_output)

    au_qo = translate_qo(au, b)
    c_qo = translate_qo(c, b)
    Ls(t) = gu_t(t) * au_qo + √(γ_) * c_qo

    ## --- Two-time correlation ---

    SUITE["Correlations"]["two-time"] = BenchmarkGroup()

    SUITE["Correlations"]["two-time"]["single photon cavity"] =
        @benchmarkable two_time_corr_matrix($T, $ρt, $input_output, $Ls)

    return nothing
end
