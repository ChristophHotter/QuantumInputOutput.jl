"""
Translation benchmarks — converting symbolic expressions to numeric QuantumOptics operators.
"""
function benchmark_translation!(SUITE)
    SUITE["Translation"] = BenchmarkGroup()

    ## Symbolic setup: atom-cavity system
    @variables κ_R::Real κ_L::Real Δ_sym::Real
    @variables E::Number

    hc = FockSpace(:cavity)
    ha_ = NLevelSpace(Symbol("a"), 2)
    h = hc ⊗ ha_

    a = Destroy(h, :a, 1)
    σ(i, j) = Transition(h, Symbol("σ"), i, j, 2)

    ## Derive H and L from cascade
    @variables gu_sym::Real γ_sym::Real gv_sym::Real
    hu_ = FockSpace(:u)
    hc_ = FockSpace(:c)
    hv_ = FockSpace(:v)
    h3 = hu_ ⊗ hc_ ⊗ hv_
    au_ = Destroy(h3, :a_u, 1)
    c_ = Destroy(h3, :c, 2)
    av_ = Destroy(h3, :a_v, 3)

    G_u = SLH(1, gu_sym * au_, 0)
    G_c = SLH(1, √(γ_sym) * c_, Δ_sym * c_'c_)
    G_v = SLH(1, gv_sym * av_, 0)
    G_cas = ▷(G_u, G_c, G_v)
    H_sym = hamiltonian(G_cas)
    L_sym = lindblad(G_cas)[1]

    ## --- Static translation (no time dependence) ---

    SUITE["Translation"]["static"] = BenchmarkGroup()

    bc = FockBasis(4)
    ba = NLevelBasis(2)
    b_multi = bc ⊗ ba

    dict_p_static = Dict([κ_R, κ_L, Δ_sym] .=> [1.5, 1.0, 0.2])

    expr_composite = a * 3 + Δ_sym * σ(2, 2)

    SUITE["Translation"]["static"]["atom-cavity"] =
        @benchmarkable to_numeric($expr_composite, $b_multi; parameter = $dict_p_static)

    ## --- Time-dependent translation ---

    SUITE["Translation"]["time-dependent"] = BenchmarkGroup()

    E_t(t) = 2 * t + 1im
    dict_p_t = Dict(E => E_t)

    expr_td = a * 3 * conj(E) + Δ_sym * σ(2, 2)

    SUITE["Translation"]["time-dependent"]["atom-cavity"] = @benchmarkable to_numeric(
        $expr_td,
        $b_multi;
        parameter = $dict_p_static,
        time_parameter = $dict_p_t,
    )

    # Full cascade H and L translation
    γ_ = 1.0
    σ_pulse = 1 / γ_
    T = [0:0.002:1;] * 12σ_pulse
    u_pulse(t) = 1 / (sqrt(σ_pulse) * π^(1 / 4)) * exp(-(t - 4σ_pulse)^2 / (2 * σ_pulse^2))
    gu_t = coupling_input(u_pulse, T)
    gv_t = coupling_output(u_pulse, T)

    bu = FockBasis(20)
    bc3 = FockBasis(6)
    bv = FockBasis(6)
    b_cav = bu ⊗ bc3 ⊗ bv

    dict_p_cav = Dict([γ_sym, Δ_sym, gv_sym] .=> [γ_, 0.0, 0])
    dict_p_t_cav = Dict([gu_sym, gv_sym] .=> [gu_t, gv_t])

    SUITE["Translation"]["time-dependent"]["3-cavity H+L"] = @benchmarkable begin
        to_numeric($H_sym, $b_cav; parameter = $dict_p_cav, time_parameter = $dict_p_t_cav)
        to_numeric(
            $L_sym,
            $b_cav;
            parameter = $dict_p_cav,
            time_parameter = $dict_p_t_cav,
        )
    end

    ## --- Closure evaluation (the ODE hot loop) ---

    SUITE["Translation"]["closure evaluation"] = BenchmarkGroup()

    H_QO = to_numeric(H_sym, b_cav; parameter = dict_p_cav, time_parameter = dict_p_t_cav)
    L_QO = to_numeric(L_sym, b_cav; parameter = dict_p_cav, time_parameter = dict_p_t_cav)

    t_mid = T[length(T)÷2]

    SUITE["Translation"]["closure evaluation"]["3-cavity H(t)"] =
        @benchmarkable $H_QO($t_mid)

    SUITE["Translation"]["closure evaluation"]["3-cavity L(t)"] =
        @benchmarkable $L_QO($t_mid)

    return nothing
end
