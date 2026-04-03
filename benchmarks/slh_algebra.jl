"""
SLH algebra benchmarks — symbolic and numeric composition of quantum networks.
"""
function benchmark_slh_algebra!(SUITE)
    SUITE["SLH Algebra"] = BenchmarkGroup()

    ## --- Symbolic SLH (SecondQuantizedAlgebra operators) ---

    SUITE["SLH Algebra"]["symbolic"] = BenchmarkGroup()

    # Setup: 3-cavity cascade (example 01-1)
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

    SUITE["SLH Algebra"]["symbolic"]["3-cavity cascade"] =
        @benchmarkable begin
            G_cas = ▷($G_u, $G_c, $G_v)
            get_hamiltonian(G_cas)
            get_lindblad(G_cas)
        end

    SUITE["SLH Algebra"]["symbolic"]["concatenate + cascade"] =
        @benchmarkable begin
            G_R = $G_u ▷ $G_c
            G_L = $G_v
            G_R ⊞ G_L
        end

    # Feedback: coherent-feedback OPO loop (from test_feedback.jl)
    hs = FockSpace(:s)
    a_fb = Destroy(hs, :a, 1)
    κ_fb = 0.7
    ϵ_fb = 0.45
    η_fb = 0.65
    r_fb = simplify(√(1 - η_fb^2))

    G_opo = SLH(1, √(κ_fb) * a_fb, 1im * ϵ_fb * (a_fb'^2 - a_fb^2))
    G_bs = SLH([-r_fb η_fb; η_fb r_fb], [0, 0], 0)
    G_fb_unconnected = G_opo ⊞ G_bs

    SUITE["SLH Algebra"]["symbolic"]["feedback OPO loop"] =
        @benchmarkable feedback($G_fb_unconnected, 1 => 2, 2 => 1)

    ## --- Numeric SLH (SLHqo with QuantumOptics operators) ---

    SUITE["SLH Algebra"]["numeric (SLHqo)"] = BenchmarkGroup()

    # Setup: 2-QD bidirectional waveguide (example 05-2, realistic basis)
    N = 2
    bu = FockBasis(20)
    ba = NLevelBasis(2)
    b_qds = tensor([ba for _ in 1:N]...)
    b = bu ⊗ b_qds

    σ(i, j, k) = embed(b, i + 1, transition(ba, j, k))
    a_u = destroy(bu) ⊗ one(b_qds)

    γ_val = 1.0
    β = 0.9
    γRn = fill(γ_val * β / 2, N)
    γLn = fill(γ_val * β / 2, N)
    ϕn = [π / 10]

    σt = 0.8
    t0 = 4σt
    T = [0:0.005:1;] * 3t0
    u1(t) = 1 / (sqrt(σt) * π^(1 / 4)) * exp(-(t - t0)^2 / (2 * σt^2))
    gu_t = u_to_gu(u1, T)

    G_u_qo = SLHqo(1, t -> gu_t(t) * a_u, 0 * one(b))
    G_ϕ_12 = SLHqo(exp(1im * ϕn[1]), 0 * one(b), 0 * one(b))
    G_R1 = SLHqo(1, √(γRn[1]) * σ(1, 1, 2), 0 * one(b))
    G_R2 = SLHqo(1, √(γRn[2]) * σ(2, 1, 2), 0 * one(b))
    G_L1 = SLHqo(1, √(γLn[1]) * σ(1, 1, 2), 0 * one(b))
    G_L2 = SLHqo(1, √(γLn[2]) * σ(2, 1, 2), 0 * one(b))

    SUITE["SLH Algebra"]["numeric (SLHqo)"]["2-QD waveguide composition"] =
        @benchmarkable begin
            G_R_t = $G_u_qo ▷ $G_R1 ▷ $G_ϕ_12 ▷ $G_R2
            G_L_t = $G_L2 ▷ $G_ϕ_12 ▷ $G_L1
            G_t = G_R_t ⊞ G_L_t
            get_hamiltonian(G_t)
            get_lindblad(G_t)
        end

    ## --- Time-dependent closure evaluation (the ODE hot loop) ---

    SUITE["SLH Algebra"]["closure evaluation"] = BenchmarkGroup()

    # Build the composed network once, then benchmark calling H(t) and L(t)
    G_R_t = G_u_qo ▷ G_R1 ▷ G_ϕ_12 ▷ G_R2
    G_L_t = G_L2 ▷ G_ϕ_12 ▷ G_L1
    G_t = G_R_t ⊞ G_L_t

    H_f = get_hamiltonian(G_t)
    L_f = get_lindblad(G_t)

    t_mid = T[length(T) ÷ 2]

    SUITE["SLH Algebra"]["closure evaluation"]["2-QD waveguide H(t)"] =
        @benchmarkable $H_f($t_mid)

    SUITE["SLH Algebra"]["closure evaluation"]["2-QD waveguide L(t)"] =
        @benchmarkable begin
            for Li in $L_f
                Li($t_mid)
            end
        end

    return nothing
end
