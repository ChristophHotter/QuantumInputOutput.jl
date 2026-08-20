"""
SLH algebra benchmarks — symbolic and numeric composition of quantum networks.
"""
function benchmark_slh_algebra!(SUITE)
    SUITE["SLH Algebra"] = BenchmarkGroup()

    ## --- Symbolic SLH (SecondQuantizedAlgebra operators) ---

    SUITE["SLH Algebra"]["symbolic"] = BenchmarkGroup()

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

    SUITE["SLH Algebra"]["symbolic"]["3-cavity cascade"] = @benchmarkable begin
        G_cas = ▷($G_u, $G_c, $G_v)
        hamiltonian(G_cas)
        jump_operator(G_cas)
    end

    SUITE["SLH Algebra"]["symbolic"]["concatenate + cascade"] = @benchmarkable begin
        G_R = $G_u ▷ $G_c
        G_L = $G_v
        G_R ⊞ G_L
    end

    # Feedback: coherent-feedback OPO loop
    hs = FockSpace(:s)
    a_fb = Destroy(hs, :a)
    κ_fb = 0.7
    ϵ_fb = 0.45
    η_fb = 0.65
    r_fb = simplify(√(1 - η_fb^2))

    G_opo = SLH(1, √(κ_fb) * a_fb, 1im * ϵ_fb * (a_fb'^2 - a_fb^2))
    G_bs = SLH([-r_fb η_fb; η_fb r_fb], [0, 0], 0)
    G_fb_unconnected = G_opo ⊞ G_bs

    SUITE["SLH Algebra"]["symbolic"]["feedback OPO loop"] =
        @benchmarkable feedback($G_fb_unconnected, 1 => 2, 2 => 1)

    ## --- Numeric SLH (unified SLH with QuantumOptics operators) ---

    SUITE["SLH Algebra"]["numeric"] = BenchmarkGroup()

    N = 2
    bu = FockBasis(20)
    ba = NLevelBasis(2)
    b_qds = tensor([ba for _ = 1:N]...)
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
    gu_t = coupling_input(u1, T)

    G_u_qo = SLH(1, t -> gu_t(t) * a_u, 0 * one(b))
    G_ϕ_12 = SLH(exp(1im * ϕn[1]), 0 * one(b), 0 * one(b))
    G_R1 = SLH(1, √(γRn[1]) * σ(1, 1, 2), 0 * one(b))
    G_R2 = SLH(1, √(γRn[2]) * σ(2, 1, 2), 0 * one(b))
    G_L1 = SLH(1, √(γLn[1]) * σ(1, 1, 2), 0 * one(b))
    G_L2 = SLH(1, √(γLn[2]) * σ(2, 1, 2), 0 * one(b))

    SUITE["SLH Algebra"]["numeric"]["2-QD waveguide composition"] = @benchmarkable begin
        G_R_t = $G_u_qo ▷ $G_R1 ▷ $G_ϕ_12 ▷ $G_R2
        G_L_t = $G_L2 ▷ $G_ϕ_12 ▷ $G_L1
        G_t = G_R_t ⊞ G_L_t
        hamiltonian(G_t)
        jump_operator(G_t)
    end

    ## --- Time-dependent closure evaluation (the ODE hot loop) ---

    SUITE["SLH Algebra"]["closure evaluation"] = BenchmarkGroup()

    G_R_t = G_u_qo ▷ G_R1 ▷ G_ϕ_12 ▷ G_R2
    G_L_t = G_L2 ▷ G_ϕ_12 ▷ G_L1
    G_t = G_R_t ⊞ G_L_t

    H_f = hamiltonian(G_t)
    L_f = jump_operator(G_t)

    t_mid = T[length(T)÷2]

    _callable(x) = x isa Union{Function,FunctionWrappers.FunctionWrapper}
    Hf = _callable(H_f) ? H_f : (t -> H_f)
    L_callables = [_callable(Li) ? Li : (t -> Li) for Li in L_f]

    SUITE["SLH Algebra"]["closure evaluation"]["2-QD waveguide H(t)"] =
        @benchmarkable $Hf($t_mid)

    SUITE["SLH Algebra"]["closure evaluation"]["2-QD waveguide L(t)"] = @benchmarkable begin
        for Li in $L_callables
            Li($t_mid)
        end
    end

    return nothing
end
