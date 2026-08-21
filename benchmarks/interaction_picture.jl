"""
Interaction picture benchmarks — computing coefficient matrices for mode transformations.
"""
function benchmark_interaction_picture!(SUITE)
    SUITE["Interaction Picture"] = BenchmarkGroup()

    ## Setup: u -> TLS -> v system
    γ_ = 1.0
    τ = 1 / γ_
    t_p = 4 / γ_
    u(t) = 1 / (sqrt(τ) * π^(1 / 4)) * exp(-0.5 * ((t - t_p) / τ)^2)

    T_end = 12.0
    T = [0:0.005:1;] * T_end

    gu_t = coupling_input(u, T)
    gv_t = coupling_output(u, T)
    A_uv = coupling_matrix((gu_t, gv_t))

    ## --- Coupling matrix A(t) evaluation (called every ODE step) ---

    SUITE["Interaction Picture"]["coupling matrix evaluation"] = BenchmarkGroup()

    t_mid = T[length(T)÷2]

    # 2-mode A(t)
    SUITE["Interaction Picture"]["coupling matrix evaluation"]["2 modes"] =
        @benchmarkable $A_uv($t_mid)

    # 4-mode A(t)
    u2(t) = 1 / (sqrt(τ) * π^(1 / 4)) * exp(-0.5 * ((t - t_p * 1.5) / τ)^2)
    g3_t = coupling_input(u2, T)
    g4_t = coupling_output(u2, T)
    A_4m = coupling_matrix((gu_t, gv_t, g3_t, g4_t))

    SUITE["Interaction Picture"]["coupling matrix evaluation"]["4 modes"] =
        @benchmarkable $A_4m($t_mid)

    ## --- Coefficient matrix M(t) ---

    SUITE["Interaction Picture"]["coefficient matrix M"] = BenchmarkGroup()

    # Allocation-heavy ODE solve: GC before each sample to stabilize the timing
    # (see the note in runbenchmarks.jl).
    SUITE["Interaction Picture"]["coefficient matrix M"]["numerical (ODE)"] =
        @benchmarkable solve_mode_evolution($A_uv, $T)

    SUITE["Interaction Picture"]["coefficient matrix M"]["analytical (2 equal modes)"] =
        @benchmarkable solve_mode_evolution_symmetric($u, $T)

    ## --- Symbolic operator substitution ---

    hu = FockSpace(:u)
    hs = NLevelSpace(:s, 2)
    hv = FockSpace(:v)
    h = hu ⊗ hs ⊗ hv

    au_sym = Destroy(h, :a_u, 1)
    av_sym = Destroy(h, :a_v, 3)
    σ_sym = Transition(h, :σ, 1, 2, 2)

    @variables gu_sym::Real γ_sym::Real gv_sym::Real

    G_u = SLH(1, gu_sym' * au_sym, 0)
    G_s = SLH(1, sqrt(γ_sym) * σ_sym, 0)
    G_v = SLH(1, gv_sym' * av_sym, 0)
    G_cas = ▷(G_u, G_s, G_v)

    H = hamiltonian(G_cas)
    J = jump_operator(G_cas)[1]
    H_uv = hamiltonian(▷(G_u, G_v))
    H_int_ = simplify(H - H_uv)

    M_sym(i, j) = Symbolics.variable(Symbol("M_{$(i)$(j)}"); T = Complex{Real})
    a0_ls = [au_sym, av_sym]
    la = length(a0_ls)
    a_int_ls = [sum(M_sym(i, j) * a0_ls[j] for j = 1:la) for i = 1:la]
    int_dict = Dict(a0_ls .=> a_int_ls)

    SUITE["Interaction Picture"]["operator substitution"] = BenchmarkGroup()

    SUITE["Interaction Picture"]["operator substitution"]["TLS cascade"] =
        @benchmarkable begin
            simplify(substitute($H_int_, $int_dict))
            simplify(substitute($J, $int_dict))
        end

    return nothing
end
