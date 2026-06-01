using QuantumInputOutput
using QuantumInputOutput: dagger
using SecondQuantizedAlgebra
using QuantumOptics
using FunctionWrappers: FunctionWrapper
using Test

@testset "compare_QDs_examples_05_1_05_2" begin
    N = 2

    # Common pulse and time grid
    σt = 0.8
    α0 = √(0.05)
    t0 = 4σt
    Tend = 3t0
    T = collect(0.0:0.005:1.0) .* Tend
    u1(t) = 1/(sqrt(σt)*π^(1/4)) * exp(-(t - t0)^2 / (2*σt^2))
    gu_t = coupling_input(u1, T)

    # -------- Example 05-1 style (symbolic -> numeric) --------
    ha(i) = NLevelSpace(Symbol("a$(i)"), 2)
    h = tensor([ha(i) for i = 1:N]...)
    σ(α, i, j) = Transition(h, Symbol("σ_$(α)"), i, j, α)

    γR(i) = real_var("γ^{($(i))}_R")
    γL(i) = real_var("γ^{($(i))}_L")
    Δ(i) = real_var("Δ_{$(i)}")
    ϕ(i, j) = real_var("ϕ_{$(i)$(j)}")
    Ein = real_var("E_{in}")

    G_d = SLH(1, Ein, 0)
    G_ϕ(i, j) = SLH(exp(1im * ϕ(i, j)), 0, 0)
    G_R(i) = SLH(1, √(γR(i)) * σ(i, 1, 2), -Δ(i) * σ(i, 2, 2))
    G_L(i) = SLH(1, √(γL(i)) * σ(i, 1, 2), 0)

    G_R_t = G_d ▷ G_R(1) ▷ G_ϕ(1, 2) ▷ G_R(2)
    G_L_t = G_L(2) ▷ G_ϕ(1, 2) ▷ G_L(1)
    G_t = G_R_t ⊞ G_L_t

    H = hamiltonian(G_t)
    L = lindblad(G_t)
    L_R = L[1]
    L_L = L[2]

    γ_ = 1.0
    β = 0.9
    γRn = fill(γ_ * β / 2, N)
    γLn = fill(γ_ * β / 2, N)
    γ_add = fill(γ_ * (1 - β), N)
    Δn = fill(0.0, N)
    ϕn = [0.0]

    p_sym = [γR(1); γR(2); γL(1); γL(2); Δ(1); Δ(2); ϕ(1, 2)]
    p_num = [γRn[1]; γRn[2]; γLn[1]; γLn[2]; Δn[1]; Δn[2]; ϕn[1]]
    dict_p = Dict(p_sym .=> p_num)
    dict_p_t = Dict(Ein => (t -> α0 * u1(t)))

    ba = NLevelBasis(2)
    b = tensor([ba for _ = 1:N]...)
    H_QO = translate_qo(H, b; parameter = dict_p, time_parameter = dict_p_t)
    L_R_QO = translate_qo(L_R, b; parameter = dict_p, time_parameter = dict_p_t)
    L_L_QO = translate_qo(L_L, b; parameter = dict_p, time_parameter = dict_p_t)

    σ_qo(α, i, j) = translate_qo(σ(α, i, j), b)
    J_add = [√(γ_add[i]) * σ_qo(i, 1, 2) for i = 1:N]

    function input_output_1(t, ρ)
        Ht = H_QO(t)
        J = [L_R_QO(t), L_L_QO(t), J_add...]
        return Ht, J, dagger.(J)
    end

    ψ0 = tensor([nlevelstate(ba, 1) for _ = 1:N]...)
    t1, ρt1 = timeevolution.master_dynamic(T, ψ0, input_output_1)

    I_R_1 = [real(expect(L_R_QO(ti)' * L_R_QO(ti), ρt1[i])) for (i, ti) in enumerate(t1)]
    I_L_1 = [real(expect(L_L_QO(ti)' * L_L_QO(ti), ρt1[i])) for (i, ti) in enumerate(t1)]

    # -------- Example 05-2 style (numeric SLH, quantum pulse) --------
    bu = FockBasis(4)
    bq = tensor([ba for _ = 1:N]...)
    b2 = bu ⊗ bq

    σ_qds(α, i, j) = embed(b2, 1 + α, transition(ba, i, j))
    a_u = destroy(bu) ⊗ one(bq)

    # Now uses unified SLH instead of SLHqo
    G_u = SLH(1, t -> gu_t(t) * a_u, 0 * one(b2))
    G_ϕ_qo(i, j) = SLH(exp(1im * ϕn[i]), 0 * one(b2), 0 * one(b2))
    G_R_qo(i) = SLH(1, √(γRn[i]) * σ_qds(i, 1, 2), -Δn[i] * σ_qds(i, 2, 2))
    G_L_qo(i) = SLH(1, √(γLn[i]) * σ_qds(i, 1, 2), 0 * one(b2))

    G_R_t_qo = G_u ▷ G_R_qo(1) ▷ G_ϕ_qo(1, 2) ▷ G_R_qo(2)
    G_L_t_qo = G_L_qo(2) ▷ G_ϕ_qo(1, 2) ▷ G_L_qo(1)
    G_t_qo = G_R_t_qo ⊞ G_L_t_qo

    H_qo = hamiltonian(G_t_qo)
    L_qo = lindblad(G_t_qo)
    L_R_qo = L_qo[1]
    L_L_qo = L_qo[2]

    # FunctionWrapper is callable but not <: Function
    _callable(x) = x isa Union{Function,FunctionWrapper}
    Hf = _callable(H_qo) ? H_qo : (t -> H_qo)
    L_R_f = _callable(L_R_qo) ? L_R_qo : (t -> L_R_qo)
    L_L_f = _callable(L_L_qo) ? L_L_qo : (t -> L_L_qo)

    J_add_qo = [√(γ_add[i]) * σ_qds(i, 1, 2) for i = 1:N]

    function input_output_2(t, ρ)
        Ht = Hf(t)
        J = [L_R_f(t), L_L_f(t), J_add_qo...]
        return Ht, J, dagger.(J)
    end

    ψ0_qo = coherentstate(bu, α0) ⊗ tensor([nlevelstate(ba, 1) for _ = 1:N]...)
    t2, ρt2 = timeevolution.master_dynamic(T, ψ0_qo, input_output_2)

    I_R_2 = [real(expect(L_R_f(ti)' * L_R_f(ti), ρt2[i])) for (i, ti) in enumerate(t2)]
    I_L_2 = [real(expect(L_L_f(ti)' * L_L_f(ti), ρt2[i])) for (i, ti) in enumerate(t2)]

    @test maximum(abs.(I_R_1 .- I_R_2)) < 5e-5
    @test maximum(abs.(I_L_1 .- I_L_2)) < 5e-5
    @test maximum(abs.(I_R_1)) > 1e-3
    @test maximum(abs.(I_L_1)) > 1e-3
end
