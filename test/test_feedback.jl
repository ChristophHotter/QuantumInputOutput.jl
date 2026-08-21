using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using SymbolicUtils
using Symbolics: Symbolics
using QuantumOpticsBase: dagger
using FunctionWrappers: FunctionWrapper
using LinearAlgebra
using StaticArrays
using Test

@testset "feedback reduction" begin
    hs = FockSpace(:s)
    a = Destroy(hs, :a)

    # # TODO: problem with simplify of conj(conj(x)), conj((0+1im)*x) and fractions
    # s11_r = rnumber("s11_r")
    # s12_r = rnumber("s12_r")
    # s21_r = rnumber("s21_r")
    # s22_r = rnumber("s22_r")
    # j1_r = rnumber("j1_r")
    # j2_r = rnumber("j2_r")
    # #
    # s11_i = rnumber("s11_i")
    # s12_i = rnumber("s12_i")
    # s21_i = rnumber("s21_i")
    # s22_i = rnumber("s22_i")
    # j1_i = rnumber("j1_i")
    # j2_i = rnumber("j2_i")
    # #
    # s11 = s11_r + 1im*s11_i
    # s12 = s12_r + 1im*s12_i
    # s21 = s21_r + 1im*s21_i
    # s22 = s22_r + 1im*s22_i
    # j1 = j1_r + 1im*j1_i
    # j2 = j2_r + 1im*j2_i

    # h0 = rnumber("h0")

    s11 = 0.1 + 1im*0.2
    s12 = 0.3 + 1im*0.24
    s21 = 0.5 + 1im*0.23
    s22 = 0.12 + 1im*0.29
    j1 = 0.18 + 1im*0.21
    j2 = 0.15 + 1im*0.25
    h0 = 1.0

    G = SLH([s11 s12; s21 s22], [j1, j2], h0)
    G_red = feedback(G, 1, 1)

    loop_gain = (1 - s11)^(-1)
    expected_S = simplify(s22 + s21 * loop_gain * s12)
    expected_J = simplify(j2 + s21 * loop_gain * j1)
    expected_term = simplify((j1' * s11 + j2' * s21) * loop_gain * j1)
    expected_H = simplify(h0 + (expected_term - expected_term') / (2im))

    @test scattering(G_red) isa SMatrix{1,1}
    @test abs(scattering(G_red)[1, 1] - expected_S) < 1e-10
    @test abs(jump_operator(G_red)[1] - expected_J) < 1e-10
    @test abs(hamiltonian(G_red) - expected_H) < 1e-10

    @testset "coherent-feedback OPO loop" begin
        # TODO: problem with simplify of conj(conj(x)), conj((0+1im)*x) and fractions
        # κ = rnumber("κ")
        # ϵ = rnumber("ϵ")
        # η = rnumber("η")
        κ = 0.7
        ϵ = 0.45
        η = 0.65

        r = simplify(√(1 - η^2))

        G_opo = SLH(1, √(κ) * a, 1im * ϵ * (a'^2 - a^2))
        G_bs = SLH([-r η; η r], [0, 0], 0)
        G_unconnected = G_opo ⊞ G_bs
        G_loop = feedback(G_unconnected, 1 => 2, 2 => 1)
        l = simplify(η / (1 + r))

        @test scattering(G_loop) isa SMatrix{1,1}
        @test iszero(simplify(jump_operator(G_loop)[1] - simplify(l * √(κ) * a)))
        @test iszero(simplify(hamiltonian(G_loop) - hamiltonian(G_opo)))
    end

    @testset "bidirectional waveguide matches cascade model" begin
        N = 2
        ha(i) = NLevelSpace(Symbol("a$(i)"), 2)
        h = tensor([ha(i) for i = 1:N]...)
        σ(α, i, j) = Transition(h, Symbol("σ_$(α)"), i, j, α)

        γR(i) = Symbolics.variable(Symbol("γ^{($(i))}_R"); T = Real)
        γL(i) = Symbolics.variable(Symbol("γ^{($(i))}_L"); T = Real)
        Δqd(i) = Symbolics.variable(Symbol("Δqd_{$(i)}"); T = Real)
        ϕ(i, j) = Symbolics.variable(Symbol("ϕ_{$(i)$(j)}"); T = Real)
        @variables Ein::Real

        G_d = SLH(1, Ein, 0)
        G_ϕ(i, j) = SLH(exp(1im * ϕ(i, j)), 0, 0)
        G_R(i) = SLH(1, √(γR(i)) * σ(i, 1, 2), -Δqd(i) * σ(i, 2, 2))
        G_L(i) = SLH(1, √(γL(i)) * σ(i, 1, 2), 0)

        G_manual = (G_d ▷ G_R(1) ▷ G_ϕ(1, 2) ▷ G_R(2)) ⊞ (G_L(2) ▷ G_ϕ(1, 2) ▷ G_L(1))

        I2 = Matrix{Int}(I, 2, 2)
        G_in = SLH(I2, [Ein, 0], 0)
        G_qd(i) =
            SLH(I2, [√(γR(i)) * σ(i, 1, 2), √(γL(i)) * σ(i, 1, 2)], -Δqd(i) * σ(i, 2, 2))
        G_phase(i, j) = SLH([exp(1im * ϕ(i, j)) 0; 0 exp(1im * ϕ(i, j))], [0, 0], 0)

        G_network = G_in ⊞ G_qd(1) ⊞ G_phase(1, 2) ⊞ G_qd(2)
        G_feedback = feedback(G_network, 1 => 3, 3 => 5, 5 => 7, 8 => 6, 6 => 4, 4 => 2)

        @test isequal(jump_operator(G_feedback)[1], jump_operator(G_manual)[2])
        @test isequal(jump_operator(G_feedback)[2], jump_operator(G_manual)[1])
        @test iszero(simplify(hamiltonian(G_feedback) - hamiltonian(G_manual)))
    end

    @testset "feedback preserves FunctionWrapper" begin
        bc = FockBasis(4)
        a_op = destroy(bc)
        H_s = sparse(0.5 * dagger(a_op) * a_op)
        J_s = sparse(sqrt(1.0) * a_op)
        gu_f(t) = exp(-t^2) * sparse(a_op)
        gv_f(t) = exp(-(t - 2)^2) * sparse(a_op)

        G_cat = SLH(1, gu_f, H_s) ⊞ SLH(1, gv_f, H_s)
        G_fb = feedback(G_cat, 1, 1)

        JT = eltype(jump_operator(G_fb))
        @test JT <: FunctionWrapper
        @test JT !== Any
        @inferred jump_operator(G_fb)[1](0.5)
    end
end
