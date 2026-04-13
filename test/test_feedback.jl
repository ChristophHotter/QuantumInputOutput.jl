using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using SymbolicUtils
using FunctionWrappers: FunctionWrapper
using LinearAlgebra
using StaticArrays
using Test

@testset "feedback reduction" begin
    hs = FockSpace(:s)
    a = Destroy(hs, :a, 1)

    # # TODO: problem with simplify of conj(conj(x)), conj((0+1im)*x) and fractions
    # s11_r = rnumber("s11_r")
    # s12_r = rnumber("s12_r")
    # s21_r = rnumber("s21_r")
    # s22_r = rnumber("s22_r")
    # l1_r = rnumber("l1_r")
    # l2_r = rnumber("l2_r")
    # #
    # s11_i = rnumber("s11_i")
    # s12_i = rnumber("s12_i")
    # s21_i = rnumber("s21_i")
    # s22_i = rnumber("s22_i")
    # l1_i = rnumber("l1_i")
    # l2_i = rnumber("l2_i")
    # #
    # s11 = s11_r + 1im*s11_i
    # s12 = s12_r + 1im*s12_i
    # s21 = s21_r + 1im*s21_i
    # s22 = s22_r + 1im*s22_i
    # l1 = l1_r + 1im*l1_i
    # l2 = l2_r + 1im*l2_i

    # h0 = rnumber("h0")

    s11 = 0.1 + 1im*0.2
    s12 = 0.3 + 1im*0.24
    s21 = 0.5 + 1im*0.23
    s22 = 0.12 + 1im*0.29
    l1 = 0.18 + 1im*0.21
    l2 = 0.15 + 1im*0.25
    h0 = 1.0

    G = SLH([s11 s12; s21 s22], [l1, l2], h0)
    G_red = feedback(G, 1, 1)

    loop_gain = (1 - s11)^(-1)
    expected_S = simplify(s22 + s21 * loop_gain * s12)
    expected_L = simplify(l2 + s21 * loop_gain * l1)
    expected_term = simplify((l1' * s11 + l2' * s21) * loop_gain * l1)
    expected_H = simplify(h0 + (expected_term - expected_term') / (2im))

    @test scattering(G_red) isa SMatrix{1,1}
    @test abs(scattering(G_red)[1, 1] - expected_S) < 1e-10
    @test abs(lindblad(G_red)[1] - expected_L) < 1e-10
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
        @test iszero(simplify(lindblad(G_loop)[1] - simplify(l * √(κ) * a)))
        @test isequal(simplify(hamiltonian(G_loop) - hamiltonian(G_opo)), 0)
    end

    @testset "bidirectional waveguide matches cascade model" begin
        N = 2
        ha(i) = NLevelSpace("a$(i)", 2)
        h = tensor([ha(i) for i = 1:N]...)
        σ(α, i, j) = Transition(h, "σ_$(α)", i, j, α)

        γR(i) = rnumber("γ^{($(i))}_R")
        γL(i) = rnumber("γ^{($(i))}_L")
        Δqd(i) = rnumber("Δqd_{$(i)}")
        ϕ(i, j) = rnumber("ϕ_{$(i)$(j)}")
        Ein = rnumber("E_{in}")

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

        @test isequal(lindblad(G_feedback)[1], lindblad(G_manual)[2])
        @test isequal(lindblad(G_feedback)[2], lindblad(G_manual)[1])
        @test isequal(simplify(hamiltonian(G_feedback) - hamiltonian(G_manual)), 0)
    end

    @testset "feedback preserves FunctionWrapper" begin
        bc = FockBasis(4)
        a_op = destroy(bc)
        H_s = sparse(0.5 * dagger(a_op) * a_op)
        L_s = sparse(sqrt(1.0) * a_op)
        gu_f(t) = exp(-t^2) * sparse(a_op)
        gv_f(t) = exp(-(t - 2)^2) * sparse(a_op)

        G_cat = SLH(1, gu_f, H_s) ⊞ SLH(1, gv_f, H_s)
        G_fb = feedback(G_cat, 1, 1)

        LT = eltype(lindblad(G_fb))
        @test LT <: FunctionWrapper
        @test LT !== Any
        @inferred lindblad(G_fb)[1](0.5)
    end
end
