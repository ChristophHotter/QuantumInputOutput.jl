using QuantumInputOutput
using SecondQuantizedAlgebra
using SymbolicUtils
using LinearAlgebra
using Test

@testset "feedback reduction" begin
    hs = FockSpace(:s)
    a = Destroy(hs, :a, 1)

    s11 = cnumber("s11")
    s12 = cnumber("s12")
    s21 = cnumber("s21")
    s22 = cnumber("s22")
    l1 = cnumber("l1")
    l2 = cnumber("l2")
    h0 = rnumber("h0")

    G = SLH([s11 s12; s21 s22], [l1, l2], h0)
    G_red = feedback(G, 1, 1)

    loop_gain = 1 / (1 - s11)
    expected_S = reshape([simplify(s22 + s21 * loop_gain * s12)], 1, 1)
    expected_L = [simplify(l2 + s21 * loop_gain * l1)]
    expected_term = simplify((l1' * s11 + l2' * s21) * loop_gain * l1)
    expected_H = simplify(h0 + (expected_term - expected_term') / (2im))

    @test isequal(G_red.scattering, expected_S)
    @test isequal(G_red.lindblad, expected_L)
    @test isequal(G_red.hamiltonian, expected_H)

    @testset "coherent-feedback OPO loop" begin
        κ = rnumber("κ")
        ϵ = rnumber("ϵ")
        η = rnumber("η")
        r = simplify(√(1 - η^2))

        G_opo = SLH(1, √(κ) * a, 1im * ϵ * (a'^2 - a^2))
        G_bs = SLH([-r η; η r], [0, 0], 0)
        G_unconnected = G_opo ⊞ G_bs
        G_loop = feedback(G_unconnected, 1 => 2, 2 => 1)

        l = simplify(η / (1 + r))

        @test isequal(G_loop.scattering, ones(1, 1))
        @test isequal(G_loop.lindblad, [simplify(l * √(κ) * a)])
        @test isequal(G_loop.hamiltonian, get_hamiltonian(G_opo))
    end

    @testset "bidirectional waveguide matches cascade model" begin
        N = 2
        ha(i) = NLevelSpace("a$(i)", 2)
        h = tensor([ha(i) for i = 1:N]...)
        σ(α, i, j) = Transition(h, "σ_$(α)", i, j, α)

        γR(i) = rnumber("γ^{($(i))}_R")
        γL(i) = rnumber("γ^{($(i))}_L")
        Δ(i) = rnumber("Δ_{$(i)}")
        ϕ(i, j) = rnumber("ϕ_{$(i)$(j)}")
        Ein = rnumber("E_{in}")

        G_d = SLH(1, Ein, 0)
        G_ϕ(i, j) = SLH(exp(1im * ϕ(i, j)), 0, 0)
        G_R(i) = SLH(1, √(γR(i)) * σ(i, 1, 2), -Δ(i) * σ(i, 2, 2))
        G_L(i) = SLH(1, √(γL(i)) * σ(i, 1, 2), 0)

        G_manual = (G_d ▷ G_R(1) ▷ G_ϕ(1, 2) ▷ G_R(2)) ⊞ (G_L(2) ▷ G_ϕ(1, 2) ▷ G_L(1))

        I2 = Matrix{Int}(I, 2, 2)
        G_in = SLH(I2, [Ein, 0], 0)
        G_qd(i) = SLH(I2, [√(γR(i)) * σ(i, 1, 2), √(γL(i)) * σ(i, 1, 2)], -Δ(i) * σ(i, 2, 2))
        G_phase(i, j) = SLH(
            [exp(1im * ϕ(i, j)) 0; 0 exp(1im * ϕ(i, j))],
            [0, 0],
            0,
        )

        G_network = G_in ⊞ G_qd(1) ⊞ G_phase(1, 2) ⊞ G_qd(2)
        G_feedback = feedback(G_network, 1 => 3, 3 => 5, 5 => 7, 8 => 6, 6 => 4, 4 => 2)

        @test isequal(G_feedback.scattering, get_scattering(G_manual))
        @test isequal(G_feedback.lindblad[1], get_lindblad(G_manual)[2])
        @test isequal(G_feedback.lindblad[2], get_lindblad(G_manual)[1])
        @test isequal(G_feedback.hamiltonian, get_hamiltonian(G_manual))
    end
end
