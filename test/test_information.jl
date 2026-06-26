using Test
using LinearAlgebra
using QuantumInputOutput
using QuantumOptics
using SecondQuantizedAlgebra

@testset "quantum_fisher_information" begin
    ρ = Diagonal([0.3, 0.7])
    dρ = Diagonal([1.0, -1.0])
    @test quantum_fisher_information(Matrix(ρ), Matrix(dρ); regularization = 0.0) ≈
          1 / 0.3 + 1 / 0.7

    θ = 0.37
    c, s = cos(θ), sin(θ)
    ρ_pure = [c^2 c*s; c*s s^2]
    dρ_pure = [-2c*s c^2-s^2; c^2-s^2 2c*s]
    @test quantum_fisher_information(ρ_pure, dρ_pure) ≈ 4 atol = 1e-8

    b = NLevelBasis(2)
    ρ_op = Operator(b, b, Matrix(ρ))
    dρ_op = Operator(b, b, Matrix(dρ))
    @test quantum_fisher_information(ρ_op, dρ_op; regularization = 0.0) ≈
          quantum_fisher_information(Matrix(ρ), Matrix(dρ); regularization = 0.0)

    @test_throws ArgumentError quantum_fisher_information(ones(2, 3), ones(2, 3))
    @test_throws ArgumentError quantum_fisher_information(ones(2, 2), ones(3, 3))
end

@testset "classical_fisher_information" begin
    ρ = [0.25 0; 0 0.75]
    dρ = [0.5 0; 0 -0.5]
    M = projective_measurement([0.0 0; 0 1.0])

    @test povm_probabilities(ρ, M) ≈ [0.25, 0.75]
    @test classical_fisher_information(ρ, dρ, M) ≈ 0.5^2 / 0.25 + 0.5^2 / 0.75

    σx = [0 1; 1 0]
    P = projective_measurement(σx)
    @test length(P) == 2
    @test sum(P) ≈ I(2)
    @test all(Pi -> Pi * Pi ≈ Pi, P)

    @test classical_fisher_information([1.0 0; 0 0], [0.0 0; 0 1.0], M) == Inf
    @test_throws DomainError classical_fisher_information([-0.1 0; 0 1.1], dρ, M)
end

@testset "SLH parameter_derivative" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    @variables Δ::Real

    b = FockBasis(1)
    ψ0 = (fockstate(b, 0) + fockstate(b, 1)) / sqrt(2)
    ρ0 = dm(ψ0)
    G = SLH(1, 0 * a, Δ * a' * a)

    t = 0.7
    θ = 0.2
    ρ, dρ = parameter_derivative(G, b, ρ0, t; estimate = Δ, parameter = Dict(Δ => θ))

    phase = exp(1im * θ * t)
    expected_ρ = [0.5 0.5 * phase; 0.5 * conj(phase) 0.5]
    expected_dρ = [0 0.5im * t * phase; -0.5im * t * conj(phase) 0]
    @test Matrix(ρ.data) ≈ expected_ρ atol = 1e-8
    @test Matrix(dρ.data) ≈ expected_dρ atol = 1e-8
    @test quantum_fisher_information(ρ, dρ) ≈ t^2 atol = 1e-6
    @test classical_fisher_information(ρ, dρ, projective_measurement([0.0 0; 0 1.0])) ≈ 0 atol = 1e-10

    @test_throws ArgumentError parameter_derivative(G, b, ρ0, t; estimate = Δ, parameter = Dict())
end
