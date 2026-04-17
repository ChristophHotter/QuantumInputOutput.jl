using QuantumInputOutput
using QuantumOptics
using LinearAlgebra
using Test

@testset "correlations_single_tls_decay" begin
    γ = 1.0
    ba = NLevelBasis(2)
    σm = transition(ba, 1, 2)
    ψ0 = nlevelstate(ba, 2) # initially excited

    H = 0.0 * one(ba)
    J = [sqrt(γ) * σm]
    Ls = J[1]

    T_end = 10 / γ
    T = collect(0.0:0.01:T_end)
    ΔT = T[2] - T[1]

    function f_const(t, ρ)
        return H, J, dagger.(J)
    end

    _, ρt_dyn = timeevolution.master_dynamic(T, ψ0, f_const)
    g1_dyn = correlation_matrix(T, ρt_dyn, f_const, Ls)

    _, ρt_static = timeevolution.master(T, ψ0, H, J)
    g1_static = correlation_matrix(T, ρt_static, H, J, Ls)

    @test maximum(abs.(g1_dyn .- g1_static)) < 1e-8

    F = eigen(Hermitian((g1_static + g1_static') / 2))
    n_modes = real.(F.values) * ΔT

    expected_n = 1 - exp(-γ * T_end)
    @test abs(n_modes[end] - expected_n) < 1e-2
    @test abs(n_modes[end-1]) < 1e-2

    mode_num = F.vectors[:, end] / sqrt(ΔT)
    mode_exp = sqrt.(γ) .* exp.(-0.5 * γ .* T)
    mode_exp ./= sqrt(sum(abs2.(mode_exp)) * ΔT)

    overlap = sum(conj.(mode_exp) .* mode_num) * ΔT
    mode_num_aligned = mode_num * exp(-1im * angle(overlap))
    overlap_aligned = abs(sum(conj.(mode_exp) .* mode_num_aligned) * ΔT)
    @test overlap_aligned > 0.99
end
