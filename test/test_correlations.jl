using QuantumInputOutput
using QuantumOptics
using QuantumOpticsBase: dagger, TimeDependentSum
using LinearAlgebra
using Test

@testset "correlations_single_tls_decay" begin
    γ = 1.0
    ba = NLevelBasis(2)
    σm = transition(ba, 1, 2)
    ψ0 = nlevelstate(ba, 2) # initially excited

    H = 0.0 * one(ba)
    J = [sqrt(γ) * σm]
    Js = J[1]

    T_end = 10 / γ
    T = collect(0.0:0.01:T_end)
    ΔT = T[2] - T[1]

    function f_const(t, ρ)
        return H, J, dagger.(J)
    end

    _, ρt_dyn = timeevolution.master_dynamic(T, ψ0, f_const)
    g1_dyn = correlation_matrix(T, ρt_dyn, f_const, Js)

    _, ρt_static = timeevolution.master(T, ψ0, H, J)
    g1_static = correlation_matrix(T, ρt_static, H, J, Js)

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

    @testset "time-dependent Js rejected" begin
        td_Js = TimeDependentSum([t -> 1.0 + 0im], J)
        @test_throws ArgumentError correlation_matrix(T, ρt_static, H, J, td_Js)
    end

    @testset "time-dependent operator path matches static" begin
        H_td = TimeDependentSum([t -> 0.0 + 0im], [one(ba)])
        J_td = [TimeDependentSum([t -> sqrt(γ) + 0im], [σm])]
        _, ρt_td = timeevolution.master_dynamic(T, ψ0, H_td, J_td)
        g1_td = correlation_matrix(T, ρt_td, H_td, J_td, Js)
        @test maximum(abs.(g1_td .- g1_static)) < 1e-8
    end
end
