using QuantumInputOutput
using Test

@testset "utils" begin
    τ = 3.0
    σ = 0.5
    δ = 0#0.3

    u(t) = 1 / (sqrt(σ) * π^(1/4)) * exp(-0.5 * ((t - τ) / σ)^2) * exp(1im * δ * t)
    v(t) = u(t)

    T = [0.0:0.0001:1;]*6

    gu_num = u_to_gu(u, T)
    gv_num = v_to_gv(v, T)

    gu_num2 = u_to_gu(u.(T), T)
    gv_num2 = v_to_gv(v.(T), T)

    @test sum(abs.(gu_num2.(T) - gu_num.(T))) < 1e-9
    @test sum(abs.(gv_num2.(T) - gv_num.(T))) < 1e-9

    gu_ana = u_to_gu_Gauss(τ, σ; δ = δ)
    gv_ana = v_to_gv_Gauss(τ, σ; δ = δ)

    gv_err = maximum(abs.(gv_num.(T[2:end]) .- gv_ana.(T[2:end])))
    gu_err = maximum(abs.(gu_num.(T[2:end]) .- gu_ana.(T[2:end])))

    @test gv_err < 5e-4
    @test gu_err < 5e-4

    v1(t) = 1 / (sqrt(σ) * π^(1/4)) * exp(-0.5 * ((t - (τ - 0.4)) / σ)^2)
    v2(t) = 1 / (sqrt(σ) * π^(1/4)) * exp(-0.5 * ((t - (τ + 0.4)) / σ)^2)
    u1(t) = 1 / (sqrt(σ) * π^(1/4)) * exp(-0.5 * ((t - (τ - 0.4)) / σ)^2)
    u2(t) = 1 / (sqrt(σ) * π^(1/4)) * exp(-0.5 * ((t - (τ + 0.4)) / σ)^2)

    v_eff_f = v_eff([v1, v2], T, 2)
    v_eff_data = v_eff([v1.(T), v2.(T)], T, 2)
    @test maximum(abs.(v_eff_f.(T) .- v_eff_data.(T))) < 1e-9

    gv1 = v_to_gv(v1, T)
    gv2 = v_to_gv(v2, T)
    v_eff_data2 = v_eff([v1.(T), v2.(T)], [gv1.(T), gv2.(T)], T, 2)
    @test maximum(abs.(v_eff_f.(T) .- v_eff_data2.(T))) < 1e-9

    u_eff_f = u_eff([u1, u2], T, 2)
    u_eff_data = u_eff([u1.(T), u2.(T)], T, 2)
    @test maximum(abs.(u_eff_f.(T) .- u_eff_data.(T))) < 1e-9

    gu1 = u_to_gu(u1, T)
    gu2 = u_to_gu(u2, T)
    u_eff_data2 = u_eff([u1.(T), u2.(T)], [gu1.(T), gu2.(T)], T, 2)
    @test maximum(abs.(u_eff_f.(T) .- u_eff_data2.(T))) < 1e-9
end
