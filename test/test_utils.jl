using QuantumInputOutput
using QuantumInputOutput: dagger
using Test

@testset "utils" begin
    τ = 3.0
    σ = 0.5
    δ = 0#0.3

    u(t) = 1 / (sqrt(σ) * π^(1/4)) * exp(-0.5 * ((t - τ) / σ)^2) * exp(1im * δ * t)
    v(t) = u(t)

    T = [0.0:0.0001:1;]*6

    gu_num = coupling_input(u, T)
    gv_num = coupling_output(v, T)

    gu_num2 = coupling_input(u.(T), T)
    gv_num2 = coupling_output(v.(T), T)

    @test sum(abs.(gu_num2.(T) - gu_num.(T))) < 1e-9
    @test sum(abs.(gv_num2.(T) - gv_num.(T))) < 1e-9

    gu_ana = coupling_input(Gaussian(τ, σ; δ = δ))
    gv_ana = coupling_output(Gaussian(τ, σ; δ = δ))

    gv_err = maximum(abs.(gv_num.(T[2:end]) .- gv_ana.(T[2:end])))
    gu_err = maximum(abs.(gu_num.(T[2:end]) .- gu_ana.(T[2:end])))

    @test gv_err < 5e-4
    @test gu_err < 5e-4

    v1(t) = 1 / (sqrt(σ) * π^(1/4)) * exp(-0.5 * ((t - (τ - 2σ)) / σ)^2)
    v2(t) = 1 / (sqrt(σ) * π^(1/4)) * exp(-0.5 * ((t - (τ + 2σ)) / σ)^2)
    u1(t) = 1 / (sqrt(σ) * π^(1/4)) * exp(-0.5 * ((t - (τ - 2σ)) / σ)^2)
    u2(t) = 1 / (sqrt(σ) * π^(1/4)) * exp(-0.5 * ((t - (τ + 2σ)) / σ)^2)

    v_eff_f = effective_output_mode([v1, v2], T, 2)
    v_eff_data = effective_output_mode([v1.(T), v2.(T)], T, 2)
    @test maximum(abs.(v_eff_f.(T) .- v_eff_data.(T))) < 1e-7

    gv1 = coupling_output(v1, T)
    gv2 = coupling_output(v2, T)
    v_eff_data2 = effective_output_mode([v1.(T), v2.(T)], [gv1.(T), gv2.(T)], T, 2)
    @test maximum(abs.(v_eff_f.(T) .- v_eff_data2.(T))) < 1e-7

    u_eff_f = effective_input_mode([u1, u2], T, 2)
    u_eff_data = effective_input_mode([u1.(T), u2.(T)], T, 2)
    @test_broken maximum(abs.(u_eff_f.(T) .- u_eff_data.(T))) /
          maximum(abs.(u_eff_f.(T) .+ u_eff_data.(T))) < 1e-5

    gu1 = coupling_input(u1, T)
    gu2 = coupling_input(u2, T)
    u_eff_data2 = effective_input_mode([u1.(T), u2.(T)], [gu1.(T), gu2.(T)], T, 2)
    @test_broken maximum(abs.(u_eff_f.(T) .- u_eff_data2.(T))) /
          maximum(abs.(u_eff_f.(T) .+ u_eff_data2.(T))) < 1e-5
end
