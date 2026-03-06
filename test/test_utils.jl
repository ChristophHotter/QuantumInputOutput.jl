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
end
