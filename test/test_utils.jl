using QuantumInputOutput
using SpecialFunctions: erf
using Test

@testset "utils" begin
    τ = 2.0
    σ = 0.5
    δ = 0.3

    u(t) = 1 / (sqrt(σ) * π^(1/4)) * exp(-0.5 * ((t - τ) / σ)^2) * exp(1im * δ * t)
    v(t) = u(t)

    T = [0.0:0.0001:1;]*6

    gu_num = u_to_gu(u, T)
    gv_num = v_to_gv(v, T)

    gu_ana = u_to_gu_Gauss(τ, σ; δ=δ)
    gv_ana = v_to_gv_Gauss(τ, σ; δ=δ)

    # sum(abs.(gu_num.(T) .- gu_ana.(T)))
    # sum(abs.(gu_num.(T) .- gu_ana.(T)))

    ∫_2(t) = 0.5 * (erf((t - τ) / σ) + erf(τ / σ))
    T_gv = [t for t in T if ∫_2(t) > 1e-6]
    T_gu = [t for t in T if (1 - ∫_2(t)) > 1e-6]

    gv_err = maximum(abs.(gv_num.(T_gv) .- gv_ana.(T_gv)))
    gu_err = maximum(abs.(gu_num.(T_gu) .- gu_ana.(T_gu)))

    @test gv_err < 5e-3
    @test gu_err < 5e-3
end
