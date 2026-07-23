using QuantumInputOutput
using DataInterpolations: LinearInterpolation
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

    # Every coupling constructor returns a `PulseCoupling`: callable and broadcastable.
    @test gu_num isa PulseCoupling
    @test gv_num isa PulseCoupling
    @test gu_ana isa PulseCoupling
    @test gv_ana isa PulseCoupling
    @test gu_num(T[5]) isa ComplexF64
    @test gu_num.(T[1:3]) isa AbstractVector{ComplexF64}

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
    @test maximum(abs.(u_eff_f.(T) .- u_eff_data.(T))) /
          maximum(abs.(u_eff_f.(T) .+ u_eff_data.(T))) < 1e-2

    gu1 = coupling_input(u1, T)
    gu2 = coupling_input(u2, T)
    u_eff_data2 = effective_input_mode([u1.(T), u2.(T)], [gu1.(T), gu2.(T)], T, 2)
    @test maximum(abs.(u_eff_f.(T) .- u_eff_data2.(T))) /
          maximum(abs.(u_eff_f.(T) .+ u_eff_data2.(T))) < 1e-2

    # zero coupling ⇒ no distortion ⇒ effective mode equals the bare mode (exercises i = 3)
    zero_c(t) = 0.0 + 0.0im
    v3(t) = 1 / (sqrt(σ) * π^(1/4)) * exp(-0.5 * ((t - τ) / σ)^2)
    v_eff3 = effective_output_mode([v1, v2, v3], [zero_c, zero_c, zero_c], T, 3)
    u_eff3 = effective_input_mode([u1, u2, v3], [zero_c, zero_c, zero_c], T, 3)
    @test maximum(abs.(v_eff3.(T) .- v3.(T))) < 1e-8
    @test maximum(abs.(u_eff3.(T) .- v3.(T))) < 1e-8

    # constant modes give the closed form g_out(tᵢ) = u / sqrt((|v|² - |u|²) i Δt)
    a_c = 0.5
    b_c = 0.8
    Td = collect(0.0:0.01:2.0)
    dt_c = Td[2] - Td[1]
    uc(t) = a_c
    vc(t) = b_c
    g_out = coupling_delay_out(uc, vc, Td)
    g_in = coupling_delay_in(uc, vc, Td)
    exp_out = [a_c / sqrt((b_c^2 - a_c^2) * (i * dt_c)) for i = 2:length(Td)]
    exp_in = [-b_c / sqrt((b_c^2 - a_c^2) * (i * dt_c)) for i = 2:length(Td)]
    @test g_out isa PulseCoupling
    @test maximum(abs.(g_out.(Td[2:end]) .- exp_out)) < 1e-6
    @test maximum(abs.(g_in.(Td[2:end]) .- exp_in)) < 1e-6
    # the shared denominator cancels, so in/out = -v/u exactly
    @test maximum(abs.(g_in.(Td[2:end]) ./ g_out.(Td[2:end]) .- (-b_c / a_c))) < 1e-12

    # LinearInterpolation dispatch must reproduce the analytic references
    u_itp = LinearInterpolation(u.(T), T)
    @test maximum(abs.(coupling_input(u_itp, T).(T[2:end]) .- gu_ana.(T[2:end]))) < 5e-4
    @test maximum(abs.(coupling_output(u_itp, T).(T[2:end]) .- gv_ana.(T[2:end]))) < 5e-4
    uc_itp = LinearInterpolation(uc.(Td), Td)
    vc_itp = LinearInterpolation(vc.(Td), Td)
    @test maximum(abs.(coupling_delay_out(uc_itp, vc_itp, Td).(Td[2:end]) .- exp_out)) < 1e-6
    @test maximum(abs.(coupling_delay_in(uc_itp, vc_itp, Td).(Td[2:end]) .- exp_in)) < 1e-6
end
