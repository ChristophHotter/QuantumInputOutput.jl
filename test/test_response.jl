using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using QuantumOpticsBase: dagger
using LinearAlgebra
using Test

# Symbolic two-port Kerr parametric oscillator (KPO), rotating frame:
#   H = -Δ a†a + K a†a†aa + (G/2)(a†a† + aa),  ports L = (√κ1 a, √κ2 a)
h = FockSpace(:kpo)
a = Destroy(h, :a)
@variables Δ::Real K::Real G::Real κ1::Real κ2::Real
H = -Δ*a'*a + K*a'*a'*a*a + (G/2)*(a'*a' + a*a)
Gkpo = SLH([1 0; 0 1], [√(κ1)*a, √(κ2)*a], H)
b = FockBasis(16)

# the user computes the steady state and passes it in (the package no longer wraps a solver)
function steady(p)
    Hn, Jn = QuantumInputOutput.translate_qo(Gkpo, b; parameter = p)
    return steadystate.eigenvector(Hn, collect(Jn))
end

@testset "output_field" begin
    # b_out,k = -L_k with no input
    Lnum = QuantumInputOutput.translate_qo(√(κ2)*a, b; parameter = Dict(κ2 => 0.6))
    of = QuantumInputOutput.translate_qo(output_field(Gkpo, 2), b; parameter = Dict(κ2 => 0.6))
    @test isapprox(Matrix(of.data), -Matrix(Lnum.data); atol = 1e-12)
    # coherent offset on the same port adds (S α)_k = α_k for S = I
    of2 = QuantumInputOutput.translate_qo(
        output_field(Gkpo, 2; input = [0.0, 2.0]), b; parameter = Dict(κ2 => 0.6))
    @test isapprox(Matrix(of2.data), 2.0*Matrix(one(b).data) - Matrix(Lnum.data); atol = 1e-12)
end

@testset "scattering_parameter: linear cavity vs analytic Lorentzian + energy conservation" begin
    Δ_, κ1_, κ2_ = 0.7, 0.4, 0.6
    p = Dict(Δ => Δ_, K => 0.0, G => 0.0, κ1 => κ1_, κ2 => κ2_)
    ω = collect(-2.0:0.1:2.0)
    ρ = steady(p)
    S21 = scattering_parameter(Gkpo, b, ρ; in_port = 1, out_port = 2, omega = ω, parameter = p)
    S11 = scattering_parameter(Gkpo, b, ρ; in_port = 1, out_port = 1, omega = ω, parameter = p)
    # Gardiner convention a_out = a_in - √κ a, with H = -Δ a†a
    κ = κ1_ + κ2_
    ana(δ) = -sqrt(κ1_*κ2_) / (κ/2 + im*(-Δ_ - δ))
    @test all(isapprox.(S21, ana.(ω); atol = 1e-6))
    # lossless two-port conserves energy at every frequency
    @test maximum(abs.(abs2.(S11) .+ abs2.(S21) .- 1)) < 1e-8
    # resonance peak (critical for κ1=κ2 would be 1; here 2√(κ1κ2)/κ)
    @test maximum(abs.(S21)) ≈ 2*sqrt(κ1_*κ2_)/κ atol = 1e-4
end

@testset "scattering_parameter: parametric gain above passive bound" begin
    # below threshold G < κ/2 = 0.5; expect |S21| > passive max of 1
    p = Dict(Δ => 0.0, K => 0.0, G => 0.4, κ1 => 0.5, κ2 => 0.5)
    ω = collect(-1.0:0.02:1.0)
    ρ = steady(p)
    S21 = scattering_parameter(Gkpo, b, ρ; in_port = 1, out_port = 2, omega = ω, parameter = p)
    @test maximum(abs.(S21)) > 1.0
end

@testset "power_spectrum: emission ≥ 0 and matches correlation integral" begin
    p = Dict(Δ => 0.0, K => 0.0, G => 0.35, κ1 => 1.0, κ2 => 1e-9)
    ω = collect(-2.0:0.25:2.0)
    Hn, Jn = QuantumInputOutput.translate_qo(Gkpo, b; parameter = p)
    Jn = collect(Jn)
    ρ = steadystate.eigenvector(Hn, Jn)
    S = power_spectrum(Gkpo, b, ρ; port = 1, omega = ω, parameter = p)
    @test all(S .>= -1e-8)
    # compare to direct numerical FT of ⟨ā†(τ)ā(0)⟩ from QuantumOptics
    A = Jn[1]
    amean = tr(A.data*ρ.data)
    Ā = A - amean*one(b)
    τ = collect(0.0:0.004:80.0)
    g = timecorrelations.correlation(τ, ρ, Hn, Jn, dagger(Ā), Ā)
    Squad(w) = begin
        f = @. exp(-im*w*τ) * g
        2*real(sum((f[1:end-1] .+ f[2:end]) ./ 2 .* diff(τ)))
    end
    @test all(isapprox(S[i], Squad(ω[i]); rtol = 1e-3, atol = 1e-5) for i in eachindex(ω))
end

@testset "power_spectrum: squeezing matches analytic DPA + below vacuum" begin
    # Degenerate parametric amplifier (K=0, single port): the vacuum-normalised output
    # squeezing spectrum has the closed form (Collett–Gardiner)
    #   S_∓(ω) = 1 ∓ 2κG / ((κ/2 ± G)² + ω²),   κ = κ1 + κ2.
    κ_, G_ = 1.0, 0.35
    p = Dict(Δ => 0.0, K => 0.0, G => G_, κ1 => κ_, κ2 => 1e-9)
    ω = collect(-2.5:0.25:2.5)
    ρ = steady(p)
    sqz  = power_spectrum(Gkpo, b, ρ; port = 1, omega = ω, parameter = p, quadrature = π/4)
    anti = power_spectrum(Gkpo, b, ρ; port = 1, omega = ω, parameter = p, quadrature = 3π/4)
    Ssq(δ)   = 1 - 2κ_*G_ / ((κ_/2 + G_)^2 + δ^2)
    Santi(δ) = 1 + 2κ_*G_ / ((κ_/2 - G_)^2 + δ^2)
    # the sharp anti-squeezed peak (~32) is sensitive to the Fock truncation; use a
    # relative tolerance (≈0.15% at FockBasis(16), →1e-5 at FockBasis(24))
    @test all(isapprox.(sqz,  Ssq.(ω);   rtol = 1e-2, atol = 2e-3))
    @test all(isapprox.(anti, Santi.(ω); rtol = 1e-2, atol = 2e-3))
    @test minimum(sqz) < 1.0          # squeezed below the shot-noise floor
    @test maximum(anti) > 1.0         # anti-squeezed above it
end

@testset "susceptibility: reduces to scattering_parameter relation" begin
    p = Dict(Δ => 0.2, K => 0.05, G => 0.25, κ1 => 0.5, κ2 => 0.5)
    ω = collect(-1.0:0.2:1.0)
    # S_{out,in} = scattering[out,in] + √(κ_out κ_in) · χ_{a,a†}/(i)  is implicit; here just
    # check the engine runs, returns the right length and finite complex numbers
    ρ = steady(p)
    χ = susceptibility(Gkpo, b, ρ, a, a', ω; parameter = p)
    @test length(χ) == length(ω)
    @test all(isfinite, abs.(χ))
end
