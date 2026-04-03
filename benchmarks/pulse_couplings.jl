"""
Pulse coupling benchmarks — converting temporal mode shapes to virtual cavity couplings.
"""
function benchmark_pulse_couplings!(SUITE)
    SUITE["Pulse Couplings"] = BenchmarkGroup()

    ## Setup: Gaussian pulse on a realistic time grid
    τ = 3.0
    σ = 0.5
    δ = 0.3
    T = [0.0:0.0001:1;] * 6  # 60001 points

    u(t) = 1 / (sqrt(σ) * π^(1 / 4)) * exp(-0.5 * ((t - τ) / σ)^2) * exp(1im * δ * t)
    v(t) = u(t)

    ## --- Single-pulse coupling computation ---

    SUITE["Pulse Couplings"]["single-pulse"] = BenchmarkGroup()

    SUITE["Pulse Couplings"]["single-pulse"]["input"] =
        @benchmarkable u_to_gu($u, $T)

    SUITE["Pulse Couplings"]["single-pulse"]["output"] =
        @benchmarkable v_to_gv($v, $T)

    ## --- Effective multi-pulse couplings ---

    SUITE["Pulse Couplings"]["multi-pulse"] = BenchmarkGroup()

    v1(t) = 1 / (sqrt(σ) * π^(1 / 4)) * exp(-0.5 * ((t - (τ - 2σ)) / σ)^2)
    v2(t) = 1 / (sqrt(σ) * π^(1 / 4)) * exp(-0.5 * ((t - (τ + 2σ)) / σ)^2)
    u1(t) = 1 / (sqrt(σ) * π^(1 / 4)) * exp(-0.5 * ((t - (τ - 2σ)) / σ)^2)
    u2(t) = 1 / (sqrt(σ) * π^(1 / 4)) * exp(-0.5 * ((t - (τ + 2σ)) / σ)^2)

    SUITE["Pulse Couplings"]["multi-pulse"]["output 2 modes"] =
        @benchmarkable v_eff([$v1, $v2], $T, 2)

    SUITE["Pulse Couplings"]["multi-pulse"]["input 2 modes"] =
        @benchmarkable u_eff([$u1, $u2], $T, 2)

    return nothing
end
