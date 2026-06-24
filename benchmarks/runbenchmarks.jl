using BenchmarkTools
using QuantumInputOutput
using SecondQuantizedAlgebra
using Symbolics: Symbolics
using QuantumOptics
using QuantumOpticsBase
using SymbolicUtils
using LinearAlgebra
using FunctionWrappers

const SUITE = BenchmarkGroup()

include("slh_algebra.jl")
include("translation.jl")
include("pulse_couplings.jl")
include("interaction_picture.jl")
include("correlations.jl")

benchmark_slh_algebra!(SUITE)
benchmark_translation!(SUITE)
benchmark_pulse_couplings!(SUITE)
benchmark_interaction_picture!(SUITE)
benchmark_correlations!(SUITE)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose = true)

# Report the minimum rather than the median. The minimum is BenchmarkTools'
# recommended estimator for tracking: measurement noise (GC pauses, scheduler
# preemption, frequency scaling) is strictly additive, so the minimum is the
# most reproducible estimate of the underlying cost and the least sensitive to
# cross-runner variance. Allocation-heavy benchmarks additionally set
# `gcsample=true` so each sample starts from a clean heap (see the individual
# `@benchmarkable`s).
display(minimum(results))

BenchmarkTools.save("benchmarks_output.json", minimum(results))
