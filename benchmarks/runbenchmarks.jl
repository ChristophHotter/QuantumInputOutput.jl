using BenchmarkTools
using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using QuantumOpticsBase
using SymbolicUtils
using LinearAlgebra

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
results = BenchmarkTools.run(SUITE; verbose=true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
