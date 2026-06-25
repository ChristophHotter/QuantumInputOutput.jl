using BenchmarkTools
using QuantumInputOutput
using SecondQuantizedAlgebra
using Symbolics: Symbolics
using QuantumOptics
using QuantumOpticsBase
using SymbolicUtils
using LinearAlgebra
using FunctionWrappers

# Keep the benchmarks single-threaded for reproducibility. Julia's task and GC
# threads are pinned via `--threads=1 --gcthreads=1` in the workflow; OpenBLAS
# multithreads the dense linear algebra in the numeric ODE/correlation
# benchmarks by default, so pin it here too. Otherwise the runner's core count
# and scheduler add cross-run noise.
LinearAlgebra.BLAS.set_num_threads(1)

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
# cross-runner variance.
display(minimum(results))

BenchmarkTools.save("benchmarks_output.json", minimum(results))
