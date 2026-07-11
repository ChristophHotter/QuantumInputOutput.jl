module QuantumInputOutput

using SecondQuantizedAlgebra: SecondQuantizedAlgebra, to_numeric
using QuantumOpticsBase: QuantumOpticsBase, expect, dagger
using QuantumOptics: QuantumOptics, timeevolution
using SymbolicUtils: SymbolicUtils, substitute, BasicSymbolic, simplify
using Symbolics: Symbolics
using SpecialFunctions: erf
using DataInterpolations: LinearInterpolation, ExtrapolationType
using NumericalIntegration: cumul_integrate
using LinearAlgebra: LinearAlgebra, I, mul!
using OrdinaryDiffEq: OrdinaryDiffEq, ODEProblem, Tsit5, solve
using StaticArrays: StaticArrays, SMatrix, SVector
using FunctionWrappers: FunctionWrappers, FunctionWrapper

const SQA = SecondQuantizedAlgebra

export SLH,
    Gaussian,
    scattering,
    lindblad,
    hamiltonian,
    # Composition
    ▷,
    cascade,
    ⊞,
    concatenate,
    feedback,
    # Translation
    to_numeric,
    substitute,
    # Pulse coupling
    coupling_input,
    coupling_output,
    coupling_delay_in,
    coupling_delay_out,
    # Effective modes
    effective_input_mode,
    effective_output_mode,
    # Interaction picture
    coupling_matrix,
    solve_mode_evolution,
    solve_mode_evolution_symmetric,
    # Correlations
    correlation_matrix

include("SLH.jl")
include("pulses.jl")
include("correlations.jl")
include("interaction_picture.jl")

end
