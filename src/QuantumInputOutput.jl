module QuantumInputOutput

using SecondQuantizedAlgebra: SecondQuantizedAlgebra, QSym, to_numeric
using QuantumOpticsBase: QuantumOpticsBase, expect, basis, dagger, sparse
using QuantumOptics: QuantumOptics, timeevolution
using SymbolicUtils: SymbolicUtils, substitute, BasicSymbolic, simplify, expand
using Symbolics: Symbolics
using SpecialFunctions: erf
using DataInterpolations: LinearInterpolation, ExtrapolationType
using NumericalIntegration: cumul_integrate
using LinearAlgebra: LinearAlgebra, I, mul!
using OrdinaryDiffEq: OrdinaryDiffEq, ODEProblem, Tsit5, solve
using StaticArrays: StaticArrays, SMatrix, SVector
using FunctionWrappers: FunctionWrappers, FunctionWrapper

const SQA = SecondQuantizedAlgebra

# TODO: move to SQA
complex_var(name::AbstractString) = Symbolics.variable(String(name))
real_var(name::AbstractString) = Symbolics.variable(String(name); T = Real)

function complex_vars(names::AbstractString)
    parts = filter(!isempty, strip.(split(replace(names, ',' => ' '))))
    return Tuple(complex_var(part) for part in parts)
end

function real_vars(names::AbstractString)
    parts = filter(!isempty, strip.(split(replace(names, ',' => ' '))))
    return Tuple(real_var(part) for part in parts)
end

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
    translate_qo,
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
    correlation_matrix,
    # Operators
    substitute_operators,
    dagger,
    # Scalar symbol helpers
    complex_var,
    real_var,
    complex_vars,
    real_vars

include("SLH.jl")
include("translate.jl")
include("utils.jl")
include("pulses.jl")
include("correlations.jl")
include("interaction_picture.jl")

end
