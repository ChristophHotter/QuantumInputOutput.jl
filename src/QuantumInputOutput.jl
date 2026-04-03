module QuantumInputOutput

using SecondQuantizedAlgebra
using QuantumOpticsBase
using QuantumOptics
using SymbolicUtils:
    SymbolicUtils, substitute, BasicSymbolic, operation, arguments, iscall, simplify, expand
using Symbolics
using SpecialFunctions: erf
using DataInterpolations: LinearInterpolation, ExtrapolationType
using NumericalIntegration: cumul_integrate
using LinearAlgebra: LinearAlgebra, I, mul!
using OrdinaryDiffEq
using StaticArrays
using FunctionWrappers

import QuantumOpticsBase: expect
import SecondQuantizedAlgebra: numeric_average
const SQA = SecondQuantizedAlgebra

export SLH,
    Gaussian,
    scattering,
    lindblad,
    hamiltonian,
    # Backward compat
    get_scattering,
    get_lindblad,
    get_hamiltonian,
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
    substitute_operators

include("SLH.jl")
include("translate.jl")
include("utils.jl")
include("pulses.jl")
include("correlations.jl")
include("interaction_picture.jl")

end
