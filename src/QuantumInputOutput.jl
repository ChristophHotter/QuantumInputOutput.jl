module QuantumInputOutput

using SecondQuantizedAlgebra
using QuantumCumulants
using QuantumOpticsBase
using QuantumOptics
using SymbolicUtils:
    SymbolicUtils, substitute, BasicSymbolic, operation, arguments, iscall, simplify, expand
using Symbolics #?
# using ModelingToolkit #?
using SpecialFunctions: erf
using DataInterpolations: LinearInterpolation, ExtrapolationType
using NumericalIntegration: cumul_integrate
using LinearAlgebra: LinearAlgebra, I
using OrdinaryDiffEq # for u_eff

import QuantumOpticsBase: expect
import SecondQuantizedAlgebra: numeric_average
const SQA = SecondQuantizedAlgebra

export SLH, # SLH.jl
    SLHqo,
    get_scattering,
    get_lindblad,
    get_hamiltonian,
    ▷,
    cascade,
    ⊞,
    concatenate,
    translate, # translate.jl
    u_to_gu, # utils.jl
    v_to_gv,
    u_to_gu_Gauss,
    v_to_gv_Gauss,
    u_eff,
    v_eff,
    uv_to_gout,
    uv_to_gin,
    two_time_corr_matrix,
    interaction_picture_M,
    interaction_picture_M_2modes_equal,
    interaction_picture_A_2modes,
    interaction_picture_A_3modes,
    interaction_picture_A_4modes,
    interaction_picture_A_uv,
    interaction_picture_A_ucv,
    interaction_picture_A_uuvv,
    substitute_operators

include("SLH.jl")
include("translate.jl")
include("utils.jl")
include("pulses.jl")
include("correlations.jl")
include("interaction_picture.jl")

end
