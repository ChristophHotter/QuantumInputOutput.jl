module QuantumInputOutput

using SecondQuantizedAlgebra
using QuantumCumulants
using QuantumOpticsBase
# using QuantumOptics
using SymbolicUtils: SymbolicUtils, substitute, BasicSymbolic, operation, arguments, iscall, simplify
# using Symbolics #?
# using ModelingToolkit #?
using SpecialFunctions: erf
using DataInterpolations: LinearInterpolation, ExtrapolationType
using NumericalIntegration: cumul_integrate
using LinearAlgebra: LinearAlgebra, I
using OrdinaryDiffEq # for ui_to_u_i_im1

# import QuantumOpticsBase: expect
import QuantumCumulants: numeric_average, QNumber
const QC = QuantumCumulants

export SLH, # SLH.jl
    get_scattering,
    get_lindblad,
    get_hamiltonian,
    ▷,
    cascade, 
    ⊞,
    concatenation, 
    translate, # translate.jl
    u_to_gu, # utils.jl
    v_to_gv, 
    u_to_gu_Gauss,
    v_to_gv_Gauss,
    ui_to_u_i_im1,
    vi_to_v_i_im1,
    two_time_corr_matrix

include("SLH.jl")
include("translate.jl")
include("utils.jl")

end