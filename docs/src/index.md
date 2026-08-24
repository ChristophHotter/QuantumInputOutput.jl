# Introduction

[![Documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://qojulia.github.io/QuantumInputOutput.jl/stable/)
[![Documentation dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://qojulia.github.io/QuantumInputOutput.jl/dev/)
[![Tests](https://github.com/qojulia/QuantumInputOutput.jl/actions/workflows/Tests.yml/badge.svg)](https://github.com/qojulia/QuantumInputOutput.jl/actions/workflows/Tests.yml)
[![codecov](https://codecov.io/gh/qojulia/QuantumInputOutput.jl/graph/badge.svg)](https://codecov.io/gh/qojulia/QuantumInputOutput.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![JET](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)

**QuantumInputOutput.jl** is a Julia framework for modeling the input-output theory with quantum pulses, using SLH (scattering, Lindblad, Hamiltonian) elements and rules. It combines:
- a symbolic layer based on [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl) to build models using SLH rules
- a numerical layer based on [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) and [QuantumCumulants.jl](https://github.com/qojulia/QuantumCumulants.jl) to simulate time dynamics and observables in a full quantum or higher-order meanfield approach, respectively

The typical workflow is:
1. Build the SLH model symbolically
2. Translate to numerical operators
3. Evolve the system in time
4. Analyze output modes and correlations

## Key Features

- SLH modeling with cascade, concatenate, and feedback reduction rules
- Symbolic-to-numeric translation (including time-dependent couplings)
- Utilities for pulse modes, virtual cavities, interaction picture and pulse delay
- Two-time correlation functions and output-mode extraction
- Compatibility with [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) and [QuantumCumulants.jl](https://github.com/qojulia/QuantumCumulants.jl)

## Installation

```julia
|pkg> add QuantumInputOutput
```

## Where to Go

- [Tutorial](@ref) for a complete walkthrough of cavity scattering
- [Theory](@ref) for the input-output formalism with quantum pulses
- [Implementation](@ref) for the symbolic-to-numeric pipeline
- [API](@ref) for the full list of functions
- [Examples](examples/01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3.md) for multiple different usage illustrations
