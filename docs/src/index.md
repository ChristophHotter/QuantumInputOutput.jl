# Introduction

**QuantumInputOutput.jl** is a Julia framework for modeling the input-output theory with quantum pulses, using SLH (scattering, Lindblad, Hamiltonian) elements and rules. It combines:
- a symbolic layer based on [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl) to build models using SLH rules
- a numerical layer based on [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) and [QuantumCumulants.jl](https://github.com/qojulia/QuantumCumulants.jl) to simulate time dynamics and observables in a full quantum or higher-order meanfield approach, respectively

The typical workflow is:
1. Build the SLH model symbolically
2. Translate to numerical operators
3. Evolve the system in time
4. Analyze output modes and correlations

## Key Features

- SLH modeling with cascade and concatenate rules
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
