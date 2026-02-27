# Introduction

**QuantumInputOutput.jl** is a Julia framework for modeling the input-output theory with quantum pulses. It combines:
- a symbolic layer based on `SecondQuantizedAlgebra.jl` to build models using SLH rules, and
- a numerical layer based on `QuantumOptics.jl` and `QuantumCumulants.jl` to simulate time dynamics and observables in a full quantum or higher-order meanfield approach, respectively.

The typical workflow is:
1. Build the SLH model symbolically.
2. Translate to numerical operators.
3. Evolve the system in time.
4. Analyze output modes and correlations.

## Key Features

- SLH modeling with cascade and concatenate rules
- Symbolic-to-numeric translation (including time-dependent couplings)
- Utilities for pulse modes and virtual cavities
- Two-time correlation functions and output-mode extraction
- Compatibility with `QuantumOptics.jl` solvers

## Installation

```julia
|pkg> add QuantumInputOutput
```

## Where to Go

- [Tutorial](@ref) for a complete walkthrough of cavity scattering
- [Theory](@ref) for the input-output formalism with quantum pulses
- [Implementation](@ref) for the symbolic-to-numeric pipeline
- [API](@ref) for the full list of functions
- `Examples` for multiple different usage illustrations # TODO: reference?
