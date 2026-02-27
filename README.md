# QuantumInputOutput.jl

**QuantumInputOutput.jl** is a Julia framework for modeling the input-output formalism with quantum pulses.
It combines a symbolic layer for the SLH formalism based on `SecondQuantizedAlgebra.jl` with numerical solvers from `QuantumOptics.jl` and `QuantumCumulants.jl`, enabling simulation of systems interacting quantum pulses.

## Key Features

- SLH modeling with cascade (`▷`) and concatenation (`⊞`) rules
- Symbolic-to-numeric translation (including time-dependent couplings)
- Virtual-cavity tools for temporal input and output modes
- Two-time correlation functions and output-mode extraction
- Compatibility with `QuantumOptics.jl` solvers

## Installation

```julia
|pkg> add QuantumInputOutput
```

## Documentation

The documentation is built with Documenter.jl and includes:
- Theory: background on SLH and virtual-cavity pulse modeling
- Tutorial: a full cavity-scattering walkthrough
- Implementation notes
- API reference
- Examples

To build locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
```

## Contributing

Contributions are welcome. Please open an issue or a pull request with:
- a brief description of the problem or feature,
- minimal reproducible example or tests if applicable.

