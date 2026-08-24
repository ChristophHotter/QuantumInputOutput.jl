<p align="center">
  <img src="assets/QuantumInputOutput-logo-gaussian.svg" alt="QuantumInputOutput.jl logo" width="320">
</p>

# QuantumInputOutput.jl

**QuantumInputOutput.jl** is a Julia framework for modeling the input-output formalism with quantum pulses.
It combines a symbolic layer for the SLH formalism based on [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl) with numerical solvers from [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) and [QuantumCumulants.jl](https://github.com/qojulia/QuantumCumulants.jl), enabling simulation of systems interacting with quantum pulses.

## Key Features

- SLH modeling with cascade (`▷`), concatenate (`⊞`), and feedback reduction rules
- Symbolic-to-numeric translation (including time-dependent couplings)
- Virtual-cavity tools for temporal input and output modes
- Two-time correlation functions and output-mode extraction
- Compatibility with [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) and [QuantumCumulants.jl](https://github.com/qojulia/QuantumCumulants.jl) solvers

## Development status

[![Documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://qojulia.github.io/QuantumInputOutput.jl/stable/)
[![Documentation dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://qojulia.github.io/QuantumInputOutput.jl/dev/)
[![Tests](https://github.com/qojulia/QuantumInputOutput.jl/actions/workflows/Tests.yml/badge.svg)](https://github.com/qojulia/QuantumInputOutput.jl/actions/workflows/Tests.yml)
[![codecov](https://codecov.io/gh/qojulia/QuantumInputOutput.jl/graph/badge.svg)](https://codecov.io/gh/qojulia/QuantumInputOutput.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![JET](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)

Note that **QuantumInputOutput.jl** is still at an early stage of development.

## Installation

```julia
|pkg> add QuantumInputOutput
```

## Documentation

Read the [stable documentation](https://qojulia.github.io/QuantumInputOutput.jl/stable/)
or the [development documentation](https://qojulia.github.io/QuantumInputOutput.jl/dev/).

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
- a brief description of the problem or feature
- minimal reproducible example or tests if applicable
