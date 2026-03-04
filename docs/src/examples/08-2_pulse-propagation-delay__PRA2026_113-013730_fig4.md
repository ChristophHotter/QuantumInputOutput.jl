```@meta
EditURL = "../../../examples/08-2_pulse-propagation-delay__PRA2026_113-013730_fig4.jl"
```

# TODO

# Ramsey Interference with Delayed Fock Pulses

In this example, we reproduces the Ramsey-like interference pattern for a partially delayed input Fock state with `n = 9`, studied in  [V. R. Christiansen and K. Mølmer, Phys. Rev. A 113, 013730 (2026)](https://doi.org/10.1103/PhysRevA.113.013730).
We model a single input Fock pulse that is split on a 50/50 beam splitter (with vacuum padding on the second port),
creating two channels that interact with the two-level system. We then scan the detuning to extract the excited-state
population at the evaluation time used in the paper and compare with a coherent-state drive of the same mean photon number.

````@example 08-2_pulse-propagation-delay__PRA2026_113-013730_fig4
using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using LinearAlgebra
using PyPlot
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

