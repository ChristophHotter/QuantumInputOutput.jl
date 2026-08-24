# TODOs

- [ ] examples
  - [x] SUPER example 
  - [ ] pulse delay advanced example
  - [x] feedback reduction (see SLH paper)  
  - [ ] cumulants correlation matrix modes 
  - [ ] H_loop in example 09-1 is not simplified correctly on the docu (locally it is fine; pkg versions)
  - [ ] Multi-mode Dicke state superradiance (https://journals.aps.org/pra/abstract/10.1103/PhysRevA.103.033713)
  
- [ ] create more tests (cov +80%)
  - [x] simple full example
  - [ ] QC.jl: comparison cavity drive with full quantum (example)
  - [x] translate.jl
  - [x] utils.jl
  - [x] codecov
  - [ ] interaction picture: test general A(t)
  - [ ] cumulants: comparison example 04-1 and 04-2

- [ ] more functionalities
  - [ ] pulse delay 
  - [ ] better method for u_eff and v_eff (Victor)

## Additional TODOs

- [x] u_eff and v_eff equation (PRA 2020, appendix) theory, API, implementation
- [x] rename translate_qo (now `to_numeric`, from SQA v0.9; `substitute_operators` now `substitute`)


### QIO maintainability wins
- [ ] Examples in CI as smoke tests (highest leverage). Examples were untested and rotted: stale API, several did not run. Real bugs (04-2 MTK flow, `get_solution` new signature, complex-coupling crash, `QAdd + BasicSymbolic` scalar drive) were only caught by running them by hand. A CI job running each example's symbolic setup plus a tiny ODE (skip plots/long solves) would catch API drift.
- [ ] Document the `::Real` vs `::Complex` contract for couplings. The `conj` in `g'` is load-bearing: a coupling that is physically complex must be `::Complex` or the physics is silently wrong. Old `rnumbers` made couplings real; that was a latent modeling bug. Consider a runtime check ("real parameter received a complex value").
- [ ] Harden SLH `_add` / `_mul` for scalar couplings: `QAdd + BasicSymbolic` (a classical scalar drive like `Et(t)`) was not handled.
- [ ] Rendering: `Complex{Num}` decomposes into `real(...)` / `imag(...)` inside products, so complex-coupling Hamiltonians print messy. This is an SQA display representation concern, not fixable in QIO; raise with SQA if it matters for the docs.

## DONE
- [x] rm QC.jl 
- [x] types
- [x] speed up derivation

- [x] Translate function: 
    - [x] Try TimeDependentSums again (after solving the gu(t) mistake/problem) [~factor 2 slower!]
    - [x] Be careful with other operations for time-dependent functions (power, conj, ...) 
    - [x] Use gu(t) in SLH expressions (distinguish with iscall -> arguments -> time )
  
- [x] Parameter equal to 1.0 
- [x] to_numerics() Lazy or "normal" kwarg #(after SQA.jl PR merge)
- [x] Kwarg to provide substitution operators 
- [x] level_map 
- [x] output_functions name - test if name matters [names matter!]

- [x] rename functions:
  - [x] concatenation -> concatenate
  - [x] u_eff, v_eff
  - [x] time_dep_param -> time_parameter

- [x] change README.me
- [x] SLH for QO.jl objects directly

- [x] Docu review
  - [x] Implementation.md more code snippet 
  - [x] concatenate docstring
  - [x] check interaction_picture_A_4modes docstrings

- [x] examples
  - [x] g2 for waveguide example
  - [x] .jl file conversion (see QC.jl, Orjan)
  - [x] interaction picture example (Victor paper)
  - [x] beam splitter loss
  - [x] beam splitter HOM
- [x] theory
- [x] introduction
- [x] API
- [x] tutorial
- [x] implementation

- [x] more functionalities
  - [x] conj() and sqrt() automatically (maybe also ()^x, exp(), sin(), cos(), tan(), etc)
  - [x] directly QO.jl objects
  - [x] padding
  - [x] gu_to_Gauss
    - [x] add Δ as kwarg
    - [x] arg u not needed?! 
  - [x] analytic expression for u=v
  - [x] kwarg for adjoint replacement in substitute_operators
  - [x] feedback reduction
  - [x] interaction picture 
    - [x] general expression for A(t)
  - [x] two_time_corr_matrix for Js const (see example 02-2)

- [x] Formatter/SpellCheck
- [x] JET
- [x] Aqua
  
- [x] docu
  - [x] u_eff and v_eff theory 
  - [x] development status
