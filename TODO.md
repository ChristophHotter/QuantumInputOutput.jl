# TODOs

- [ ] examples
  - [ ] SUPER example 
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
- [ ] rename translate_qo (to to_numerics ?)



## Architecture & maintainability (review after SQA 0.6 / QC 0.5 migration)

### Layering: SQA vs QC vs QIO
Clean intended split: SQA owns symbolic operator algebra plus the numeric backend; QC owns dynamics derivation (meanfield/cumulants) plus the ModelingToolkit ODE bridge; QIO owns SLH network composition plus the pulse/input-output physics. The split is mostly real. Every leak is in one spot: the "symbolic to numeric coefficient" bridge.

Every hard bug in the migration was the same bug rediscovered at a different layer ("reduce or compile a `Complex{Num}` coefficient under Symbolics 7"):
- `to_numeric` not constant-folding `sqrt(1.0)` (worked around in QIO `_translate_numeric`).
- `build_function` emitting `complex(re, im)` which throws on complex inputs (QC already solved this with `_im_form` / `_qc_maketerm`; QIO re-solved it with `_compile_coeff`).
- Substituting a `Complex{Num}` parameter key not matching (QIO `_expand_parameter`).
- `substitute_operators`: SQA `substitute` stack-overflows on self-referential operator replacement (`a => 2a`); QIO hand-rolled dict-lookup; QC already has robust `rewrite`.

Conclusion: QC and QIO each carry a private copy of the coefficient-codegen machinery and each independently ate the Symbolics 7 upgrade pain. That duplication is the core smell.

### translate_qo really belongs in SQA (completing the QO bridge)
SQA already translates to QuantumOptics (`to_numeric` in `SQA/src/numeric.jl`, already depends on QuantumOpticsBase). `translate_qo` is not a new capability; it is `to_numeric` completed with the three things `to_numeric` punts on: parameter substitution, constant-folding, and time-dependence. Symmetry: SQA to QuantumOptics is one numeric backend (currently shipped half in SQA, half stranded in QIO); QC to ModelingToolkit is the other (fully owned by QC). The complete QO backend is SQA's to own.

- [ ] Move the generic bridge into SQA (ideally a QuantumOptics package extension so SQA core stays algebra-only): `to_numeric(op, basis; parameters, time_parameters)` returning an operator or a `t -> op(t)`; parameter substitution incl. `Complex{Num}` keys; constant-folding; the `im`-recombination coefficient compiler (one robust copy).
- [ ] Keep in QIO only the domain-specific bits: the `SLH` overload, the `operators` kwarg (substitute symbolic ops with user QO ops on another basis), the `adjoint_ops` convenience, and pulse/interaction-picture/correlation physics. QIO `translate_qo` then becomes a thin SLH wrapper over the SQA primitive.
- [ ] This collapses QC `_im_form` and QIO `_compile_coeff` into one shared path; fix the Symbolics-version pain once, not N times.
- Caveats: `t -> op(t)` is solver-flavored (be deliberate that the abstraction is "time-dependent operator", not "thing for master_dynamic"); time path pulls `FunctionWrappers` + `build_function` (small, QOBase already a dep); this is a qojulia-wide scope decision for the SQA maintainers.
- Related: the existing "rename translate_qo (to to_numerics ?)" item above.

### Things to fix in SQA (or lift from QC)
- [ ] `to_numeric` should accept a parameter dict and constant-fold residual symbolic constants (`sqrt`, `exp`), instead of erroring in `_to_complex`.
- [ ] Robust operator substitution: fix SQA `substitute` for self-referential replacement, or lift QC `rewrite`/`_qc_maketerm` into SQA; then QIO `substitute_operators` is a thin wrapper.

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
  - [x] two_time_corr_matrix for Ls const (see example 02-2)

- [x] Formatter/SpellCheck
- [x] JET
- [x] Aqua
  
- [x] docu
  - [x] u_eff and v_eff theory 
  - [x] development status
