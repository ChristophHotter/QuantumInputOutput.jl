# MAIN TODO - LIST:

- [ ] clean up 
  - [x] rm QC.jl 
  - [ ] types
  - [ ] rename translate_qo (to to_numerics ?)
  - [ ] speed up derivation

- [ ] examples
  - [ ] Multi-mode Dicke state superradiance (https://journals.aps.org/pra/abstract/10.1103/PhysRevA.103.033713)
  - [ ] SUPER example (Johannes)
  - [ ] pulse delay advanced example (Victor paper)
  - [x] feedback reduction (see SLH paper)  
  - [ ] cumulants correlation matrix modes 
  
- [ ] create more tests (cov +80%)
  - [x] simple full example
  - [ ] QC.jl: comparison cavity drive with full quantum (example)
  - [ ] translate.jl
  - [ ] utils.jl
  - [ ] codecov
  - [ ] interaction picture: test general A(t)
  - [ ] cumulants: comparison example 04-1 and 04-2
  - [ ] u_Gaussian, v_Gaussian; slower?
  - [ ] look at feedback reduction tests

- [ ] more functionalities
  - [x] feedback reduction
  - [ ] better method for u_eff and v_eff (Victor)
  - [ ] interaction picture (Victor paper)
    - [ ] general expression for A(t)
  - [ ] pulse delay (Victor paper)
  - [ ] two_time_corr_matrix for Ls const (see example 02-2)

- [x] Formatter/SpellCheck

- [ ] JET

- [ ] Aqua


## Additional TODOs

- [ ] SLH more specific types
  - [ ] AbstractVector{<:QTerm}
- [ ] directly use gu(t) without substitution? (think about iscall)
- [ ] u_eff and v_eff equation (PRA 2020, appendix) theory,  API, implementation



## DONE

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
