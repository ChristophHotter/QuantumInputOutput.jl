# MAIN TODO - LIST:

- [ ] clean up code (types, ...)
- [ ] docu
  - [ ] examples
    - [x] .jl file conversion (see QC.jl, Orjan)
    - [x] interaction picture example (Victor paper)
    - [ ] Multi-mode Dicke state superradiance
    - [ ] SUPER example (Johannes)
    - [x] g2 for waveguide example
    - [ ] pulse delay (Victor paper)
    - [ ] beam splitter 
  - [x] theory
  - [x] introduction
  - [ ] API
  - [x] tutorial
  - [ ] implementation
- [ ] create more tests (cov +80%)
  - [x] simple full example
  - [ ] translate.jl
  - [ ] utils.jl
- [ ] more functionalities
  - [ ] conj() and sqrt() automatically (maybe also ()^x, exp(), sin(), cos(), tan(), etc)
  - [ ] interaction picture (Victor)
    - [ ] general expression for A(t)
    - [ ] analytic expression for u=v
    - [ ] kwarg for adjoint replacement in substitute_operators
  - [ ] pulse delay (Victor paper)
  - [x] directly QO.jl objects
  - [ ] padding
  - [ ] feedback reduction
- [ ] More tests
  - [ ] codecov
  - [ ] interaction picture: gamma=0; comparison analytic M; comparison (atom) with "normal" system
- [ ] Formatter/SpellCheck
- [ ] JET
- [ ] Aqua


## Additional TODOs

- [x] Translate function: 
    - [x] Try TimeDependentSums again (after solving the gu(t) mistake/problem) [~factor 2 slower!]
    - [x] Be careful with other operations for time-dependent functions (power, conj, ...) 
    - [x] Use gu(t) in SLH expressions (destinguish with iscall -> arguments -> time )
- [x] Parameter equal to 1.0 
- [x] to_numerics() Lazy or "normal" kwarg #(afer SQA.jl PR merge)
- [x] Kwarg to provide substitution operators 
- [x] level_map 
- [x] output_functions name - test if name matters [names matter!]


- [ ] rename functions:
  - [x] concatenation -> concatenate
  - [ ] ui_to_u_i_im1, vi_to_v_i_im1
  - [x] time_dep_param -> time_parameter
- [x] change README.me
- [ ] SLH more specific types
  - [ ] AbstractVector{<:QTerm}
- [x] SLH for QO.jl objects directly
- [ ] Is QC.jl needed for the package? (or only SQA.jl?)
- [ ] time_parameter allow for (linear) interpolations (not only functions)
- [ ] Get rid of "using SecondQuantizedAlgebra" in the examples: 
  - [ ] Import and export all SQA.jl 
