using QuantumInputOutput
using SecondQuantizedAlgebra
using SymbolicUtils
using Symbolics
QIO = QuantumInputOutput
using BenchmarkTools

@cnumbers g t
@variables gt2(..) # symbolic function placeholder

expr = 3 * sqrt(2 * conj(g))
# expr = 3 * sqrt(2 * g)
gt2_numeric(t::Real) = t>1 ? (7t + 1im) : 9

### works 
time_parameter = Dict(g => gt2(t))
expr_sub = substitute(expr, time_parameter)
f_ = build_function(expr_sub, t, gt2; expression=Val(true)) # works!! # TODO
result = f_(2.0, gt2_numeric)
3 * sqrt(2 * conj(7*2 + 1im))
f_new_sym = t -> f_(t, gt2_numeric)
f_new_sym(2.0)
@benchmark f_new_sym(2.0)
###

f = build_function(expr, t, g; expression=Val(true))
f_new = t -> f(t, gt2_numeric(t))
f_new(2.0)

f_def(t) = 3 * sqrt(t * conj(gt2_numeric(t)))
f_def(2.0)


@benchmark f_new(2.0)

@benchmark f_def(2.0)


####
@rnumbers gg

f2 = build_function(expr, t, g, gg; expression=Val(true))
f_new2 = t -> f2(t, gt2_numeric(t), 0)
f_new2(2.0)

@benchmark f_new2(2.0)

time_parameter