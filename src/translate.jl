# ── Numeric conversion helpers ──
# `to_numeric` requires the substitution dict to be `AbstractDict{<:QSym}`; the
# empty default `Dict()` is not, so route through the dict-less method instead.
_to_numeric(op, b, operators) =
    isempty(operators) ? to_numeric(op, b) : to_numeric(op, b, operators)

_translate_numeric(
    op,
    b;
    parameter = Dict(),
    operators = Dict{QSym,Any}(),
    op_type = sparse,
) = op_type(_to_numeric(substitute(op, parameter), b, operators))

# Static numeric translation of a full operator sum. SQA's `to_numeric` reduces
# coefficients via `_to_complex`, which fails on constant-but-symbolic values such
# as `sqrt(1.0)` (the result of substituting a number into `√(symbol)`) because
# Symbolics 7 does not constant-fold them. Fold each term's coefficient through
# `_const_coeff` (compiled evaluation) before assembling the operator.
function _translate_numeric(
    op::SQA.QAdd,
    b;
    parameter = Dict(),
    operators = Dict{QSym,Any}(),
    op_type = sparse,
)
    op_ = substitute(op, parameter)
    iszero(op_) && return op_type(0 * _one_op(b, operators))
    result = nothing
    for (term, c_) in op_.arguments
        c = _coeff_num(c_)
        _coeff_is_const(c) || throw(
            ArgumentError(
                "cannot translate `$op` to a static operator: coefficient `$c` still " *
                "depends on a free variable; supply it via `parameter` or `time_parameter`.",
            ),
        )
        contrib = _const_coeff(c) * _numeric_product(term.ops, b, operators, op_type)
        result = result === nothing ? contrib : result + contrib
    end
    return op_type(result)
end

_translate_numeric_raw(op, b; operators = Dict{QSym,Any}(), op_type = sparse) =
    op_type(_to_numeric(op, b, operators))

# Numeric operator product of a `QTerm.ops` vector.
function _numeric_product(ops, b, operators, op_type)
    isempty(ops) && return op_type(_one_op(b, operators))
    return op_type(prod(_to_numeric(o, b, operators) for o in ops))
end

_one_op(b, operators) = isempty(operators) ? one(b) : one(basis(first(values(operators))))

_translate_one(b, operators, op_type) = op_type(_one_op(b, operators))

# Build the operator-substitution dict, optionally adding adjoint entries.
# Always returns a `Dict{QSym,Any}` so it dispatches to `to_numeric(_, _, ::AbstractDict{<:QSym})`.
function _translate_operator_dict(operators, adjoint_ops::Bool)
    out = Dict{QSym,Any}()
    for (k, v) in operators
        out[k] = v
    end
    if adjoint_ops
        for (k, v) in operators
            k_adj = Base.adjoint(k)
            haskey(out, k_adj) || (out[k_adj] = Base.adjoint(v))
        end
    end
    return out
end

# Complex-valued parameters (`Complex{Num}`) are stored in coefficients as separate
# `real`/`imag` parts, so a whole-variable substitution `gv => 0.0` never matches.
# Expand each complex key into independent real- and imaginary-part substitutions;
# the imaginary leg is dropped when it is a literal (e.g. a real variable that was
# promoted to `Complex{Num}` inside a mixed key vector).
function _expand_parameter(parameter)
    isempty(parameter) && return parameter
    out = Dict{Any,Any}()
    for (k, v) in parameter
        if k isa Complex
            out[real(k)] = real(v)
            ik = imag(k)
            SymbolicUtils.unwrap(ik) isa Number || (out[ik] = imag(v))
        else
            out[k] = v
        end
    end
    return out
end

# ── Time-parameter handling ──

# Wrap plain-number values as constant functions; pass functions through.
function _normalize_time_parameter(time_parameter)
    if isempty(time_parameter)
        return time_parameter
    end
    KT = keytype(time_parameter)
    out = Dict{KT,Any}()
    for (k, v) in time_parameter
        out[k] = v isa Number ? (t -> v) : v
    end
    return out
end

# Reduce the time-parameter keys to the underlying symbolic base variables and
# the runtime value of each. Complex coefficients are stored in `real`/`imag`
# form in v0.6, so `build_function` only accepts bare variables as arguments;
# an expression key such as `conj(E)` is inverted to its base variable `E` with
# the value function conjugated accordingly.
function _time_basis(time_parameter)
    basevars = Any[]
    valuefuncs = Any[]
    for (k, f) in time_parameter
        uk = SymbolicUtils.unwrap(k)
        if SymbolicUtils.issym(uk)
            push!(basevars, k)
            push!(valuefuncs, f)
            continue
        end
        vs = collect(Symbolics.get_variables(k))
        length(vs) == 1 || error(
            "time_parameter key `$k` must depend on exactly one variable, got $(length(vs)).",
        )
        wv = Symbolics.wrap(vs[1])
        if isequal(uk, SymbolicUtils.unwrap(conj(wv)))
            push!(basevars, wv)
            push!(valuefuncs, t -> conj(f(t)))
        else
            error(
                "unsupported time_parameter key `$k`; only a bare variable `v` or `conj(v)` are supported.",
            )
        end
    end
    return basevars, Tuple(valuefuncs)
end

# SQA v0.7 stores each `QTerm` prefactor as a `Coeff` (native/poly/symbolic forms);
# lower it to a `Complex{Num}` at the symbolic boundary so the helpers below operate
# on plain Symbolics expressions. A `Complex{Num}` passes through unchanged.
_coeff_num(c) = SQA.to_num(c)
_coeff_num(c::Complex{Num}) = c

# A coefficient with no free variables is concrete.
_coeff_is_const(c::Complex{Num}) =
    isempty(Symbolics.get_variables(real(c))) && isempty(Symbolics.get_variables(imag(c)))

# Reduce a concrete `Complex{Num}` to a plain Julia number. Coefficients that are
# already numeric take a fast path; constant symbolic expressions (e.g. `exp(0.5im)`
# produced by substituting a numeric value into `exp(im*ϕ)`) are compiled and evaluated.
function _const_coeff(c::Complex{Num})
    re = Symbolics.value(SymbolicUtils.unwrap(real(c)))
    im = Symbolics.value(SymbolicUtils.unwrap(imag(c)))
    if re isa Number && im isa Number
        return iszero(im) ? re : complex(re, im)
    end
    # `build_function` returns either a single function or an (out-of-place, in-place)
    # tuple; take the out-of-place callable and evaluate the (variable-free) expression.
    f = build_function(c; expression = Val(false))
    g = f isa Tuple ? first(f) : f
    return g()
end

_as_cnum(x::Complex) = Complex{Num}(Num(real(x)), Num(imag(x)))
_as_cnum(x) = Complex{Num}(Num(x), Num(false))

# Compile a `Complex{Num}` coefficient into a runtime function of the time-parameter
# base variables. A naive `build_function(c, vars...)` emits a `complex(re, im)` call,
# which requires both arguments to be `<:Real` and therefore throws `MethodError:
# complex(::Int64, ::ComplexF64)` when a base variable is fed a complex-valued time
# function. Following QuantumCumulants (`_im_form`/`_qc_maketerm`), we instead compile
# the real and imaginary parts separately and recombine them with `im` in Julia, so the
# generated code uses ordinary `+`/`*` and stays valid for real- and complex-valued
# arguments alike. This is behaviour-preserving for the previously-working cases.
function _compile_coeff(c::Complex{Num}, vars...)
    f_re = build_function(real(c), vars...; expression = Val(false))
    f_im = build_function(imag(c), vars...; expression = Val(false))
    g_re = f_re isa Tuple ? first(f_re) : f_re
    g_im = f_im isa Tuple ? first(f_im) : f_im
    return (vals...) -> g_re(vals...) + im * g_im(vals...)
end

# ── Per-term translation ──
# Translate a single `(ops, coefficient)` term into a time-dependent function
# `t -> op`. A concrete coefficient yields a constant function; a symbolic
# coefficient is compiled against the time-parameter base variables.
function _translate_term(ops, c::Complex{Num}, b, time_parameter, operators, op_type)
    prodop = _numeric_product(ops, b, operators, op_type)
    if _coeff_is_const(c)
        op = _const_coeff(c) * prodop
        return t -> op
    end
    basevars, valuefuncs = _time_basis(time_parameter)
    pref = _compile_coeff(c, basevars...)
    return t -> pref(map(g -> g(t), valuefuncs)...) * prodop
end

"""
    translate_qo(op, b::QuantumOpticsBase.Basis; parameter=Dict(), time_parameter=Dict(),
              operators=Dict(), adjoint_ops=true, op_type=sparse)

Translate a symbolic operator `op` into a numeric QuantumOptics.jl operator with the corresponding basis `b`.
The dictionary `parameter` substitutes symbolic parameters with numbers. Time-dependent functions can be provided
with the dictionary `time_parameter`.
If `time_parameter` is non-empty, the result is a time-dependent function `t -> op(t)`.
The operator type which should be returned can be set with the kwarg `op_type=sparse` and
a list of user-defined operators (e.g. on a different basis than `b`) can be provided with the
dictionary `operators=Dict()`. These operators will then be used to replace the symbolic expressions.
If `adjoint_ops=true`, adjoint entries are added automatically for operators missing from the
dictionary, e.g. `a => A` also provides `a' => A'`.
"""
function translate_qo(
    op,
    b::QuantumOpticsBase.Basis;
    parameter = Dict(),
    time_parameter = Dict(),
    operators = Dict(),
    adjoint_ops = true,
    op_type = sparse,
)
    param = _expand_parameter(parameter)
    tp = _normalize_time_parameter(time_parameter)
    operators_ = _translate_operator_dict(operators, adjoint_ops)
    return _translate_qo(
        op,
        b;
        parameter = param,
        time_parameter = tp,
        operators = operators_,
        op_type,
    )
end

# ── Internal dispatch: time_parameter already normalized ──

function _translate_qo(
    op::SQA.QAdd,
    b::QuantumOpticsBase.Basis;
    parameter = Dict(),
    time_parameter = Dict(),
    operators = Dict{QSym,Any}(),
    op_type = sparse,
)
    if isempty(time_parameter)
        return _translate_numeric(op, b; parameter, operators, op_type)
    end

    op_ = substitute(op, parameter)
    iszero(op_) && return t -> op_type(0 * one(b))

    pairs = collect(op_.arguments)
    if length(pairs) == 1
        (term, c) = pairs[1]
        return _translate_term(
            term.ops,
            _coeff_num(c),
            b,
            time_parameter,
            operators,
            op_type,
        )
    end

    # Multi-term: combine via FunctionWrappers with a common concrete type so the
    # returned closure is type-stable. Inner products use `sparse` for a uniform type.
    first_res = _translate_term(
        pairs[1].first.ops,
        _coeff_num(pairs[1].second),
        b,
        time_parameter,
        operators,
        sparse,
    )
    OpType = typeof(first_res(0.0))
    FW = FunctionWrapper{OpType,Tuple{Float64}}

    args_wrapped = ntuple(length(pairs)) do k
        res =
            k == 1 ? first_res :
            _translate_term(
                pairs[k].first.ops,
                _coeff_num(pairs[k].second),
                b,
                time_parameter,
                operators,
                sparse,
            )
        FW(res)
    end

    return t -> begin
        result = args_wrapped[1](t)
        for i = 2:length(args_wrapped)
            result = result + args_wrapped[i](t)
        end
        result
    end
end

function _translate_qo(
    op::QSym,
    b::QuantumOpticsBase.Basis;
    parameter = Dict(),
    time_parameter = Dict(),
    operators = Dict{QSym,Any}(),
    op_type = sparse,
)
    if isempty(time_parameter)
        return _translate_numeric(op, b; parameter, operators, op_type)
    end
    output_op = _translate_numeric_raw(op, b; operators, op_type)
    return t -> output_op
end

function _translate_qo(
    arg_c_,
    b::QuantumOpticsBase.Basis;
    parameter = Dict(),
    time_parameter = Dict(),
    operators = Dict{QSym,Any}(),
    op_type = sparse,
)
    arg_c = substitute(arg_c_, parameter)
    one_b = _translate_one(b, operators, op_type)
    c = _as_cnum(arg_c)

    if isempty(time_parameter)
        _coeff_is_const(c) || throw(
            ArgumentError(
                "cannot translate symbolic scalar `$arg_c` without a value: supply it via " *
                "`parameter` (numeric) or `time_parameter` (time-dependent).",
            ),
        )
        return _const_coeff(c) * one_b
    end

    if _coeff_is_const(c)
        val = _const_coeff(c)
        return t -> val * one_b
    end
    basevars, valuefuncs = _time_basis(time_parameter)
    pref = _compile_coeff(c, basevars...)
    return t -> pref(map(g -> g(t), valuefuncs)...) * one_b
end

function translate_qo(ops::Vector, b::QuantumOpticsBase.Basis; kwargs...)
    return [translate_qo(op, b; kwargs...) for op in ops]
end

function translate_qo(G::SLH, b::QuantumOpticsBase.Basis; kwargs...)
    L = lindblad(G)
    H = hamiltonian(G)
    H_QO = translate_qo(H, b; kwargs...)
    L_QO = [translate_qo(L_, b; kwargs...) for L_ in L]
    return H_QO, L_QO
end
