_translate_numeric(
    op,
    b;
    parameter = Dict(),
    level_map = nothing,
    operators = Dict(),
    op_type = sparse,
) = op_type(to_numeric(substitute(op, parameter), b, operators))

_translate_numeric_raw(op, b; level_map = nothing, operators = Dict(), op_type = sparse) =
    op_type(to_numeric(op, b, operators))

_translate_one(b, operators, op_type) =
    isempty(operators) ? op_type(one(b)) : op_type(one(basis(first(values(operators)))))

function _translate_operator_dict(operators::Dict, adjoint_ops::Bool)
    if isempty(operators) || !adjoint_ops
        return Dict{QSym,Any}(operators)
    end

    operators_ = Dict{QSym,Any}(operators)
    for (k, v) in operators
        k_adj = Base.adjoint(k)
        if !haskey(operators_, k_adj)
            operators_[k_adj] = Base.adjoint(v)
        end
    end
    return operators_
end

# simplified, no applicable(), values are Number → wrap, anything else → pass through
function _normalize_time_parameter(time_parameter)
    if isempty(time_parameter)
        return time_parameter
    end
    KT = keytype(time_parameter)
    out = Dict{KT,Any}()
    for (k, v) in time_parameter
        vf = v isa Number ? (t -> v) : v
        out[k] = vf
        k_conj = conj(k)
        if !isequal(k_conj, k) && !haskey(out, k_conj)
            out[k_conj] = t -> conj(vf(t))
        end
    end
    return out
end

_numeric_prefactor(x::Number) = x
function _materialize_symbolic_number(x)
    xu = SymbolicUtils.unwrap(x)
    try
        return Core.eval(@__MODULE__, SymbolicUtils.Code.toexpr(xu))
    catch
        return xu
    end
end
_numeric_prefactor(x::BasicSymbolic) = _materialize_symbolic_number(x)
_numeric_prefactor(x::Symbolics.Num) = _materialize_symbolic_number(x)
_numeric_prefactor(x::Complex) = complex(_numeric_prefactor(real(x)), _numeric_prefactor(imag(x)))
_numeric_prefactor(x) = x

_unwrap_subs_dict(subs) = Dict(SymbolicUtils.unwrap(k) => v for (k, v) in subs)

_eval_prefactor(x::Real, subs) = x
_eval_prefactor(x::BasicSymbolic, subs) =
    _numeric_prefactor(SymbolicUtils.substitute(SymbolicUtils.unwrap(x), subs))
_eval_prefactor(x::Symbolics.Num, subs) =
    _numeric_prefactor(SymbolicUtils.substitute(SymbolicUtils.unwrap(x), subs))
_eval_prefactor(x::Complex, subs) =
    _numeric_prefactor(
        SymbolicUtils.substitute(
            SymbolicUtils.unwrap(real(x)) + im * SymbolicUtils.unwrap(imag(x)),
            subs,
        ),
    )
_eval_prefactor(x, subs) = _numeric_prefactor(substitute(x, subs))

_is_concrete_prefactor(x::Real) = !(x isa BasicSymbolic || x isa Symbolics.Num)
_is_concrete_prefactor(x::Complex) =
    _is_concrete_prefactor(real(x)) && _is_concrete_prefactor(imag(x))
_is_concrete_prefactor(x) = x isa Number

# use Tuple + map instead of generator splat to avoid per-call allocation
function _translate_prefactor(arg_c, time_parameter)
    if _is_concrete_prefactor(arg_c)
        return false, arg_c
    end
    pref_f = t -> begin
        subs = _unwrap_subs_dict(Dict(k => v(t) for (k, v) in time_parameter))
        _eval_prefactor(arg_c, subs)
    end
    return true, pref_f
end

_to_numeric_qsym(arg::QSym, b, operators) =
    isempty(operators) ? to_numeric(arg, b) : (haskey(operators, arg) ? operators[arg] : to_numeric(arg, b))

function _translate_qadd_term(term_ops, arg_c, b, time_parameter, level_map, operators, op_type)
    numeric_term = if isempty(term_ops)
        _translate_one(b, operators, op_type)
    else
        op_type(prod(_to_numeric_qsym(arg, b, operators) for arg in term_ops))
    end

    is_func, pref = _translate_prefactor(arg_c, time_parameter)
    if is_func
        return t -> pref(t) * numeric_term
    end
    return pref * numeric_term
end

"""
    translate_qo(op, b::QuantumOpticsBase.Basis; parameter=Dict(), time_parameter=Dict(),
              level_map=nothing, operators=Dict(), adjoint_ops=true, op_type=sparse)

Translate a symbolic operator `op` into a numeric QuantumOptics.jl operator with the corresponding basis `b`.
The dictionary `parameter` substitutes symbolic parameters with numbers. Time-dependent functions can be provide
with the dictionary `time_parameter`.
If `time_parameter` is non-empty, the result is a time-dependent function `t -> op(t)`.
The kwarg `level_map=nothing` is used to provide the names of levels for `transition` operators.
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
    level_map = nothing,
    operators = Dict(),
    adjoint_ops = true,
    op_type = sparse,
)
    tp = _normalize_time_parameter(time_parameter)
    operators_ = _translate_operator_dict(operators, adjoint_ops)
    return _translate_qo(
        op,
        b;
        parameter,
        time_parameter = tp,
        level_map,
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
    level_map = nothing,
    operators = Dict(),
    op_type = sparse,
)
    op = expand(substitute(op, parameter))
    if isempty(time_parameter)
        return _translate_numeric_raw(op, b; level_map, operators, op_type)
    end

    terms = collect(op.arguments)
    isempty(terms) && return op_type(0 * one(b))

    # Translate first arg to determine concrete operator type
    first_translated = _translate_qadd_term(
        first(terms).first.ops,
        first(terms).second,
        b,
        time_parameter,
        level_map,
        operators,
        sparse,
    )
    OpType =
        typeof(first_translated isa Function ? first_translated(0.0) : first_translated)
    FW = FunctionWrapper{OpType,Tuple{Float64}}

    args_wrapped = ntuple(length(terms)) do k
        a_k = if k == 1
            first_translated
        else
            _translate_qadd_term(
                terms[k].first.ops,
                terms[k].second,
                b,
                time_parameter,
                level_map,
                operators,
                sparse,
            )
        end
        FW(a_k isa Function ? a_k : (_ -> a_k))
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
    level_map = nothing,
    operators = Dict(),
    op_type = sparse,
)
    if isempty(time_parameter)
        return _translate_numeric(op, b; parameter, level_map, operators, op_type)
    end
    return _translate_qo(
        substitute(op, parameter),
        b;
        time_parameter,
        level_map,
        operators,
        op_type,
    )
end

function _translate_qo(
    arg_c_,
    b::QuantumOpticsBase.Basis;
    parameter = Dict(),
    time_parameter = Dict(),
    level_map = nothing,
    operators = Dict(),
    op_type = sparse,
)
    arg_c = substitute(arg_c_, parameter)
    one_b = _translate_one(b, operators, op_type)

    if isempty(time_parameter)
        return arg_c * one_b
    end

    is_func, pref = _translate_prefactor(arg_c, time_parameter)
    if is_func
        return t -> pref(t) * one_b
    end
    prod_c_nc_num = pref * one_b
    return t -> prod_c_nc_num
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
