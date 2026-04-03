_translate_numeric(
    op,
    b;
    parameter = Dict(),
    level_map = nothing,
    operators = Dict(),
    op_type = sparse,
) = op_type(to_numeric(substitute(op, parameter), b, operators; level_map = level_map))

_translate_numeric_raw(op, b; level_map = nothing, operators = Dict(), op_type = sparse) =
    op_type(to_numeric(op, b, operators; level_map = level_map))

_translate_one(b, operators, op_type) =
    isempty(operators) ? op_type(one(b)) :
    op_type(one(basis(collect(values(operators))[1])))

# [I4+I5 fix]: simplified, no applicable(), values are Number → wrap, anything else → pass through
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

# [C5 fix]: use Tuple + map instead of generator splat to avoid per-call allocation
function _translate_prefactor(arg_c, time_parameter)
    keys_ = collect(keys(time_parameter))
    values_tuple = Tuple(values(time_parameter))
    pref_build = build_function(arg_c, keys_...; expression = Val(false))
    pref_f = t -> pref_build(map(v -> v(t), values_tuple)...)
    return true, pref_f
end
function _translate_prefactor(arg_c::Number, time_parameter)
    return false, arg_c
end

"""
    translate_qo(op, b::QuantumOpticsBase.Basis; parameter=Dict(), time_parameter=Dict(),
              level_map=nothing, operators=Dict(), op_type=sparse)

Translate a symbolic operator `op` into a numeric QuantumOptics.jl operator.
Time-dependent parameters are normalized ONCE at this entry point.
"""
function translate_qo(
    op,
    b::QuantumOpticsBase.Basis;
    parameter = Dict(),
    time_parameter = Dict(),
    level_map = nothing,
    operators = Dict(),
    op_type = sparse,
)
    tp = _normalize_time_parameter(time_parameter)
    return _translate_qo(op, b; parameter, time_parameter = tp, level_map, operators, op_type)
end

# ── Internal dispatch: time_parameter already normalized ──

function _translate_qo(
    op::SQA.QMul,
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

    op_ = substitute(op, parameter)

    if isa(op_, QSym)
        output_op = _translate_numeric_raw(op_, b; level_map, operators, op_type)
        return t -> output_op
    elseif iszero(op_)
        return op_type(0 * one(b))
    else
        arg_c = op_.arg_c
        args_nc = op_.args_nc
        prod_args_nc = op_type(
            prod((to_numeric(arg, b, operators; level_map)) for arg in args_nc),
        )
        is_func, pref = _translate_prefactor(arg_c, time_parameter)
        if is_func
            return t -> pref(t) * prod_args_nc
        end
        prod_c_nc_num = pref * prod_args_nc
        return t -> prod_c_nc_num
    end
end

function _translate_qo(
    op::SQA.QAdd,
    b::QuantumOpticsBase.Basis;
    parameter = Dict(),
    time_parameter = Dict(),
    level_map = nothing,
    operators = Dict(),
    op_type = sparse,
)
    op = expand(op)
    if isempty(time_parameter)
        return _translate_numeric(op, b; parameter, level_map, operators, op_type)
    end

    args = arguments(substitute(op, parameter))
    args_translated = Tuple(
        _translate_qo(arg, b; parameter, time_parameter, level_map, operators, op_type = sparse)
        for arg in args
    )

    n = length(args_translated)
    return t -> begin
        result = args_translated[1](t)
        for i in 2:n
            result = result + args_translated[i](t)
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
    output_op = _translate_numeric_raw(op, b; level_map, operators, op_type)
    return t -> output_op
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
