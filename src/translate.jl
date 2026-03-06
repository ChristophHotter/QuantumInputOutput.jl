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

function _normalize_time_parameter(time_parameter)
    if isempty(time_parameter)
        return time_parameter
    end
    out = Dict()
    for (k, v) in time_parameter
        if isa(v, Function)
            out[k] = v
        elseif isa(v, Number)
            out[k] = t -> v
        else
            error(
                "time_parameter values must be callable or numeric, got $(typeof(v)) for key $(k)",
            )
        end
    end
    return out
end

function _translate_prefactor(arg_c, time_parameter)
    keys_ = collect(keys(time_parameter))
    values_ = collect(values(time_parameter))

    pref_build = build_function(arg_c, keys_...; expression = Val(false))
    pref_f = t -> pref_build((v(t) for v in values_)...)
    return true, pref_f
end
function _translate_prefactor(arg_c::Number, time_parameter)
    return false, arg_c
end

"""
    translate(op, b::QuantumOpticsBase.Basis; parameter=Dict(), time_parameter=Dict(),
              level_map=nothing, operators=Dict(), op_type=sparse)

Translate a symbolic operator `op` into a numeric QuantumOptics.jl operator with the corresponding basis `b`. 
The dictionary `parameter` substitutes symbolic parameters with numbers. Time-dependent functions can be provide 
with the dictionary `time_parameter`. 
If `time_parameter` is non-empty, the result is a time-dependent function `t -> op(t)`. 
The kwarg `level_map=nothing` is used to provide the names of levels for `transition` operators. 
The operator type which should be returned can be set with the kwarg `op_type=sparse` and  
a list of user-defined operators (e.g. on a different basis than `b`) can be provided with the 
dictionary `operators=Dict()`. These operators will then be used to replace the symbolic expressions.  
"""
function translate(
    op::SQA.QMul,
    b::QuantumOpticsBase.Basis;
    parameter = Dict(),
    time_parameter = Dict(),
    level_map = nothing,
    operators = Dict(),
    op_type = sparse,
)

    time_parameter = _normalize_time_parameter(time_parameter)
    if isempty(time_parameter)
        return _translate_numeric(
            op,
            b;
            parameter = parameter,
            level_map = level_map,
            operators = operators,
            op_type = op_type,
        ) # TODO: test
    end

    op_ = substitute(op, parameter) # First substitute all symbolic parameters with the numerical values -> only time-dependent parameters left

    if isa(op_, QSym) # special case if parameter is 1.0: after substitute p*op -> 1.0*op -> op (has no fields arg_c, args_nc)
        output_op_QMul_QO_ = _translate_numeric_raw(
            op_,
            b;
            level_map = level_map,
            operators = operators,
            op_type = op_type,
        ) # this is faster!
        output_op_QMul_QO = t -> output_op_QMul_QO_
        return output_op_QMul_QO
    elseif iszero(op_) # special case if parameter is 0.0
        return op_type(0 * one(b))
    else
        arg_c = op_.arg_c
        args_nc = op_.args_nc
        prod_args_nc = op_type(
            prod((to_numeric(arg, b, operators; level_map = level_map)) for arg in args_nc),
        )

        is_func, pref = _translate_prefactor(arg_c, time_parameter)
        if is_func
            return t -> pref(t) * prod_args_nc
        end
        prod_c_nc_num = pref * prod_args_nc
        return t -> prod_c_nc_num
    end
end
#
function translate(
    op::SQA.QAdd,
    b::QuantumOpticsBase.Basis;
    parameter = Dict(),
    time_parameter = Dict(),
    level_map = nothing,
    operators = Dict(),
    op_type = sparse,
)

    time_parameter = _normalize_time_parameter(time_parameter)
    op = expand(op) # should ensure that only multiplications are in arguments
    ### only static operators - returns a time-independent Hamiltonian and Lindblad
    if isempty(time_parameter)
        return _translate_numeric(
            op,
            b;
            parameter = parameter,
            level_map = level_map,
            operators = operators,
            op_type = op_type,
        ) # TODO: test
    end

    ### static and time-dependent operators
    args = arguments(substitute(op, parameter))
    args_translated = [
        translate(
            arg,
            b;
            parameter = parameter,
            time_parameter = time_parameter,
            level_map = level_map,
            operators = operators,
            op_type = sparse,
        ) for arg in args
    ]

    output_op_QAdd_QO = t -> sum(args_translated[i](t) for i = 1:length(args_translated))
    return output_op_QAdd_QO
end

function translate(
    op::QSym,
    b::QuantumOpticsBase.Basis;
    parameter = Dict(),
    time_parameter = Dict(),
    level_map = nothing,
    operators = Dict(),
    op_type = sparse,
)
    # QSym are only fundamental symbolic operators, e.g. a, σ_-, x, ...

    time_parameter = _normalize_time_parameter(time_parameter)
    if isempty(time_parameter)
        return _translate_numeric(
            op,
            b;
            parameter = parameter,
            level_map = level_map,
            operators = operators,
            op_type = op_type,
        ) # TODO: test
    end

    output_op_QSym_QO_ = _translate_numeric_raw(
        op,
        b;
        level_map = level_map,
        operators = operators,
        op_type = op_type,
    ) # this is faster!
    output_op_QSym_QO = t -> output_op_QSym_QO_
    return output_op_QSym_QO
end

function translate(
    arg_c_,
    b::QuantumOpticsBase.Basis;
    parameter = Dict(),
    time_parameter = Dict(),
    level_map = nothing,
    operators = Dict(),
    op_type = sparse,
)
    # should only be needed for numbers and symbolic parameters

    arg_c = substitute(arg_c_, parameter)
    one_b = _translate_one(b, operators, op_type)

    time_parameter = _normalize_time_parameter(time_parameter)
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

function translate(ops::Vector, b::QuantumOpticsBase.Basis; kwargs...)
    return [translate(op, b; kwargs...) for op in ops]
end
function translate(G::SLH, b::QuantumOpticsBase.Basis; kwargs...)
    L = get_lindblad(G);
    H = get_hamiltonian(G)
    H_QO = translate(H, b; kwargs...)
    L_QO = [translate(L_, b; kwargs...) for L_ in L]
    return H_QO, L_QO
end
