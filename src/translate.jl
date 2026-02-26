_translate_numeric(op, b; parameter=Dict(), level_map=nothing, operators=Dict(), op_type=sparse) =
    op_type(to_numeric(substitute(op, parameter), b, operators; level_map=level_map))

_translate_numeric_raw(op, b; level_map=nothing, operators=Dict(), op_type=sparse) =
    op_type(to_numeric(op, b, operators; level_map=level_map))

_translate_one(b, operators, op_type) =
    isempty(operators) ? op_type(one(b)) : op_type(one(basis(collect(values(operators))[1])))

function _normalize_time_dep_param(time_dep_param)
    if isempty(time_dep_param)
        return time_dep_param
    end
    out = Dict()
    for (k, v) in time_dep_param
        if isa(v, Function) 
            out[k] = v
        elseif isa(v, Number)
            out[k] = t -> v
        else
            error("time_dep_param values must be callable or numeric, got $(typeof(v)) for key $(k)")
        end
    end
    return out
end

function _translate_prefactor(arg_c, time_dep_param)
    if !iscall(arg_c) # only one symbolic function in arg_c_ #TODO: time-dep parameter, e.g. g(t)
        arg_c_sub = substitute(arg_c, time_dep_param)
        if isa(arg_c_sub, Function)
            return true, arg_c_sub
        else
            return false, arg_c_sub
        end
    end

    arg_c_args = arguments(arg_c)
    arg_c_op = operation(arg_c)

    # TODO: directly use gu(t) without substitution? (think about iscall)
    if length(arg_c_args) == 1  # conj(g)*a, sqrt(g)*a # TODO: g^2*a [length(arg_c_args) == 2]
        if arg_c_op ∈ [conj, sqrt]
            arg_c_1 = substitute(arg_c, time_dep_param)
            return true, arg_c_1
        end
    else
        !(arg_c_op == *) && error("The operation of the prefactor is $(arg_c_op). 
        Currently, if you want to use a function like conj() or sqrt() on a time dependent 
        parameter, you need to include it in the time-dependent function dictionary, see e.g.: [REF-example]
        At the moment we are also limited to functions with only one argument.") # TODO
    end
    arg_c_all = [substitute(arg, time_dep_param) for arg ∈ arg_c_args]
    arg_c_numbers = filter(x -> isa(x, Number), arg_c_all)
    arg_c_functions = filter(x -> !isa(x, Number), arg_c_all)
    pref_num = prod(arg_c_numbers)
    if isempty(arg_c_functions)
        return false, pref_num
    end
    pref_f = t -> prod(arg_c_f(t) for arg_c_f in arg_c_functions) * pref_num
    return true, pref_f
end

function translate(op::QuantumCumulants.QMul, b::QuantumOpticsBase.Basis;
    parameter=Dict(), time_dep_param=Dict(), level_map=nothing, operators=Dict(), op_type=sparse)

    time_dep_param = _normalize_time_dep_param(time_dep_param)
    if isempty(time_dep_param)
        return _translate_numeric(op, b; parameter=parameter, level_map=level_map, operators=operators, op_type=op_type) # TODO: test
    end

    op_ = substitute(op, parameter) # First substitute all symbolic paramters with the numerical values -> only time-dependent parameters

    if isa(op_, QuantumCumulants.QSym) # special case if parameter is 1.0: after substitute p*op -> 1.0*op -> op (has no fields arg_c, args_nc)
        output_op_QMul_QO_ = _translate_numeric_raw(op_, b; level_map=level_map, operators=operators, op_type=op_type) # this is faster!
        output_op_QMul_QO = t -> output_op_QMul_QO_
        return output_op_QMul_QO
    elseif iszero(op_) # special case if parameter is 0.0
        return op_type(0 * one(b))
    else
        arg_c = op_.arg_c
        args_nc = op_.args_nc
        prod_args_nc = op_type(prod((to_numeric(arg, b, operators; level_map=level_map)) for arg in args_nc)) #TODO: dense, sparse, etc

        is_func, pref = _translate_prefactor(arg_c, time_dep_param)
        if is_func
            return t -> pref(t) * prod_args_nc
        end
        prod_c_nc_num = pref * prod_args_nc
        return t -> prod_c_nc_num
    end
end
#
function translate(op::QuantumCumulants.QAdd, b::QuantumOpticsBase.Basis; 
        parameter=Dict(), time_dep_param=Dict(), level_map=nothing, operators=Dict(), op_type=sparse)  

    time_dep_param = _normalize_time_dep_param(time_dep_param)
    op = expand(op) # should ensure that only multiplications are in arguments
    ### only static operators - returns a time-independent Hamiltonian and Lindblad
    if isempty(time_dep_param)
        return _translate_numeric(op, b; parameter=parameter, level_map=level_map, operators=operators, op_type=op_type) # TODO: test
    end

    ### static and time-dependent operators
    args = arguments(substitute(op, parameter)) 
    args_translated = [translate(arg, b; parameter=parameter, time_dep_param=time_dep_param, level_map=level_map, operators=operators, op_type=sparse) for arg in args]

    output_op_QAdd_QO = t -> sum( args_translated[i](t) for i=1:length(args_translated) )
    return output_op_QAdd_QO
end

function translate(op::QSym, b::QuantumOpticsBase.Basis; 
    parameter=Dict(), time_dep_param=Dict(), level_map=nothing, operators=Dict(), op_type=sparse)  
    # QSym are only fundamental symbolic operators, e.g. a, σ_-, x, ...

    time_dep_param = _normalize_time_dep_param(time_dep_param)
    if isempty(time_dep_param)
        return _translate_numeric(op, b; parameter=parameter, level_map=level_map, operators=operators, op_type=op_type) # TODO: test
    end

    output_op_QSym_QO_ = _translate_numeric_raw(op, b; level_map=level_map, operators=operators, op_type=op_type) # this is faster!
    output_op_QSym_QO = t -> output_op_QSym_QO_
    return output_op_QSym_QO
end
# # Default just use to_numeric() ? 
# function translate(op, b::QuantumOpticsBase.Basis; 
#     parameter=Dict(), time_dep_param=Dict(), level_map=nothing)  
    
#     if isempty(time_dep_param)
#         return to_numeric(substitute(op, parameter), b; level_map=level_map) # TODO: test
#     end

#     output_op_QSym_QO_ = to_numeric(op, b; level_map=level_map) 
#     output_op_QSym_QO = t -> output_op_QSym_QO_
#     return output_op_QSym_QO
# end
#

function translate(arg_c_, b::QuantumOpticsBase.Basis; 
    parameter=Dict(), time_dep_param=Dict(), level_map=nothing, operators=Dict(), op_type=sparse)  
    # should only be needed for numbers and symbolic paramters

    arg_c = substitute(arg_c_, parameter)
    one_b = _translate_one(b, operators, op_type)

    time_dep_param = _normalize_time_dep_param(time_dep_param)
    if isempty(time_dep_param)
        return arg_c * one_b
    end

    is_func, pref = _translate_prefactor(arg_c, time_dep_param)
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
    L = get_lindblad(G); H = get_hamiltonian(G)
    H_QO = translate(H,b; kwargs...)
    L_QO = [translate(L_,b; kwargs...) for L_ in L]
    return H_QO, L_QO
end
