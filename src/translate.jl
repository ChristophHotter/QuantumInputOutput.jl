function translate(op::QuantumCumulants.QMul, b::QuantumOpticsBase.Basis; 
    parameter=Dict(), time_dep_param=Dict(), level_map=nothing, operators=Dict(), op_type=sparse)
    
    if isempty(time_dep_param)
        return op_type(to_numeric(substitute(op, parameter), b, operators; level_map=level_map)) # TODO: test
    end

    op_ = substitute(op, parameter) # First substitute all symbolic paramters with the numerical values -> only time-dependent parameters 

    if isa(op_, QuantumCumulants.QSym) # special case if parameter is 1.0: after substitute p*op -> 1.0*op -> op (has no fields arg_c, args_nc)
        output_op_QMul_QO_ = op_type(to_numeric(op_, b, operators; level_map=level_map)) # this is faster!
        output_op_QMul_QO = t -> output_op_QMul_QO_
        return output_op_QMul_QO
    elseif iszero(op_) # special case if parameter is 0.0
        return op_type(0*one(b))
    else 
        arg_c = op_.arg_c
        args_nc = op_.args_nc
        prod_args_nc = op_type(prod((to_numeric(arg, b, operators; level_map=level_map)) for arg in args_nc)) #TODO: dense, sparse, etc

        if !iscall(arg_c) # only one symbolic function in arg_c_ #TODO: time-dep parameter, e.g. g(t)
            arg_c_sub = substitute(arg_c, time_dep_param)
            if isa(arg_c_sub, Function)
                output_op_QMul_QO = t -> arg_c_sub(t)*prod_args_nc
                return output_op_QMul_QO 
            else
                prod_c_nc_num = arg_c_sub*prod_args_nc
                output_op_QMul_QO = t -> prod_c_nc_num
                return output_op_QMul_QO
            end
        else # more than one expression in arg_c (numeric, symbolic, funcion), e.g. 5*x    
            arg_c_args = arguments(arg_c)
            arg_c_op = operation(arg_c)

            # TODO: directly use gu(t) without substitution? (think about iscall)            
            if length(arg_c_args) == 1  # conj(g)*a, sqrt(g)*a # TODO: g^2*a [length(arg_c_args) == 2]
                if arg_c_op ∈ [conj, sqrt]
                    arg_c_1 = substitute(arg_c, time_dep_param)
                    output_QO_num = t -> arg_c_1(t)*prod_args_nc # only one element
                    return output_QO_num
                end
            else
                !(arg_c_op == *) && error("The operation of the prefactor is $(arg_c_op). 
                Currently, if you want to use a function like conj() or sqrt() on a time dependent 
                parameter, you need to include it in the time-dependent function dictionary, see e.g.: [REF-example]
                At the moment we are also limited to functions with only one argument.") # TODO
            end
            arg_c_all = [substitute(arg, time_dep_param) for arg ∈ arg_c_args]
            arg_c_numbers = filter(x->isa(x,Number), arg_c_all)
            prod_c_nc_num = prod(arg_c_numbers)*prod_args_nc
            arg_c_functions = filter(x->!isa(x,Number), arg_c_all)
            
            output_op_QMul_QO = t -> prod(arg_c_f(t) for arg_c_f in arg_c_functions)*prod_c_nc_num
            return output_op_QMul_QO
        end 
    end
end
#
function translate(op::QuantumCumulants.QAdd, b::QuantumOpticsBase.Basis; 
        parameter=Dict(), time_dep_param=Dict(), level_map=nothing, operators=Dict(), op_type=sparse)  

    op = expand(op) # should ensure that only multiplications are in arguments
    ### only static operators - returns a time-independent Hamiltonian and Lindblad
    if isempty(time_dep_param)
        return op_type(to_numeric(substitute(op, parameter), b, operators; level_map=level_map)) # TODO: test
    end

    ### static and time-dependent operators
    args = arguments(substitute(op, parameter)) 
    args_translated = [translate(arg, b; parameter=parameter, time_dep_param=time_dep_param, level_map=level_map) for arg in args]

    output_op_QAdd_QO = t -> sum( args_translated[i](t) for i=1:length(args_translated) )
    return output_op_QAdd_QO
end

function translate(op::QSym, b::QuantumOpticsBase.Basis; 
    parameter=Dict(), time_dep_param=Dict(), level_map=nothing, operators=Dict(), op_type=sparse)  
    # QSym are only fundamental symbolic operators, e.g. a, σ_-, x, ...

    if isempty(time_dep_param)
        return op_type(to_numeric(substitute(op, parameter), b, operators; level_map=level_map)) # TODO: test
    end

    output_op_QSym_QO_ = op_type(to_numeric(op, b, operators; level_map=level_map)) # this is faster!
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
    parameter=Dict(), time_dep_param=Dict(), level_map=nothing, operators=Dict())  
    # should only be needed for numbers and symbolic paramters

    arg_c = substitute(arg_c_, parameter)
    one_b = sparse(one(b))

    if isempty(time_dep_param)
        return arg_c*one_b
    end

    if !iscall(arg_c) # only one symbolic function in arg_c_ #TODO: time-dep parameter, e.g. g(t)
        arg_c_sub = substitute(arg_c, time_dep_param)
        if isa(arg_c_sub, Function)
            output_QO_num = t -> arg_c_sub(t)*one_b
            return output_QO_num 
        else
            prod_c_nc_num = arg_c_sub*one_b
            output_QO_num = t -> prod_c_nc_num
            return output_QO_num
        end
    else # more than one expression in arg_c (numeric, symbolic, funcion), e.g. 5*x. Or e.g. conj(g)*a
        arg_c_args = arguments(arg_c)
        arg_c_op = operation(arg_c)

        # TODO: directly use gu(t) without substitution? (think about iscall)
        if length(arg_c_args) == 1  # e.g. conj(g), sqrt(g) # TODO: g^2 [length(arg_c_args) == 2]
            if arg_c_op ∈ [conj, sqrt]
                arg_c_1 = substitute(arg_c, time_dep_param)
                output_QO_num = t -> arg_c_1(t)*one_b # only one element
                return output_QO_num
            end
        else
            !(arg_c_op == *) && error("The operation of the prefactor is $(arg_c_op). 
            Currently, if you want to use a function like conj() or sqrt() on a time dependent 
            parameter, you need to include it in the time-dependent function dictionary, see e.g.: [REF-example]
            At the moment we are also limited to functions with only one argument.") # TODO
        end
        arg_c_all = [substitute(arg, time_dep_param) for arg ∈ arg_c_args]
        arg_c_numbers = filter(x->isa(x,Number), arg_c_all)
        arg_c_functions = filter(x->!isa(x,Number), arg_c_all)
        prod_c_nc_num = prod(arg_c_numbers)*one_b 

        output_QO_num = t -> prod(arg_c_f(t) for arg_c_f in arg_c_functions)*prod_c_nc_num
        return output_QO_num
    end 
end

function translate(ops::Vector, b::QuantumOpticsBase.Basis; kwargs...)
    return [translate(op, b; kwargs...) for op in ops]
end
function translate(G::SLH, b::QuantumOpticsBase.Basis; 
        parameter=Dict(), time_dep_param=Dict(), level_map=nothing, operators=Dict()) 
    L = get_lindblad(G); H = get_hamiltonian(G)
    H_QO = translate(H,b; parameter=parameter,time_dep_param=Dict(),time=nothing,level_map=level_map,operators=operators)
    L_QO = [translate(L_,b; parameter=parameter,time_dep_param=Dict(),time=nothing,level_map=level_map,operators=operators) for L_ in L]
    return H_QO, L_QO
end


