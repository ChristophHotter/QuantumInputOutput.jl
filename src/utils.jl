# move to SQA
function numeric_average(op::SQA.QNumber, state::Vector; kwargs...)
    op_num = sparse(to_numeric(op, state[1]; kwargs...))
    # TODO: sparse, dense
    return QuantumOpticsBase.expect(op_num, state)
end
function numeric_average(avg::Average, state::Vector; kwargs...)
    op = undo_average(avg)
    return numeric_average(op, state; kwargs...)
end

expect(avg::Average, state; kwargs...) = numeric_average(avg, state; kwargs...)
expect(op::SQA.QNumber, state; kwargs...) = numeric_average(op, state; kwargs...)

function substitute_operators(x, dict)
    return substitute(x, dict)
end
function substitute_operators(x::SQA.QAdd, dict)
    return SQA.QAdd([substitute_operators(arg, dict) for arg in x.arguments])
end
function substitute_operators(x::SQA.QMul, dict)
    return substitute(x.arg_c, dict)*prod([substitute(arg, dict) for arg ∈ x.args_nc])
end