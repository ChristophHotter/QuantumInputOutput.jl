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

"""
    substitute_operators(op, dict::Dict)

Like `substitute(op, dict::Dict)` but with special handling for `QMul` and `QAdd`.
This is needed if an operator is substitute by a `QMul` or `QAdd`, e.g. 
    ``a_1 -> g_2*a_2 + g_3*a_3``
"""
function substitute_operators(op, dict::Dict)
    return substitute(op, dict::Dict)
end
function substitute_operators(op::SQA.QAdd, dict::Dict)
    return SQA.QAdd([substitute_operators(arg, dict::Dict) for arg in op.arguments])
end
function substitute_operators(op::SQA.QMul, dict::Dict)
    return substitute(op.arg_c, dict::Dict)*prod([substitute(arg, dict::Dict) for arg ∈ op.args_nc])
end
