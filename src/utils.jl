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
expect(op::QNumber, state; kwargs...) = numeric_average(op, state; kwargs...)
